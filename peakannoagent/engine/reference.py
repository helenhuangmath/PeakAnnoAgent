"""Gene-model references for peak annotation.

A reference is a set of sorted numpy arrays (transcription start sites, exons, introns, UTRs,
transcript ends and genes) built once from a GTF file and cached. Later runs read the cache,
which is what makes annotation fast. Coordinates inside a reference are one-based and closed,
as in GTF.
"""
from __future__ import annotations

import gzip
import pickle
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from .common import cache_dir, download

FORMAT = 2
GENOMES = ("mm10", "mm39", "hg38", "hg19")
GENOME_ALIASES = {"mouse": "mm10", "human": "hg38", "grcm38": "mm10", "grch38": "hg38", "grcm39": "mm39",
                  "grch37": "hg19"}
GENCODE_RELEASE = {"hg38": "v46", "hg19": "v46lift37", "mm39": "vM35", "mm10": "vM25"}
REFERENCE_OPTIONS = {
    "gencode": "GENCODE basic annotation, pinned release per build (default); "
               "'gencode:v44' or 'gencode:vM25' picks a release",
    "refseq": "NCBI RefSeq (UCSC ncbiRefSeq GTF)",
    "ucsc": "UCSC knownGene (UCSC knownGene GTF)",
    "<file or URL>.gtf[.gz]": "any GTF with gene_id and transcript_id attributes (Ensembl, GENCODE, custom)",
}


def normalize_genome(genome: str) -> str:
    g = GENOME_ALIASES.get(str(genome).lower(), str(genome).lower())
    if g not in GENOMES:
        raise ValueError(f"Unsupported genome '{genome}'. Choose one of {', '.join(GENOMES)}.")
    return g


def gencode_url(release: str) -> str:
    release = re.sub(r"^gencode[:._-]?", "", release, flags=re.I)
    if not release.startswith("v"):
        release = "v" + release
    mouse = release.startswith("vM")
    number = re.sub(r"lift37$", "", release[1:])
    folder = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{'mouse' if mouse else 'human'}/release_{number}/"
    if release.endswith("lift37"):
        folder += "GRCh37_mapping/"
    return f"{folder}gencode.{release}.basic.annotation.gtf.gz"


def resolve_reference(genome: str = "mm10", reference: Optional[str] = None, gtf: Optional[str] = None) -> dict:
    """Turn (genome, reference, gtf) into a source: a label and a GTF path or URL."""
    genome = normalize_genome(genome)
    if gtf:
        reference = gtf
    reference = reference or "gencode"
    key = reference.lower()
    ucsc = f"https://hgdownload.soe.ucsc.edu/goldenPath/{genome}/bigZips/genes/{genome}"
    if key == "gencode" or re.match(r"gencode[:._-]?v", key):
        release = GENCODE_RELEASE[genome] if key == "gencode" else re.sub(r"^gencode[:._-]?", "", reference, flags=re.I)
        return {"genome": genome, "source": f"GENCODE {release} basic", "label": f"gencode.{release}.basic",
                "url": gencode_url(release)}
    if key in ("refseq", "ncbirefseq", "refgene"):
        return {"genome": genome, "source": "NCBI RefSeq (UCSC ncbiRefSeq)", "label": f"{genome}.ncbiRefSeq",
                "url": f"{ucsc}.ncbiRefSeq.gtf.gz"}
    if key in ("ucsc", "knowngene", "ucsc_knowngene"):
        return {"genome": genome, "source": "UCSC knownGene", "label": f"{genome}.knownGene",
                "url": f"{ucsc}.knownGene.gtf.gz"}
    if re.search(r"\.(gtf|gff2?)(\.gz)?$", key):
        name = reference.rsplit("/", 1)[-1]
        out = {"genome": genome, "source": f"GTF {name}", "label": re.sub(r"\.(gtf|gff2?)(\.gz)?$", "", name, flags=re.I)}
        if re.match(r"(https?|ftp)://", reference):
            out["url"] = reference
        else:
            path = Path(reference).expanduser()
            if not path.is_file():
                raise FileNotFoundError(f"GTF file not found: {path}")
            out["path"] = str(path.resolve())
        return out
    raise ValueError(f"Unknown reference '{reference}'. Use gencode, gencode:<release>, refseq, ucsc, "
                     f"or a GTF file or URL.")


@dataclass
class Reference:
    """Gene models as arrays. Interval tables are sorted by (chromosome, start)."""
    meta: dict
    chroms: List[str]                      # chromosome names; arrays hold indexes into this list
    genes: Dict[str, list]                 # gene_id, symbol, gene_type, chrom, start, end, strand (lists)
    tss: Dict[str, np.ndarray]             # chrom, pos, strand (+1/-1), gene (index), pc; tx_id is a list
    tes: Dict[str, np.ndarray]             # chrom, pos, strand, pc
    exons: Dict[str, np.ndarray]           # chrom, start, end, first, pc, first_pc
    introns: Dict[str, np.ndarray]
    utr5: Dict[str, np.ndarray]            # chrom, start, end, pc
    utr3: Dict[str, np.ndarray]
    chrom_ends: Dict[str, int] = field(default_factory=dict)

    def __repr__(self) -> str:
        m = self.meta
        return (f"<Reference {m['genome']}, {m['source']}: {m['n_genes']:,} genes "
                f"({m['n_protein_coding']:,} protein-coding), {m['n_transcripts']:,} transcripts>")

    def gene_types(self, gene_types: str = "auto") -> str:
        """"auto" means protein-coding genes: GTF sources list many pseudogenes and lncRNAs,
        which would otherwise become the nearest gene of many peaks."""
        if gene_types not in ("auto", "protein_coding", "all"):
            raise ValueError("gene_types must be 'auto', 'protein_coding' or 'all'.")
        return "protein_coding" if gene_types == "auto" else gene_types

    def symbols(self, gene_types: str = "auto") -> List[str]:
        keep = self.gene_types(gene_types) == "all"
        return sorted({s for s, t in zip(self.genes["symbol"], self.genes["gene_type"])
                       if s and (keep or t == "protein_coding")})


_ATTR = {key: re.compile(r'(?<![A-Za-z_])' + key + r' "([^"]*)"')
         for key in ("gene_id", "transcript_id", "gene_name", "gene_type", "gene_biotype")}


def _attr(text: str, key: str) -> Optional[str]:
    m = _ATTR[key].search(text)
    return m.group(1) if m else None


def _ucsc_name(chrom: str, has_prefix: bool) -> str:
    """"1", "X", "MT" (Ensembl) -> "chr1", "chrX", "chrM"; names with "chr" are kept."""
    if has_prefix:
        return chrom
    if re.fullmatch(r"[0-9]+|[XYxy]|MT|M", chrom):
        return "chr" + ("M" if chrom == "MT" else chrom)
    return chrom


def read_gtf_models(path: str) -> dict:
    """Transcript models of a GTF: exons and coding range per transcript, gene per transcript,
    and gene names and types."""
    opener = gzip.open if str(path).endswith(".gz") else open
    exons: Dict[str, list] = {}            # transcript -> [chrom, strand, [(start, end), ...]]
    cds: Dict[str, list] = {}              # transcript -> [min start, max end]
    tx_gene: Dict[str, str] = {}
    genes: Dict[str, list] = {}            # gene -> [symbol, type]
    has_prefix = None
    with opener(path, "rt") as handle:
        for line in handle:
            if line[0] == "#":
                continue
            parts = line.split("\t", 8)
            if len(parts) < 9:
                continue
            kind = parts[2]
            if kind not in ("exon", "CDS", "gene", "transcript"):
                continue
            attr = parts[8]
            if kind == "gene":
                gid = _attr(attr, "gene_id")
                if gid:
                    genes[gid] = [_attr(attr, "gene_name"), _attr(attr, "gene_type") or _attr(attr, "gene_biotype")]
                continue
            tx = _attr(attr, "transcript_id")
            if not tx:
                if kind == "exon":
                    raise ValueError(f"{path} is not a GTF with transcript_id attributes (GFF3 is not supported).")
                continue
            start, end = int(parts[3]), int(parts[4])
            if kind == "exon":
                rec = exons.get(tx)
                if rec is None:
                    if has_prefix is None:
                        has_prefix = parts[0].startswith("chr")
                    rec = exons[tx] = [_ucsc_name(parts[0], has_prefix), parts[6], []]
                rec[2].append((start, end))
            elif kind == "CDS":
                rng = cds.get(tx)
                if rng is None:
                    cds[tx] = [start, end]
                else:
                    if start < rng[0]:
                        rng[0] = start
                    if end > rng[1]:
                        rng[1] = end
            if tx not in tx_gene:
                gid = _attr(attr, "gene_id")
                if gid:
                    tx_gene[tx] = gid
                    if gid not in genes:
                        genes[gid] = [_attr(attr, "gene_name"), _attr(attr, "gene_type") or _attr(attr, "gene_biotype")]
                    elif genes[gid][0] is None:
                        genes[gid][0] = _attr(attr, "gene_name")
    if not exons:
        raise ValueError(f"{path} has no exon records; a GTF with exons is required.")
    # genePred-derived GTFs (UCSC knownGene) repeat the transcript ID as gene_id: group by name.
    if len(set(tx_gene.values())) == len(tx_gene) and len(tx_gene) > 1:
        named = {tx: genes[g][0] for tx, g in tx_gene.items() if genes[g][0]}
        if len(named) > 0.9 * len(tx_gene) and len(set(named.values())) < len(named):
            genes = {name: [name, None] for name in set(named.values())}
            tx_gene = named
    return {"exons": exons, "cds": cds, "tx_gene": tx_gene, "genes": genes}


def _table(rows: list, columns: tuple, dtypes: tuple) -> Dict[str, np.ndarray]:
    """Rows -> column arrays sorted by the first three columns, with duplicate intervals merged."""
    if not rows:
        return {c: np.zeros(0, dtype=d) for c, d in zip(columns, dtypes)}
    arr = [np.array([r[i] for r in rows], dtype=d) for i, d in enumerate(dtypes)]
    order = np.lexsort((arr[2], arr[1], arr[0]))
    arr = [a[order] for a in arr]
    key = np.ones(len(order), dtype=bool)
    key[1:] = (arr[0][1:] != arr[0][:-1]) | (arr[1][1:] != arr[1][:-1]) | (arr[2][1:] != arr[2][:-1])
    group = np.cumsum(key) - 1
    out = {columns[i]: arr[i][key] for i in range(3)}
    for i in range(3, len(columns)):                      # flags: any() within duplicates
        out[columns[i]] = np.bincount(group, weights=arr[i].astype(float), minlength=group[-1] + 1) > 0
    return out


def assemble_reference(models: dict, meta: dict) -> Reference:
    tx_gene, cds, gene_info = models["tx_gene"], models["cds"], models["genes"]
    exons = {tx: rec for tx, rec in models["exons"].items() if tx in tx_gene}
    if not exons:
        raise ValueError("The reference contains no transcripts assigned to a gene.")
    coding_genes = {tx_gene[tx] for tx in cds if tx in exons}
    gene_ids = sorted({tx_gene[tx] for tx in exons})
    gene_index = {g: i for i, g in enumerate(gene_ids)}
    gene_type = []
    for g in gene_ids:
        kind = gene_info.get(g, [None, None])[1]
        gene_type.append(kind or ("protein_coding" if g in coding_genes else "noncoding"))
    is_pc = [t == "protein_coding" for t in gene_type]
    chroms = sorted({rec[0] for rec in exons.values()})
    chrom_index = {c: i for i, c in enumerate(chroms)}

    span: Dict[str, list] = {}
    tss_rows, tes_rows, exon_rows, intron_rows, utr5_rows, utr3_rows = {}, [], [], [], [], []
    chrom_ends: Dict[str, int] = {}
    for tx, (chrom, strand, parts) in exons.items():
        gene = tx_gene[tx]
        gi, ci = gene_index[gene], chrom_index[chrom]
        pc = is_pc[gi]
        parts.sort()
        start, end = parts[0][0], max(p[1] for p in parts)
        minus = strand == "-"
        sign = -1 if minus else 1
        g = span.get(gene)
        if g is None:
            span[gene] = [chrom, strand, start, end]
        else:
            g[2], g[3] = min(g[2], start), max(g[3], end)
        chrom_ends[chrom] = max(chrom_ends.get(chrom, 0), end)
        tss_pos = end if minus else start
        tss_rows.setdefault((ci, tss_pos, sign, gi), [tx, pc])       # one row per distinct (gene, TSS)
        tes_rows.append((ci, start if minus else end, sign, pc))
        n = len(parts)
        rng = cds.get(tx)
        for k, (s, e) in enumerate(parts):
            rank = n - k if minus else k + 1
            first = rank == 1
            exon_rows.append((ci, s, e, first, pc, first and pc))
            if k:
                i_s, i_e = parts[k - 1][1] + 1, s - 1
                if i_e >= i_s:
                    i_first = (n - k if minus else k) == 1
                    intron_rows.append((ci, i_s, i_e, i_first, pc, i_first and pc))
            if rng:                                               # exonic parts outside the coding range
                if s < rng[0]:
                    (utr3_rows if minus else utr5_rows).append((ci, s, min(e, rng[0] - 1), pc))
                if e > rng[1]:
                    (utr5_rows if minus else utr3_rows).append((ci, max(s, rng[1] + 1), e, pc))

    ranked = (("chrom", "start", "end", "first", "pc", "first_pc"), (np.int32, np.int64, np.int64, bool, bool, bool))
    plain = (("chrom", "start", "end", "pc"), (np.int32, np.int64, np.int64, bool))
    keys = sorted(tss_rows)
    tss = {"chrom": np.array([k[0] for k in keys], dtype=np.int32), "pos": np.array([k[1] for k in keys], dtype=np.int64),
           "strand": np.array([k[2] for k in keys], dtype=np.int8), "gene": np.array([k[3] for k in keys], dtype=np.int32),
           "pc": np.array([tss_rows[k][1] for k in keys], dtype=bool), "tx_id": [tss_rows[k][0] for k in keys]}
    tes_set = sorted({r[:3] for r in tes_rows})
    tes_pc = {r[:3] for r in tes_rows if r[3]}
    tes = {"chrom": np.array([r[0] for r in tes_set], dtype=np.int32), "pos": np.array([r[1] for r in tes_set], dtype=np.int64),
           "strand": np.array([r[2] for r in tes_set], dtype=np.int8), "pc": np.array([r in tes_pc for r in tes_set], dtype=bool)}
    symbols = [gene_info.get(g, [None, None])[0] or None for g in gene_ids]
    if not any(symbols):
        symbols = list(gene_ids)                                  # GTF without gene_name
    genes = {"gene_id": gene_ids, "symbol": symbols, "gene_type": gene_type,
             "chrom": [span[g][0] for g in gene_ids], "strand": [span[g][1] for g in gene_ids],
             "start": [span[g][2] for g in gene_ids], "end": [span[g][3] for g in gene_ids]}
    meta = dict(meta, n_genes=len(gene_ids), n_transcripts=len(exons), n_protein_coding=sum(is_pc),
                format=FORMAT, built=time.strftime("%Y-%m-%d %H:%M:%S"))
    return Reference(meta=meta, chroms=chroms, genes=genes, tss=tss, tes=tes,
                     exons=_table(exon_rows, *ranked), introns=_table(intron_rows, *ranked),
                     utr5=_table(utr5_rows, *plain), utr3=_table(utr3_rows, *plain), chrom_ends=chrom_ends)


def build_reference(genome: str = "mm10", reference: Optional[str] = None, gtf: Optional[str] = None,
                    cache: Optional[str] = None, rebuild: bool = False) -> Reference:
    """Build the annotation reference for a genome, or read it from the cache."""
    spec = resolve_reference(genome, reference, gtf)
    root = cache_dir(cache)
    path = spec.get("path") or str(download(spec["url"], root / "gtf" / spec["url"].rsplit("/", 1)[-1],
                                            f"{spec['source']} annotation"))
    stat = Path(path).stat()
    name = re.sub(r"[^A-Za-z0-9._-]+", "-", f"{spec['genome']}__{spec['label']}__{stat.st_size}-{int(stat.st_mtime)}")
    cached = root / f"ref{FORMAT}__{name}.pkl"
    if not rebuild and cached.is_file():
        try:
            with open(cached, "rb") as handle:
                ref = pickle.load(handle)
            if isinstance(ref, Reference) and ref.meta.get("format") == FORMAT:
                ref.meta.update(cache=str(cached), from_cache=True)
                return ref
        except Exception:
            pass                                                  # unreadable cache: rebuild
    meta = {"genome": spec["genome"], "source": spec["source"], "label": spec["label"], "file": path,
            "url": spec.get("url")}
    ref = assemble_reference(read_gtf_models(path), meta)
    tmp = cached.with_name(cached.name + f".{time.time_ns()}.part")
    with open(tmp, "wb") as handle:
        pickle.dump(ref, handle, protocol=pickle.HIGHEST_PROTOCOL)
    tmp.replace(cached)
    ref.meta.update(cache=str(cached), from_cache=False)
    return ref
