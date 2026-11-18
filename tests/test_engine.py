"""The analysis engine: gene-model references, fast annotation, enrichment and networks.
Everything here runs offline, on small GTF files written by the test or shipped in tests/data."""

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from peakannoagent import engine
from peakannoagent.cli import main
from peakannoagent.engine import Peaks, annotate_peaks, build_reference, read_peak_file
from peakannoagent.engine import common, enrich, networks, reference
from peakannoagent.engine.annotate import COLUMNS

from helpers import DATA

GTF = DATA / "gencode_vM25_chr11_subset.gtf.gz"
BED = DATA / "example_peaks_mm10_chr11.bed"


def toy_gtf(tmp) -> str:
    """Ensembl-style names ("1"); gene A on the plus strand, B on the minus strand, C a lncRNA."""
    def row(kind, start, end, strand, attrs):
        return "\t".join(["1", "toy", kind, str(start), str(end), ".", strand, ".", attrs])
    a = 'gene_id "GA"; transcript_id "TA1"; gene_name "Alpha"; gene_type "protein_coding";'
    b = 'gene_id "GB"; transcript_id "TB1"; gene_name "Beta"; gene_type "protein_coding";'
    c = 'gene_id "GC"; transcript_id "TC1"; gene_name "Gamma"; gene_type "lncRNA";'
    lines = ["#!toy annotation",
             row("gene", 1000, 3500, "+", 'gene_id "GA"; gene_name "Alpha"; gene_type "protein_coding";'),
             row("transcript", 1000, 3500, "+", a),
             row("exon", 1000, 1200, "+", a), row("exon", 2000, 2300, "+", a), row("exon", 3000, 3500, "+", a),
             row("CDS", 1100, 1200, "+", a), row("CDS", 2000, 2300, "+", a), row("CDS", 3000, 3200, "+", a),
             row("gene", 10000, 11500, "-", 'gene_id "GB"; gene_name "Beta"; gene_type "protein_coding";'),
             row("transcript", 10000, 11500, "-", b),
             row("exon", 10000, 10400, "-", b), row("exon", 11000, 11500, "-", b),
             row("CDS", 10200, 10400, "-", b), row("CDS", 11000, 11300, "-", b),
             row("gene", 20000, 20500, "+", 'gene_id "GC"; gene_name "Gamma"; gene_type "lncRNA";'),
             row("transcript", 20000, 20500, "+", c), row("exon", 20000, 20500, "+", c)]
    path = Path(tmp) / "toy.gtf"
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def toy_peaks() -> Peaks:
    """BED (zero-based start) intervals chosen to hit one feature each."""
    start = [849, 1059, 2099, 1499, 2499, 3299, 3599, 5999, 11599, 10049, 20099, 99]
    end = [900, 1090, 2200, 1600, 2600, 3400, 3700, 6100, 11650, 10100, 20200, 200]
    chrom = ["chr1"] * 11 + ["chr2"]
    return Peaks([f"p{i + 1:02d}" for i in range(12)], np.array(chrom, dtype=object), np.array(start), np.array(end))


def by_peak(anno: dict, column: str) -> dict:
    return dict(zip(anno["PeakID"], anno[column]))


class ReferenceTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_reference_names_resolve_to_sources(self):
        g = reference.resolve_reference("mm10")
        self.assertTrue(g["url"].endswith("Gencode_mouse/release_M25/gencode.vM25.basic.annotation.gtf.gz"))
        self.assertIn("Gencode_human/release_44/gencode.v44.basic", reference.resolve_reference("human", "gencode:v44")["url"])
        self.assertIn("release_46/GRCh37_mapping/gencode.v46lift37", reference.resolve_reference("hg19", "gencode")["url"])
        self.assertTrue(reference.resolve_reference("hg38", "refseq")["url"].endswith("hg38.ncbiRefSeq.gtf.gz"))
        self.assertTrue(reference.resolve_reference("mm39", "ucsc")["url"].endswith("mm39.knownGene.gtf.gz"))
        gtf = toy_gtf(self.tmp)
        self.assertEqual(reference.resolve_reference("mm10", gtf=gtf)["path"], str(Path(gtf).resolve()))
        with self.assertRaises(ValueError):
            reference.resolve_reference("mm10", "nonsense")
        with self.assertRaises(ValueError):
            reference.resolve_reference("rn6")
        with self.assertRaises(FileNotFoundError):
            reference.resolve_reference("mm10", gtf=str(Path(self.tmp) / "absent.gtf"))

    def test_gtf_becomes_a_cached_reference(self):
        gtf = toy_gtf(self.tmp)
        ref = build_reference("mm10", gtf=gtf, cache=self.tmp)
        self.assertFalse(ref.meta["from_cache"])
        self.assertEqual((ref.meta["n_genes"], ref.meta["n_protein_coding"], ref.meta["n_transcripts"]), (3, 2, 3))
        self.assertEqual(ref.chroms, ["chr1"])                       # "1" converted to UCSC style
        self.assertEqual(sorted(ref.genes["symbol"]), ["Alpha", "Beta", "Gamma"])
        beta = ref.genes["gene_id"].index("GB")
        self.assertEqual(int(ref.tss["pos"][list(ref.tss["gene"]).index(beta)]), 11500)   # minus strand: TSS is the end
        self.assertEqual(list(zip(ref.utr5["start"], ref.utr5["end"])), [(1000, 1099), (11301, 11500)])
        self.assertEqual(list(zip(ref.utr3["start"], ref.utr3["end"])), [(3201, 3500), (10000, 10199)])
        introns = list(zip(ref.introns["start"], ref.introns["end"], ref.introns["first"]))
        self.assertEqual(introns, [(1201, 1999, True), (2301, 2999, False), (10401, 10999, True)])
        self.assertIn("3 genes (2 protein-coding)", repr(ref))
        again = build_reference("mm10", gtf=gtf, cache=self.tmp)
        self.assertTrue(again.meta["from_cache"])
        np.testing.assert_array_equal(again.tss["pos"], ref.tss["pos"])

    def test_gff3_and_exonless_files_are_rejected(self):
        bad = Path(self.tmp) / "bad.gtf"
        bad.write_text("1\tx\texon\t1\t10\t.\t+\t.\tID=exon1;Parent=tx1\n")
        with self.assertRaisesRegex(ValueError, "transcript_id"):
            build_reference("mm10", gtf=str(bad), cache=self.tmp)
        empty = Path(self.tmp) / "empty.gtf"
        empty.write_text("#nothing\n")
        with self.assertRaisesRegex(ValueError, "no exon records"):
            build_reference("mm10", gtf=str(empty), cache=self.tmp)

    def test_offline_never_downloads(self):
        common.OFFLINE = True
        try:
            with self.assertRaisesRegex(RuntimeError, "offline"):
                build_reference("mm10", "gencode", cache=self.tmp)
        finally:
            common.OFFLINE = False


class AnnotationTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.ref = build_reference("mm10", gtf=toy_gtf(self.tmp), cache=self.tmp)

    def test_feature_gene_and_signed_distance(self):
        anno = annotate_peaks(toy_peaks(), self.ref, promoter_upstream=200, promoter_downstream=50)
        self.assertEqual(set(anno), set(COLUMNS))
        feature = by_peak(anno, "feature")
        self.assertEqual([feature[f"p{i:02d}"] for i in range(1, 13)],
                         ["Promoter", "5' UTR", "Other Exon", "1st Intron", "Other Intron", "3' UTR",
                          "Downstream (<=300bp)", "Distal Intergenic", "Promoter", "3' UTR", "Distal Intergenic",
                          "Distal Intergenic"])
        distance, symbol = by_peak(anno, "distance_to_tss"), by_peak(anno, "nearest_symbol")
        self.assertEqual(distance["p01"], -125)                      # upstream of a plus-strand gene
        self.assertEqual(distance["p09"], -125)                      # upstream of a minus-strand gene
        self.assertEqual(distance["p08"], 5050)
        self.assertEqual([symbol[p] for p in ("p01", "p08", "p09")], ["Alpha", "Alpha", "Beta"])
        self.assertEqual(symbol["p11"], "Beta")                      # the lncRNA is not used by default
        self.assertIsNone(symbol["p12"])                             # chromosome absent from the reference
        promoter = by_peak(anno, "promoter_symbols")
        self.assertEqual((promoter["p01"], promoter["p09"], promoter["p08"]), ("Alpha", "Beta", None))
        region = by_peak(anno, "region_class")
        self.assertEqual([region[p] for p in ("p01", "p08", "p12")], ["promoter", "distal", "distal_far"])
        self.assertEqual(by_peak(anno, "feature_overlaps")["p06"], "3' UTR;Exon")
        self.assertEqual((by_peak(anno, "Start")["p01"], by_peak(anno, "End")["p01"]), (849, 900))   # BED, as given

    def test_gene_types_and_promoter_bins(self):
        everything = annotate_peaks(toy_peaks(), self.ref, 200, 50, gene_types="all")
        self.assertEqual(by_peak(everything, "feature")["p11"], "1st Exon")
        self.assertEqual(by_peak(everything, "nearest_gene_type")["p11"], "lncRNA")
        wide = annotate_peaks(toy_peaks(), self.ref, 3000, 3000)
        labels = {f for f in wide["feature"] if f.startswith("Promoter")}
        self.assertTrue(labels and labels <= {"Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)"})
        with self.assertRaises(ValueError):
            annotate_peaks(toy_peaks(), self.ref, gene_types="some")

    def test_invalid_input_is_rejected(self):
        bad = Peaks(["a"], np.array(["chr1"], dtype=object), np.array([10]), np.array([5]))
        with self.assertRaisesRegex(ValueError, "half-open"):
            annotate_peaks(bad, self.ref)
        with self.assertRaises(TypeError):
            annotate_peaks(toy_peaks(), {})

    def test_peak_files(self):
        path = Path(self.tmp) / "p.bed"
        path.write_text("chr1\t10\t20\nchr1\t30\t40\n")
        self.assertEqual(read_peak_file(path).ids, ["peak_00001", "peak_00002"])
        path.write_text("Chr\tStart\tEnd\tname\nchr1\t10\t20\ta\nchr2\t5\t9\tb\n")
        peaks = read_peak_file(path)
        self.assertEqual((peaks.ids, list(peaks.start)), (["a", "b"], [10, 5]))
        path.write_text("chr1\t1.5\t3\n")
        with self.assertRaisesRegex(ValueError, "whole numbers"):
            read_peak_file(path)
        self.assertEqual(read_peak_file(BED).ids[0][:8], "tex_peak")

    def test_run_annotation_writes_tables_gene_sets_and_all_figures(self):
        out = Path(self.tmp) / "anno"
        res = engine.run_annotation(read_peak_file(BED), out, "t", genome="mm10", gtf=str(GTF), cache=self.tmp)
        json.dumps(res)                                              # the summary must be JSON-serialisable
        self.assertEqual(res["n_peaks"], 12)
        self.assertEqual(res["reference"]["gene_types"], "protein_coding")
        self.assertAlmostEqual(sum(res["feature_pct"].values()), 100, delta=0.1)
        self.assertAlmostEqual(sum(res["region_distribution_pct"].values()), 100, delta=0.1)
        self.assertEqual(set(res["figures"]), {"feature_pie", "feature_bar", "tss_distance_bins", "feature_upset",
                                               "region_distribution", "tss_profile", "tss_distance",
                                               "chromosome_coverage"})
        for name, pdf in res["figures"].items():
            self.assertEqual(Path(pdf), out / f"t_{name}.pdf")
            self.assertTrue(Path(pdf).is_file(), name)
        self.assertFalse((out / ".preview").exists())                # annotation alone writes PDFs only
        for table in res["tables"].values():
            self.assertTrue(Path(table).is_file())
        sets = common.read_tsv(res["gene_sets"])
        universe = {r["SYMBOL"] for r in common.read_tsv(res["universe"])}
        self.assertTrue(sets and all(r["SYMBOL"] in universe for r in sets))
        self.assertEqual({r["Cluster"] for r in sets}, {"promoter", "distal", "all"})
        self.assertEqual(res["warnings"], [])

    def test_annotate_command(self):
        out = Path(self.tmp) / "cli"
        code = main(["annotate", "--peaks", str(BED), "--genome", "mm10", "--out", str(out), "--prefix", "t",
                     "--gtf", str(GTF), "--no-figures", "--no-download"])
        self.assertEqual(code, 0)
        summary = json.loads((out / "t_annotation.json").read_text())
        self.assertEqual(summary["n_peaks"], 12)
        self.assertEqual(summary["figures"], {})
        self.assertTrue((out / "t_peak_genes.tsv").is_file())
        self.assertEqual(main(["annotate", "--peaks", str(BED), "--genome", "rn6", "--out", str(out)]), 2)


class PaletteAndPreviewTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(common.set_palette, "chippeakanno")
        self.addCleanup(common.set_previews, False)

    def test_chippeakanno_is_the_default_and_other_palettes_can_be_chosen(self):
        from peakannoagent.engine.figures import feature_colors
        self.assertEqual(common.palette().name, "chippeakanno")
        self.assertEqual(common.palette().colors(2), ["#D55E00", "#0072B2"])       # ChIPpeakAnno labelColors
        features = ["Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "5' UTR", "3' UTR", "1st Exon",
                    "Other Exon", "1st Intron", "Other Intron", "Downstream (<=300bp)", "Distal Intergenic"]
        cols = feature_colors(features)
        self.assertEqual(cols["Promoter (<=1kb)"], "#FF8E32")                      # promoter-level oranges
        self.assertEqual(cols["Promoter (2-3kb)"], "#FFCA99")
        self.assertEqual((cols["5' UTR"], cols["Downstream (<=300bp)"]), ("#0072B2", "#E69F00"))
        self.assertEqual(len(set(cols.values())), len(features))
        self.assertEqual(feature_colors(["Promoter"])["Promoter"], "#D55E00")
        self.assertEqual(common.set_palette("NPG"), "npg")
        self.assertEqual(common.palette().colors(2), ["#E64B35", "#4DBBD5"])       # Nature Publishing Group colours
        cols = feature_colors(["Promoter (<=1kb)", "Promoter (1-2kb)", "3' UTR", "Distal Intergenic"])
        self.assertEqual(cols["Promoter (<=1kb)"].upper(), "#E64B35")
        self.assertEqual(len(set(cols.values())), 4)
        self.assertEqual(cols["3' UTR"], feature_colors(["3' UTR"])["3' UTR"])     # stable from run to run
        self.assertEqual(common.set_palette("Spectral"), "spectral")
        spectral = common.palette().colors(6)
        self.assertEqual((spectral[0], spectral[-1]), ("#9e0142", "#5e4fa2"))      # the ends of ColorBrewer Spectral
        self.assertEqual(len(common.palette().ordered(6)), 6)
        self.assertEqual(common.set_palette("set2"), "Set2")                       # Matplotlib names, any case
        self.assertEqual(len(set(common.palette().colors(8))), 8)
        self.assertEqual(common.set_palette("viridis"), "viridis")
        self.assertEqual(len(set(common.palette().colors(11))), 11)
        with self.assertRaisesRegex(ValueError, "Unknown palette 'rainbows'"):
            common.set_palette("rainbows")

    def test_every_palette_option_draws_the_annotation_figures(self):
        for name in common.PALETTE_OPTIONS:
            out = Path(self.tmp) / name
            res = engine.run_annotation(read_peak_file(BED), out, "t", genome="mm10", gtf=str(GTF), cache=self.tmp,
                                        palette=name)
            self.assertEqual(len(res["figures"]), 8, name)
            self.assertEqual(res["warnings"], [], name)

    def test_figures_are_pdf_only_unless_previews_are_on(self):
        out = Path(self.tmp) / "plain"
        fig = common.save_figure(lambda f: f.add_axes([0.1, 0.1, 0.8, 0.8]).plot([0, 1]), out, "a")
        self.assertEqual(Path(fig), out / "a.pdf")
        self.assertFalse((out / ".preview").exists())
        common.set_previews(True)
        common.save_figure(lambda f: f.add_axes([0.1, 0.1, 0.8, 0.8]).plot([0, 1]), out, "b")
        self.assertTrue(common.preview_path(out / "b.pdf").is_file())

    def test_unknown_palette_is_rejected_by_the_command(self):
        code = main(["annotate", "--peaks", str(BED), "--genome", "mm10", "--out", self.tmp, "--gtf", str(GTF),
                     "--palette", "rainbows", "--no-download"])
        self.assertEqual(code, 2)
        self.assertFalse(list(Path(self.tmp).glob("*.pdf")))


class EnrichmentTest(unittest.TestCase):
    def test_hypergeometric_test_and_bh(self):
        import math
        self.assertAlmostEqual(enrich.hypergeom_sf(3, 1000, 5, 4), sum(
            math.comb(5, k) * math.comb(995, 4 - k) for k in (3, 4)) / math.comb(1000, 4))
        self.assertEqual(enrich.hypergeom_sf(0, 100, 5, 4), 1.0)
        self.assertEqual(common.bh_adjust([0.01, 0.04, 0.03]), [0.03, 0.04, 0.04])

    def test_ora_backgrounds_and_case(self):
        universe = [f"G{i}" for i in range(1, 1001)]
        sets = {"A": ["g1", "G2", "G3", "G4", "G5"], "B": ["G900", "G901", "G902"]}
        rows = enrich.ora({"x": ["G1", "G2", "G3", "G10"]}, sets, universe, "custom", min_size=3,
                          annotated_background=False)
        a = next(r for r in rows if r["term"] == "A")
        self.assertEqual((a["count"], a["bg_ratio"], a["gene_ratio"], a["genes"]), (3, "5/1000", "3/4", "G1;G2;G3"))
        self.assertAlmostEqual(a["pvalue"], enrich.hypergeom_sf(3, 1000, 5, 4))
        self.assertEqual(next(r for r in rows if r["term"] == "B")["count"], 0)
        annotated = enrich.ora({"x": ["G1", "G2", "G3", "G10"]}, sets, universe, "m", min_size=3)
        self.assertEqual(annotated[0]["bg_ratio"], "5/8")            # background: genes present in the library

    def test_custom_enrichment_with_networks(self):
        tmp = Path(tempfile.mkdtemp())
        genes = [f"G{i}" for i in range(1, 41)]
        common.write_tsv(tmp / "sets.tsv", ("Cluster", "SYMBOL"), [("all", g) for g in genes[:12]])
        common.write_tsv(tmp / "universe.tsv", ("SYMBOL",), [(f"G{i}",) for i in range(1, 2001)])
        gmt = tmp / "sig.gmt"
        gmt.write_text("".join(f"{name}\tna\t" + "\t".join(members) + "\n" for name, members in (
            ("S1", genes[:8]), ("S2", genes[2:10]), ("S3", genes[4:12]), ("S4", genes[30:40]))))
        res = engine.run_enrichment(str(tmp / "sets.tsv"), tmp / "out", "t", methods=["custom"],
                                    universe=str(tmp / "universe.tsv"), custom_gmt=str(gmt))
        json.dumps(res)
        self.assertEqual(res["n_significant"], [{"method": "custom", "group": "all", "N": 3}])
        self.assertEqual(res["top_terms"][0]["count"], 8)
        self.assertIn("map_custom", res["figures"])
        self.assertIn("cnet_custom", res["figures"])
        nodes = common.read_tsv(res["tables"]["network_nodes_custom"])
        self.assertEqual({n["term"] for n in nodes}, {"S1", "S2", "S3"})
        self.assertTrue(all(n["x"] and n["y"] for n in nodes))
        self.assertTrue(res["hub_genes"]["custom"])
        with self.assertRaisesRegex(RuntimeError, "All pathway methods failed"):
            engine.run_enrichment(str(tmp / "sets.tsv"), tmp / "out", "t", methods=["custom"])

    def test_modules_group_overlapping_terms(self):
        rows = [{"method": "m", "group": "g", "term": t, "padj": p, "count": len(g.split(";")), "genes": g}
                for t, p, g in (("A1", 1e-9, "a;b;c;d"), ("A2", 1e-8, "a;b;c;e"), ("A3", 1e-7, "a;b;d;e"),
                                ("B1", 1e-6, "x;y;z;w"), ("B2", 1e-5, "x;y;z;v"))]
        em = networks.build_enrichment_map(rows)
        module = {n["term"]: n["module"] for n in em["nodes"]}
        self.assertEqual(module["A1"], module["A2"])
        self.assertEqual(module["B1"], module["B2"])
        self.assertNotEqual(module["A1"], module["B1"])
        self.assertEqual(module["A1"], 1)                            # modules are ordered by their best term
        self.assertEqual({n["module_label"] for n in em["nodes"]}, {"A1", "B1"})
        self.assertIsNone(networks.build_enrichment_map(rows[:2]))

    def test_method_names(self):
        self.assertEqual(enrich.canonical_method("msigdb_H"), "hallmark")
        self.assertEqual(enrich.canonical_method("msigdb_C7:IMMUNESIGDB"), "immunesigdb")
        self.assertEqual(enrich.canonical_method("GO_BP"), "go_bp")
        with self.assertRaisesRegex(ValueError, "Unknown enrichment method"):
            enrich.library("nonsense", "mm10")


class QcAndMotifTest(unittest.TestCase):
    def test_qc_cleans_contigs_duplicates_and_blacklist(self):
        tmp = Path(tempfile.mkdtemp())
        (tmp / "p.bed").write_text("chr1\t100\t400\ta\nchr1\t100\t400\tb\nchrUn_JH584304\t5\t50\tc\n"
                                   "chr2\t1000\t1300\td\nchr3\t500\t900\te\n")
        (tmp / "bl.bed").write_text("chr2\t1200\t5000\n")
        res = engine.run_peak_qc(read_peak_file(tmp / "p.bed"), tmp / "qc", "t", blacklist=str(tmp / "bl.bed"))
        self.assertEqual((res["n_peaks"], res["n_duplicated"], res["n_nonstandard_chrom"], res["n_blacklisted"],
                          res["n_clean_peaks"]), (5, 1, 1, 1, 2))
        self.assertTrue(Path(res["figures"]["width"]).is_file())
        self.assertTrue(any("Only 5 peaks" in w for w in res["warnings"]))

    def test_blacklist_overlap_is_half_open(self):
        from peakannoagent.engine.qc import blacklist_overlaps
        tmp = Path(tempfile.mkdtemp())
        (tmp / "bl.bed").write_text("chr1\t1\t4\n")
        peaks = Peaks(["left", "hit", "hit2", "right"], np.array(["chr1"] * 4, dtype=object),
                      np.array([0, 1, 2, 4]), np.array([1, 2, 4, 5]))
        self.assertEqual(list(blacklist_overlaps(peaks, str(tmp / "bl.bed"))), [False, True, True, False])

    def test_homer_results_are_standardised(self):
        from peakannoagent.engine.motifs import motif_tf_name, read_homer_known_results
        self.assertEqual(motif_tf_name("CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer"), "CTCF")
        self.assertEqual(motif_tf_name("Nr4a1"), "Nr4a1")
        path = Path(tempfile.mkdtemp()) / "knownResults.txt"
        path.write_text("\t".join(["Motif Name", "Consensus", "P-value", "Log P-value", "q-value (Benjamini)",
                                   "# of Target Sequences with Motif(of 300)", "% of Target Sequences with Motif",
                                   "# of Background Sequences with Motif(of 40000)",
                                   "% of Background Sequences with Motif"]) + "\n" +
                        "\t".join(["BATF(bZIP)/Th17-BATF-ChIP-Seq/Homer", "DATGASTCAT", "1e-20", "-46.05", "0.0000",
                                   "90", "30.00%", "4000", "10.00%"]) + "\n")
        row = read_homer_known_results(path)[0]
        self.assertEqual((row["tf"], row["family"], row["pct_target"]), ("BATF", "bZIP", 30.0))
        self.assertGreater(row["log2_enrichment"], 1.5)


if __name__ == "__main__":
    unittest.main()
