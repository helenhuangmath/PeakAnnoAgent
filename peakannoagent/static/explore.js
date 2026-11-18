/* PeakAnnoAgent pathway explorer: interactive top-term chart, pathway network and gene-pathway map.
   No dependencies. Used by the web run page (data from /api/runs/<id>/explore/<E>) and by the
   self-contained HTML report (data inlined as JSON).

   PeakExplore.mount(container, data, opts)
     data  = {enrichment: [{method, group, term, padj, count, genes}], networks: {key: {nodes, edges}}}
     opts  = {onGene(symbol) -> void (optional), geneInfo(symbol) -> Promise<object> (optional), title}
*/
(function () {
  const NS = "http://www.w3.org/2000/svg";
  const PALETTE = ["#18bc9c", "#3498db", "#e74c3c", "#f39c12", "#9b59b6", "#2c3e50", "#1abc9c", "#e67e22", "#34495e", "#16a085", "#c0392b", "#8e44ad"];
  const esc = (s) => String(s == null ? "" : s).replace(/[&<>"']/g, (c) => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"}[c]));
  const el = (tag, attrs, ...kids) => {
    const e = tag.startsWith("svg:") ? document.createElementNS(NS, tag.slice(4)) : document.createElement(tag);
    for (const [k, v] of Object.entries(attrs || {})) {
      if (k === "text") e.textContent = v;
      else if (k === "html") e.innerHTML = v;
      else if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
      else if (v !== null && v !== undefined) e.setAttribute(k, v);
    }
    kids.forEach((k) => k && e.appendChild(typeof k === "string" ? document.createTextNode(k) : k));
    return e;
  };
  const genesOf = (s) => String(s || "").split(/[;,\s]+/).map((g) => g.trim()).filter(Boolean);
  const nlog = (p) => { const v = parseFloat(p); return isNaN(v) || v <= 0 ? 0 : -Math.log10(v); };
  const fmtP = (p) => { const v = parseFloat(p); return isNaN(v) ? "" : (v < 1e-3 ? v.toExponential(2) : v.toFixed(3)); };
  const download = (name, text, mime) => {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([text], {type: mime || "text/plain"}));
    a.download = name; document.body.appendChild(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 100);
  };
  function svgText(svg) {
    const c = svg.cloneNode(true);
    c.setAttribute("xmlns", NS);
    c.querySelectorAll("*").forEach((n) => { n.removeAttribute("class"); });
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(c);
  }
  function svgToPng(svg, name) {
    const w = +svg.getAttribute("width"), h = +svg.getAttribute("height");
    const img = new Image();
    img.onload = () => {
      const cv = document.createElement("canvas"); cv.width = w * 2; cv.height = h * 2;
      const ctx = cv.getContext("2d"); ctx.scale(2, 2); ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, w, h); ctx.drawImage(img, 0, 0);
      cv.toBlob((b) => { const a = document.createElement("a"); a.href = URL.createObjectURL(b); a.download = name; a.click(); }, "image/png");
    };
    img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(svgText(svg));
  }
  // Inline the few styles the SVG needs so downloads look like the page.
  const SVG_STYLE = "text{font-family:Lato,Helvetica,Arial,sans-serif;fill:#2c3e50}.axis line,.axis path{stroke:#dce4ec}.bar:hover{opacity:.8}.node{cursor:pointer}.edge{stroke:#95a5a6;stroke-opacity:.5}";

  // ------------------------------------------------------------------ force layout (Fruchterman-Reingold, deterministic)
  function layout(nodes, edges, w, h, iters, o) {
    o = Object.assign({k: 0.45, gravity: 0.1, pull: 0.4}, o || {});
    const n = nodes.length; if (!n) return nodes;
    const idx = new Map(nodes.map((d, i) => [d.id, i]));
    let seed = 7; const rnd = () => (seed = (seed * 9301 + 49297) % 233280) / 233280;
    nodes.forEach((d, i) => { const a = 2 * Math.PI * i / n; const rr = 0.25 + 0.6 * rnd(); d.x = w / 2 + (w / 2.6) * rr * Math.cos(a); d.y = h / 2 + (h / 2.6) * rr * Math.sin(a); });
    const links = edges.map((e) => [idx.get(e.source), idx.get(e.target), e.weight || 1]).filter((l) => l[0] !== undefined && l[1] !== undefined && l[0] !== l[1]);
    const k = Math.sqrt((w * h) / n) * o.k;
    const total = iters || 300;
    const dx = new Float64Array(n), dy = new Float64Array(n);
    for (let t = 0; t < total; t++) {
      const temp = Math.max(1, (w / 8) * (1 - t / total));
      dx.fill(0); dy.fill(0);
      for (let i = 0; i < n; i++) {
        const a = nodes[i];
        for (let j = i + 1; j < n; j++) {
          const b = nodes[j];
          let ddx = a.x - b.x, ddy = a.y - b.y; let d = Math.sqrt(ddx * ddx + ddy * ddy) || 0.01;
          const min = a.r + b.r + 6;                       // never overlap
          const f = (k * k) / d * (d < min ? 3 : 1);
          ddx /= d; ddy /= d;
          dx[i] += ddx * f; dy[i] += ddy * f; dx[j] -= ddx * f; dy[j] -= ddy * f;
        }
      }
      for (const [i, j, wgt] of links) {
        const a = nodes[i], b = nodes[j];
        let ddx = b.x - a.x, ddy = b.y - a.y; const d = Math.sqrt(ddx * ddx + ddy * ddy) || 0.01;
        const f = (d * d) / k * (o.pull + Math.min(1, wgt));
        ddx /= d; ddy /= d;
        dx[i] += ddx * f; dy[i] += ddy * f; dx[j] -= ddx * f; dy[j] -= ddy * f;
      }
      for (let i = 0; i < n; i++) {
        const a = nodes[i]; if (a.fixed) continue;
        // gravity keeps disconnected pieces on the canvas
        const gx = w / 2 - a.x, gy = h / 2 - a.y; const gd = Math.sqrt(gx * gx + gy * gy) || 0.01;
        dx[i] += gx * o.gravity * (a.kind === "gene" ? 0.5 : 1); dy[i] += gy * o.gravity * (a.kind === "gene" ? 0.5 : 1);
        const len = Math.sqrt(dx[i] * dx[i] + dy[i] * dy[i]) || 0.01;
        const step = Math.min(len, temp);
        a.x = Math.max(a.r + 4, Math.min(w - a.r - 4, a.x + dx[i] / len * step));
        a.y = Math.max(a.r + 14, Math.min(h - a.r - 16, a.y + dy[i] / len * step));
      }
    }
    return nodes;
  }

  // ------------------------------------------------------------------ the explorer
  function mount(root, data, opts) {
    opts = opts || {};
    root.classList.add("explore");
    const enrich = (data.enrichment || []).map((r) => ({...r, genes_list: genesOf(r.genes), score: nlog(r.padj)}));
    const networks = data.networks || {};
    if (!enrich.length && !Object.keys(networks).length) { root.innerHTML = '<p class="muted">No pathway data to explore.</p>'; return; }
    const methods = [...new Set(enrich.map((r) => r.method))];
    const groups = [...new Set(enrich.map((r) => r.group))];
    const state = {view: opts.view || "terms", method: methods[0] || Object.keys(networks)[0], group: groups.includes("all") ? "all" : groups[0], top: 20, selected: null, gene: null};

    // toolbar
    const tabs = el("div", {class: "ex-tabs"});
    const views = [["terms", "Top terms"], ["network", "Pathway network"], ["genes", "Gene–pathway map"]];
    const sel = (name, values, labelFn) => el("select", {onchange: (e) => { state[name] = e.target.value; render(); }},
      ...values.map((v) => el("option", {value: v, text: labelFn ? labelFn(v) : v, selected: String(state[name]) === String(v) ? "selected" : null})));
    const controls = el("div", {class: "ex-controls"});
    const stage = el("div", {class: "ex-stage"});
    const panel = el("aside", {class: "ex-panel"});
    const dl = el("div", {class: "ex-dl"});
    root.append(el("div", {class: "ex-head"}, tabs, controls, dl), el("div", {class: "ex-body"}, stage, panel));

    function toolbar() {
      tabs.innerHTML = "";
      views.forEach(([k, label]) => {
        if (k === "network" && !Object.keys(networks).length) return;
        tabs.appendChild(el("button", {type: "button", class: "ex-tab" + (state.view === k ? " on" : ""), text: label, onclick: () => { state.view = k; render(); }}));
      });
      controls.innerHTML = "";
      const netKeys = Object.keys(networks);
      if (state.view === "network") {
        if (!netKeys.includes(state.method)) state.method = netKeys[0];
        controls.append(el("label", {}, "Collection ", sel("method", netKeys)));
      } else {
        if (!methods.includes(state.method)) state.method = methods[0];
        controls.append(el("label", {}, "Collection ", sel("method", methods)));
        if (groups.length > 1) controls.append(el("label", {}, "Peaks ", sel("group", groups)));
        controls.append(el("label", {}, "Show ", sel("top", [10, 20, 30, 50], (v) => `top ${v}`)));
      }
      dl.innerHTML = "";
      const svg = () => stage.querySelector("svg");
      const name = `${(opts.title || "pathways").replace(/\W+/g, "_")}_${state.view}_${String(state.method).replace(/\W+/g, "_")}`;
      dl.append(el("span", {class: "muted small", text: "Download"}),
        el("button", {type: "button", class: "ex-link", text: "SVG", onclick: () => svg() && download(name + ".svg", svgText(svg()), "image/svg+xml")}),
        el("button", {type: "button", class: "ex-link", text: "PNG", onclick: () => svg() && svgToPng(svg(), name + ".png")}),
        el("button", {type: "button", class: "ex-link", text: "TSV", onclick: () => download(name + ".tsv", currentTsv(), "text/tab-separated-values")}));
    }

    function rows() {
      return enrich.filter((r) => r.method === state.method && (groups.length < 2 || r.group === state.group))
        .sort((a, b) => b.score - a.score).slice(0, +state.top);
    }
    function currentTsv() {
      if (state.view === "network") {
        const net = networks[state.method] || {nodes: []};
        return ["term\tmodule\tmodule_label\tpadj\tcount\tgenes", ...net.nodes.map((n) => [n.term, n.module, n.module_label, n.padj, n.count, n.genes].join("\t"))].join("\n");
      }
      return ["method\tgroup\tterm\tpadj\tcount\tgenes", ...rows().map((r) => [r.method, r.group, r.term, r.padj, r.count, r.genes_list.join(";")].join("\t"))].join("\n");
    }

    // --- gene panel
    function showTerm(t) {
      state.selected = t; state.gene = null;
      panel.innerHTML = "";
      if (!t) { panel.append(el("p", {class: "muted", text: "Click a bar, a node or a term to see its genes."})); return; }
      const genes = t.genes_list || genesOf(t.genes);
      panel.append(el("h4", {text: t.term}),
        el("p", {class: "muted small", html: `${esc(t.method || state.method)}${t.group ? " · " + esc(t.group) + " peaks" : ""}${t.module_label ? " · module: " + esc(t.module_label) : ""}<br>padj ${esc(fmtP(t.padj))} · ${genes.length} gene${genes.length === 1 ? "" : "s"}${t.gene_ratio ? " · ratio " + esc(t.gene_ratio) : ""}`}));
      const chips = el("div", {class: "ex-genes"});
      genes.forEach((g) => chips.append(el("button", {type: "button", class: "ex-gene", text: g, title: "gene details", onclick: () => showGene(g, t)})));
      panel.append(chips);
      panel.append(el("div", {class: "ex-actions"},
        el("button", {type: "button", class: "ex-link", text: "copy genes", onclick: () => navigator.clipboard && navigator.clipboard.writeText(genes.join("\n"))}),
        el("button", {type: "button", class: "ex-link", text: "genes.tsv", onclick: () => download(`${t.term.replace(/\W+/g, "_")}_genes.tsv`, "gene\n" + genes.join("\n"))})));
      // which other shown terms share these genes
      const others = rows().filter((r) => r !== t && r.term !== t.term).map((r) => [r, r.genes_list.filter((g) => genes.includes(g)).length]).filter((x) => x[1] > 0).sort((a, b) => b[1] - a[1]).slice(0, 6);
      if (others.length) panel.append(el("h5", {text: "Shares genes with"}), el("ul", {class: "ex-list"}, ...others.map(([r, n]) => el("li", {}, el("a", {href: "#", text: r.term, onclick: (e) => { e.preventDefault(); showTerm(r); highlight(r); }}), ` (${n})`))));
    }
    async function showGene(g, fromTerm) {
      state.gene = g;
      const box = el("div", {class: "ex-geneinfo"}, el("b", {text: g}), " ");
      const inTerms = enrich.filter((r) => r.method === state.method && r.genes_list.includes(g)).map((r) => r.term);
      const uniq = [...new Set(inTerms)];
      box.append(el("span", {class: "muted small", text: `in ${uniq.length} ${state.method} term${uniq.length === 1 ? "" : "s"} shown in this run`}));
      if (opts.onGene) box.append(" ", el("a", {href: "#", class: "ex-link", text: "search this run →", onclick: (e) => { e.preventDefault(); opts.onGene(g); }}));
      const old = panel.querySelector(".ex-geneinfo"); if (old) old.remove();
      panel.insertBefore(box, panel.querySelector(".ex-genes").nextSibling);
      if (opts.geneInfo) {
        try {
          const d = await opts.geneInfo(g);
          if (d && d.found !== false) {
            const bits = [];
            if (d.peaks && d.peaks.length) bits.push(`${d.peaks.length} peak(s): ` + d.peaks.slice(0, 4).map((p) => `${p.PeakID || p.peak || ""} (${p.region_class || ""}${p.distance_to_tss !== undefined ? ", " + p.distance_to_tss + " bp" : ""})`).join("; "));
            if (d.rnaseq) bits.push(`RNA-seq ${d.rnaseq.direction || ""} log2FC ${d.rnaseq.log2FoldChange ?? d.rnaseq.log2fc ?? ""}`);
            if (d.control) bits.push("expected gene (positive control)");
            if (d.pathways && d.pathways.length) bits.push(`${d.pathways.length} pathway(s) in all collections`);
            if (bits.length) box.append(el("ul", {class: "ex-list small"}, ...bits.map((b) => el("li", {text: b}))));
          }
        } catch (e) { /* offline or no index */ }
      }
      uniq.filter((t) => !fromTerm || t !== fromTerm.term).slice(0, 8).forEach((t) => box.append(el("a", {href: "#", class: "ex-term", text: t, onclick: (e) => { e.preventDefault(); const r = enrich.find((x) => x.term === t && x.method === state.method); if (r) { showTerm(r); highlight(r); } }})));
    }
    function highlight(t) {
      stage.querySelectorAll(".sel").forEach((n) => n.classList.remove("sel"));
      if (!t) return;
      stage.querySelectorAll(`[data-term]`).forEach((n) => { if (n.dataset.term === t.term) n.classList.add("sel"); });
    }

    // --- views
    function drawTerms() {
      const rs = rows();
      stage.innerHTML = "";
      if (!rs.length) { stage.append(el("p", {class: "muted", text: "No significant terms for this selection."})); return; }
      const W = Math.max(560, Math.min(stage.clientWidth || 760, 900)), left = Math.min(360, W * 0.45), rowH = 22, H = rs.length * rowH + 50;
      const max = Math.max(...rs.map((r) => r.score), 1);
      const x = (v) => left + (v / max) * (W - left - 60);
      const svg = el("svg:svg", {width: W, height: H, viewBox: `0 0 ${W} ${H}`, class: "ex-svg"});
      svg.append(el("svg:style", {text: SVG_STYLE}));
      svg.append(el("svg:text", {x: left, y: 14, "font-size": 11, fill: "#7b8a8b", text: `−log10 adjusted p · ${state.method}${groups.length > 1 ? " · " + state.group + " peaks" : ""} · bar label = genes`}));
      const ticks = 4;
      for (let i = 0; i <= ticks; i++) {
        const v = max * i / ticks, px = x(v);
        svg.append(el("svg:line", {x1: px, x2: px, y1: 26, y2: H - 20, stroke: "#ecf0f1"}));
        svg.append(el("svg:text", {x: px, y: H - 6, "font-size": 10, "text-anchor": "middle", fill: "#7b8a8b", text: v.toFixed(0)}));
      }
      rs.forEach((r, i) => {
        const y = 30 + i * rowH;
        const g = el("svg:g", {class: "bar node", "data-term": r.term, onclick: () => { showTerm(r); highlight(r); }});
        g.append(el("svg:title", {text: `${r.term}\npadj ${fmtP(r.padj)} · ${r.count || r.genes_list.length} genes\n${r.genes_list.slice(0, 12).join(", ")}${r.genes_list.length > 12 ? " …" : ""}`}));
        const label = r.term.length > left / 6.2 ? r.term.slice(0, left / 6.2 - 1) + "…" : r.term;
        g.append(el("svg:text", {x: left - 8, y: y + rowH / 2 + 4, "font-size": 12, "text-anchor": "end", text: label}));
        g.append(el("svg:rect", {x: left, y: y + 3, width: Math.max(2, x(r.score) - left), height: rowH - 6, rx: 2, fill: r.group === "distal" ? "#3498db" : r.group === "promoter" ? "#e74c3c" : "#18bc9c"}));
        g.append(el("svg:text", {x: x(r.score) + 5, y: y + rowH / 2 + 4, "font-size": 10, fill: "#7b8a8b", text: String(r.count || r.genes_list.length)}));
        svg.append(g);
      });
      stage.append(svg);
      stage.append(el("p", {class: "muted small", text: "Click a bar to list its genes. Hover for the first genes."}));
    }

    // Module colours and layout match the R figure (SimpNetworkViz style): when the node table carries
    // x/y (inches, y up) the modules keep their packed positions, each in a tinted outline.
    const NET_PALETTE = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
      "#7E6148", "#B09C85", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#EE4C97", "#0072B5", "#BC3C29", "#925E9F", "#FDAF91", "#AD002A"];
    const mixHex = (hex, to, a) => "#" + [1, 3, 5].map((i) => { const c = parseInt(hex.slice(i, i + 2), 16); return Math.round(c + (to - c) * a).toString(16).padStart(2, "0"); }).join("");
    function convexHull(pts) {
      const p = pts.slice().sort((a, b) => a[0] - b[0] || a[1] - b[1]);
      if (p.length < 3) return p;
      const cross = (o, a, b) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
      const lo = [], up = [];
      for (const q of p) { while (lo.length >= 2 && cross(lo[lo.length - 2], lo[lo.length - 1], q) <= 0) lo.pop(); lo.push(q); }
      for (const q of p.slice().reverse()) { while (up.length >= 2 && cross(up[up.length - 2], up[up.length - 1], q) <= 0) up.pop(); up.push(q); }
      return lo.slice(0, -1).concat(up.slice(0, -1));
    }

    function drawNetwork() {
      const net = networks[state.method];
      stage.innerHTML = "";
      if (!net || !net.nodes || !net.nodes.length) { stage.append(el("p", {class: "muted", text: "No network for this collection."})); return; }
      const W = Math.max(560, Math.min(stage.clientWidth || 760, 900));
      const maxC = Math.max(...net.nodes.map((n) => +n.count || 1));
      const nodes = net.nodes.map((n) => ({...n, id: n.term, genes_list: genesOf(n.genes), method: state.method, r: 5 + 11 * Math.sqrt((+n.count || 1) / maxC)}));
      const edges = (net.edges || []).map((e) => ({source: e.from, target: e.to, weight: +e.jaccard || 0.2}));
      const sizes = {};
      nodes.forEach((n) => { sizes[n.module] = (sizes[n.module] || 0) + 1; });
      const clustered = (n) => sizes[n.module] >= 2;
      const color = (n) => clustered(n) ? NET_PALETTE[(+n.module - 1) % NET_PALETTE.length] : "#BDBDBD";
      const fixed = nodes.every((n) => n.x !== undefined && n.x !== null && n.x !== "" && isFinite(+n.x) && isFinite(+n.y));
      let H;
      if (fixed) {
        const pad = 70, xs = nodes.map((n) => +n.x), ys = nodes.map((n) => +n.y);
        const x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...ys), y1 = Math.max(...ys);
        const sc = Math.min((W - 2 * pad) / Math.max(x1 - x0, 1e-6), 95);
        H = Math.round((y1 - y0) * sc + 2 * pad + 20);
        nodes.forEach((n) => { n.x = pad + (+n.x - x0) * sc + (W - 2 * pad - (x1 - x0) * sc) / 2; n.y = pad + 20 + (y1 - +n.y) * sc; });
      } else {
        H = Math.max(420, Math.min(640, 160 + nodes.length * 12));
        layout(nodes, edges, W, H, 300, {k: 0.45, gravity: 0.1, pull: 0.4});
      }
      const svg = el("svg:svg", {width: W, height: H, viewBox: `0 0 ${W} ${H}`, class: "ex-svg"});
      svg.append(el("svg:style", {text: SVG_STYLE}));
      const byId = new Map(nodes.map((n) => [n.id, n]));
      const label = (t) => t.length > 30 ? t.slice(0, 29) + "…" : t;
      // Same wrapping as the R figure (18 characters, at most 4 lines), which the layout made room for.
      const wrap = (t) => {
        const lines = [];
        t.split(/\s+/).forEach((w) => { const l = lines.length ? lines[lines.length - 1] : null; if (l !== null && (l + " " + w).length <= 18) lines[lines.length - 1] = l + " " + w; else lines.push(w); });
        if (lines.length > 4) { lines.length = 4; lines[3] = lines[3].replace(/[ ,;:-]*$/, "") + "…"; }
        return lines;
      };
      const FS = 10, LH = 11.5;
      nodes.forEach((n) => { n.lines = fixed ? wrap(n.term) : [label(n.term)]; n.lw = Math.max(...n.lines.map((l) => l.length)) * 0.58 * FS / 2; });
      // Outlines and titles per module (only for the precomputed layout, where modules are apart).
      const mods = [...new Set(nodes.map((n) => String(n.module)))];
      const hg = el("svg:g");
      if (fixed) {
        mods.forEach((m) => {
          const ms = nodes.filter((n) => String(n.module) === m);
          if (ms.length < 2) return;
          const pts = [];
          ms.forEach((n) => { const w = Math.max(n.r, n.lw) + 8, bot = n.y + n.r + 6 + n.lines.length * LH; for (let a = 0; a < 16; a++) pts.push([n.x + (n.r + 10) * Math.cos(a * Math.PI / 8), n.y + (n.r + 10) * Math.sin(a * Math.PI / 8)]); pts.push([n.x - w, bot], [n.x + w, bot]); });
          const hull = convexHull(pts), base = color(ms[0]);
          hg.append(el("svg:polygon", {points: hull.map((q) => q.join(",")).join(" "), fill: mixHex(base, 255, 0.86), stroke: mixHex(base, 255, 0.25), "stroke-width": 1.4, "stroke-linejoin": "round"}));
          const top = Math.min(...hull.map((q) => q[1])), xs = hull.map((q) => q[0]);
          hg.append(el("svg:text", {x: (Math.min(...xs) + Math.max(...xs)) / 2, y: top - 7, "text-anchor": "middle", "font-size": 13, style: `fill:${mixHex(base, 0, 0.45)}`, text: `M${m}: ${label(ms[0].module_label || "")}`}));
        });
      }
      svg.append(hg);
      const eg = el("svg:g");
      edges.forEach((e) => {
        const a = byId.get(e.source), b = byId.get(e.target);
        if (!a || !b) return;
        const within = String(a.module) === String(b.module);
        eg.append(el("svg:line", {class: "edge", x1: a.x, y1: a.y, x2: b.x, y2: b.y, "stroke-width": within ? 0.5 + 3 * e.weight : 0.6, style: within ? "" : "stroke-opacity:.18", "data-a": a.id, "data-b": b.id}));
      });
      svg.append(eg);
      const ng = el("svg:g");
      nodes.forEach((n) => {
        const g = el("svg:g", {class: "node", "data-term": n.term, transform: `translate(${n.x},${n.y})`, onclick: () => { showTerm({...n, module_label: n.module_label}); highlight(n); }});
        g.append(el("svg:title", {text: `${n.term}\nmodule ${n.module}: ${n.module_label || ""}\npadj ${fmtP(n.padj)} · ${n.count} genes`}));
        g.append(el("svg:circle", {r: n.r, fill: color(n), stroke: mixHex(color(n), 0, 0.3), "stroke-width": 0.8}));
        if (fixed || n.r > 11 || nodes.length <= 40) {
          const t = el("svg:text", {y: n.r + 11, "font-size": FS, "text-anchor": "middle"});
          n.lines.forEach((l, i) => t.append(el("svg:tspan", {x: 0, dy: i ? LH : 0, text: l})));
          g.append(t);
        }
        drag(g, n, () => { g.setAttribute("transform", `translate(${n.x},${n.y})`); eg.querySelectorAll(`line[data-a="${CSS.escape(n.id)}"]`).forEach((l) => { l.setAttribute("x1", n.x); l.setAttribute("y1", n.y); }); eg.querySelectorAll(`line[data-b="${CSS.escape(n.id)}"]`).forEach((l) => { l.setAttribute("x2", n.x); l.setAttribute("y2", n.y); }); });
        ng.append(g);
      });
      svg.append(ng);
      stage.append(svg);
      const legend = el("div", {class: "ex-legend"});
      mods.forEach((m) => { const n0 = nodes.find((n) => String(n.module) === m) || {}; if (fixed && !clustered(n0)) return; legend.append(el("span", {}, el("i", {style: `background:${color(n0)}`}), ` M${m}: ${n0.module_label || ""}`)); });
      stage.append(legend, el("p", {class: "muted small", text: "Nodes are pathways (size = genes, colour = module); edges link pathways that share genes (width = Jaccard; faint lines join different modules). Grey nodes share too few genes to join a module. Drag to rearrange, click for genes."}));
    }

    function drawGeneMap() {
      const rs = rows();
      stage.innerHTML = "";
      if (!rs.length) { stage.append(el("p", {class: "muted", text: "No terms for this selection."})); return; }
      const terms = rs.slice(0, Math.min(rs.length, 15));
      const geneCount = new Map();
      terms.forEach((t) => t.genes_list.forEach((g) => geneCount.set(g, (geneCount.get(g) || 0) + 1)));
      const genes = [...geneCount.keys()];
      const W = Math.max(560, Math.min(stage.clientWidth || 760, 900)), H = Math.max(520, Math.min(760, 320 + genes.length * 3));
      const nodes = [
        ...terms.map((t) => ({id: "T:" + t.term, kind: "term", data: t, r: 9 + 10 * Math.sqrt(t.genes_list.length / Math.max(1, Math.max(...terms.map((x) => x.genes_list.length)))), label: t.term})),
        ...genes.map((g) => ({id: "G:" + g, kind: "gene", data: g, r: Math.min(9, 3 + 1.1 * geneCount.get(g)), label: g})),
      ];
      const edges = []; terms.forEach((t) => t.genes_list.forEach((g) => edges.push({source: "T:" + t.term, target: "G:" + g, weight: 0.05})));
      layout(nodes, edges, W, H, 300, {k: 0.7, gravity: 0.06, pull: 0.08});
      const svg = el("svg:svg", {width: W, height: H, viewBox: `0 0 ${W} ${H}`, class: "ex-svg"});
      svg.append(el("svg:style", {text: SVG_STYLE}));
      const byId = new Map(nodes.map((n) => [n.id, n]));
      const eg = el("svg:g");
      edges.forEach((e) => { const a = byId.get(e.source), b = byId.get(e.target); eg.append(el("svg:line", {class: "edge", x1: a.x, y1: a.y, x2: b.x, y2: b.y, "stroke-width": 0.7, "data-a": a.id, "data-b": b.id})); });
      svg.append(eg);
      const ng = el("svg:g");
      nodes.forEach((n, i) => {
        const isT = n.kind === "term";
        const g = el("svg:g", {class: "node" + (isT ? " term" : " gene"), "data-term": isT ? n.data.term : null, transform: `translate(${n.x},${n.y})`,
          onclick: () => { if (isT) { showTerm(n.data); highlight(n.data); } else { showTerm(terms.find((t) => t.genes_list.includes(n.data))); showGene(n.data); } }});
        g.append(el("svg:title", {text: isT ? `${n.data.term}\npadj ${fmtP(n.data.padj)} · ${n.data.genes_list.length} genes` : `${n.data}\nin ${geneCount.get(n.data)} of the ${terms.length} terms shown`}));
        g.append(el("svg:circle", {r: n.r, fill: isT ? PALETTE[terms.indexOf(n.data) % PALETTE.length] : (geneCount.get(n.data) > 1 ? "#2c3e50" : "#95a5a6"), stroke: "#fff", "stroke-width": 1}));
        if (isT || geneCount.get(n.data) > 1 || genes.length <= 60) g.append(el("svg:text", {y: n.r + 10, "font-size": isT ? 11 : 9, "font-weight": isT ? 700 : 400, "text-anchor": "middle", text: n.label.length > 30 ? n.label.slice(0, 29) + "…" : n.label}));
        drag(g, n, () => { g.setAttribute("transform", `translate(${n.x},${n.y})`); eg.querySelectorAll(`line[data-a="${CSS.escape(n.id)}"]`).forEach((l) => { l.setAttribute("x1", n.x); l.setAttribute("y1", n.y); }); eg.querySelectorAll(`line[data-b="${CSS.escape(n.id)}"]`).forEach((l) => { l.setAttribute("x2", n.x); l.setAttribute("y2", n.y); }); });
        ng.append(g);
      });
      svg.append(ng);
      stage.append(svg, el("p", {class: "muted small", text: `${terms.length} top terms (coloured) and their ${genes.length} genes. Dark genes sit in several terms. Click a gene for its terms, a term for its genes; drag to rearrange.`}));
    }

    function drag(g, n, redraw) {
      let start = null;
      g.addEventListener("pointerdown", (e) => { start = {x: e.clientX, y: e.clientY, nx: n.x, ny: n.y}; g.setPointerCapture(e.pointerId); e.preventDefault(); });
      g.addEventListener("pointermove", (e) => { if (!start) return; const svg = g.ownerSVGElement; const sc = svg.getBoundingClientRect().width / +svg.getAttribute("width"); n.x = start.nx + (e.clientX - start.x) / sc; n.y = start.ny + (e.clientY - start.y) / sc; n.fixed = true; redraw(); });
      g.addEventListener("pointerup", () => { start = null; });
    }

    function render() {
      toolbar();
      if (state.view === "network") drawNetwork(); else if (state.view === "genes") drawGeneMap(); else drawTerms();
      if (state.selected) highlight(state.selected); else showTerm(null);
    }
    render();
    return {render, state};
  }

  window.PeakExplore = {mount};
})();
