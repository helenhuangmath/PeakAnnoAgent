(function () {
  const body = document.body;
  if (!body.classList.contains("runpage")) return;
  const jobId = body.dataset.job;
  const api = (p) => `/api/runs/${jobId}${p}`;
  const $ = (s) => document.querySelector(s);
  const esc = (s) => String(s == null ? "" : s).replace(/[&<>"']/g, (c) => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"}[c]));

  // ---------------------------------------------------------------- parts (sidebar navigation)
  const nav = document.querySelectorAll(".parts a[data-part]");
  const parts = document.querySelectorAll(".part");
  function showPart(name, scrollTo) {
    if (!document.getElementById(`part-${name}`)) return false;
    parts.forEach((p) => (p.hidden = p.id !== `part-${name}`));
    nav.forEach((a) => a.classList.toggle("active", a.dataset.part === name));
    if (history.replaceState) history.replaceState(null, "", `#${name}`);
    if (name === "chat") loadChat();
    if (scrollTo) {
      const el = document.getElementById(scrollTo);
      if (el) { el.scrollIntoView({behavior: "smooth", block: "start"}); el.classList.add("flash"); setTimeout(() => el.classList.remove("flash"), 1500); }
    } else {
      window.scrollTo({top: 0});
    }
    return true;
  }
  nav.forEach((a) => a.addEventListener("click", (e) => { e.preventDefault(); showPart(a.dataset.part); }));

  // Evidence ids (E4) and section anchors anywhere on the page open the part that holds them.
  function partOf(id) {
    const el = document.getElementById(id);
    const sec = el && el.closest(".part");
    return sec ? sec.id.replace(/^part-/, "") : null;
  }
  document.addEventListener("click", (e) => {
    const a = e.target.closest("a[href^='#']");
    if (!a || a.dataset.part !== undefined && a.closest(".parts")) return;
    const id = a.getAttribute("href").slice(1);
    if (!id) return;
    if (a.dataset.part) { e.preventDefault(); showPart(a.dataset.part); return; }
    const part = partOf(id);
    if (part) { e.preventDefault(); showPart(part, id); }
  });
  const params = new URLSearchParams(location.search);
  const initial = params.get("tab") || location.hash.replace("#", "");
  if (initial) {
    const part = document.getElementById(`part-${initial}`) ? initial : partOf(initial);
    if (part) showPart(part, part === initial ? null : initial);
  }

  // ---------------------------------------------------------------- interactive tables (sort, filter, show all)
  const num = (t) => { const v = parseFloat(String(t).replace(/,/g, "")); return isNaN(v) ? null : v; };
  document.addEventListener("click", (e) => {
    const th = e.target.closest(".itable th[data-col]");
    if (th) {
      const table = th.closest("table"), tbody = table.tBodies[0], col = +th.dataset.col;
      const asc = th.dataset.dir !== "asc";
      table.querySelectorAll("th").forEach((h) => { delete h.dataset.dir; h.classList.remove("asc", "desc"); });
      th.dataset.dir = asc ? "asc" : "desc"; th.classList.add(asc ? "asc" : "desc");
      const rows = [...tbody.rows];
      const allNum = rows.every((r) => { const t = r.cells[col].textContent.trim(); return t === "" || num(t) !== null; });
      rows.sort((a, b) => {
        const x = a.cells[col].textContent.trim(), y = b.cells[col].textContent.trim();
        const c = allNum ? ((num(x) ?? Infinity) - (num(y) ?? Infinity)) : x.localeCompare(y, undefined, {numeric: true, sensitivity: "base"});
        return asc ? c : -c;
      });
      rows.forEach((r, i) => { tbody.appendChild(r); r.classList.toggle("xtra", !th.closest(".itable").classList.contains("all") && i >= +th.closest(".itable").dataset.show); });
      applyFilter(th.closest(".itable"));
      return;
    }
    const sa = e.target.closest(".itable .showall");
    if (sa) {
      const box = sa.closest(".itable"); box.classList.toggle("all");
      sa.textContent = box.classList.contains("all") ? "show fewer" : sa.textContent.replace("show fewer", "show all");
      box.querySelectorAll("tbody tr").forEach((r, i) => r.classList.toggle("xtra", !box.classList.contains("all") && i >= +box.dataset.show));
      applyFilter(box);
    }
  });
  function applyFilter(box) {
    const q = (box.querySelector(".tfilter") || {}).value || "";
    const rows = [...box.querySelectorAll("tbody tr")];
    let n = 0;
    rows.forEach((r) => { const hit = !q || r.textContent.toLowerCase().includes(q.toLowerCase()); r.classList.toggle("nomatch", !hit); if (hit) n++; });
    const c = box.querySelector(".tcount");
    if (c) c.textContent = q ? `${n} of ${rows.length} rows match` : `${rows.length} rows`;
    if (q) box.classList.add("filtering"); else box.classList.remove("filtering");
  }
  document.addEventListener("input", (e) => { const f = e.target.closest(".tfilter"); if (f) applyFilter(f.closest(".itable")); });

  // ---------------------------------------------------------------- progress
  const logEl = $("#log"), statusEl = $("#status");
  let nLog = 0;
  async function poll() {
    try {
      const r = await fetch(api(`?since=${nLog}`));
      const d = await r.json();
      if (d.log && d.log.length && logEl) {
        logEl.textContent += (logEl.textContent ? "\n" : "") + d.log.join("\n");
        nLog += d.log.length;
        logEl.scrollTop = logEl.scrollHeight;
      }
      statusEl.textContent = d.status;
      statusEl.className = `pill ${d.status}`;
      if (d.status === "running" || d.status === "queued") setTimeout(poll, 2000);
      else if (body.dataset.status !== d.status) location.reload();
    } catch (e) { setTimeout(poll, 5000); }
  }
  if (body.dataset.status === "running" || body.dataset.status === "queued") { if (logEl) logEl.textContent = ""; poll(); }

  // ---------------------------------------------------------------- search
  const results = $("#search-results"), counts = $("#search-counts");
  function hitHtml(h) {
    const eid = h.evidence ? `<a class="eid" href="#${esc(h.evidence)}" title="open the evidence card">${esc(h.evidence)}</a>` : "";
    let extra = "";
    if (h.kind === "gene" && h.data) {
      const pw = (h.data.pathways || []).slice(0, 6).map(esc).join("; ");
      const hub = (h.data.hub_in || []).map(esc).join("; ");
      extra = `<ul>${pw ? `<li>pathways: ${pw}${(h.data.pathways || []).length > 6 ? " …" : ""}</li>` : ""}${hub ? `<li>${hub}</li>` : ""}</ul>`;
    }
    if ((h.kind === "paper" || h.kind === "dataset") && h.data && h.data.url) {
      extra = `<div class="d"><a href="${esc(h.data.url)}" target="_blank" rel="noopener">${esc(h.data.url)}</a></div>`;
    }
    if (h.data && h.data.anchor) extra = `<div class="d"><a href="#" data-part="${esc(h.data.anchor)}">open this part →</a></div>`;
    return `<div class="hit"><div class="t">${esc(h.title)}${eid}</div><div class="d">${esc(h.detail)}</div>${extra}</div>`;
  }
  async function search(q) {
    if (!results) return;
    if (!q.trim()) { results.innerHTML = ""; counts.textContent = ""; return; }
    const kinds = [...document.querySelectorAll("#kinds input:checked")].map((i) => `kind=${i.value}`).join("&");
    results.innerHTML = '<p class="thinking">searching</p>';
    const r = await fetch(api(`/search?q=${encodeURIComponent(q)}${kinds ? "&" + kinds : ""}`));
    const d = await r.json();
    if (d.error) { results.innerHTML = `<p class="flag">${esc(d.error)}</p>`; return; }
    counts.textContent = `${d.total} hit(s) for “${d.query}”`;
    const groups = Object.entries(d.groups);
    results.innerHTML = groups.length ? groups.map(([k, hits]) =>
      `<div class="group"><h3>${esc(k)} <span class="muted small">${hits.length}</span></h3>${hits.map(hitHtml).join("")}</div>`).join("")
      : '<p class="muted">Nothing matched. Try a gene symbol, a peak id, part of a pathway name, a TF, or a PMID.</p>';
  }
  if ($("#searchform")) {
    $("#searchform").addEventListener("submit", (e) => { e.preventDefault(); search($("#q").value); });
    document.querySelectorAll("#kinds input").forEach((i) => i.addEventListener("change", () => search($("#q").value)));
    if (params.get("q")) { $("#q").value = params.get("q"); showPart("search"); search(params.get("q")); }
  }

  // ---------------------------------------------------------------- pathway explorer (explore.js)
  window.peakSearch = (q) => { if ($("#q")) { $("#q").value = q; showPart("search"); search(q); } };
  document.querySelectorAll(".explore[id^='explore-']").forEach(async (box) => {
    const eid = box.id.replace("explore-", "");
    try {
      const d = await (await fetch(api(`/explore/${eid}`))).json();
      if (!window.PeakExplore) throw new Error("explore.js not loaded");
      box.innerHTML = "";
      PeakExplore.mount(box, d, {title: `${jobId}_${eid}`, onGene: (g) => window.peakSearch(g),
        geneInfo: async (g) => (await fetch(api(`/gene/${encodeURIComponent(g)}`))).json()});
    } catch (e) { box.innerHTML = `<p class="muted">Interactive view unavailable: ${esc(e.message || e)}</p>`; }
  });

  // ---------------------------------------------------------------- chat
  const messages = $("#messages"), input = $("#chat-input"), chatForm = $("#chatform");
  let chatLoaded = false;
  function msgHtml(m) {
    if (m.role === "user") return `<div class="msg user">${esc(m.text)}</div>`;
    const tools = (m.tool_calls || []).map((c) => `<span title='${esc(JSON.stringify(c.args))}'>${esc(c.tool.replace(/_/g, " "))}${c.ok ? "" : " ✗"}</span>`).join("");
    const err = m.error ? `<div class="err">${esc(m.error)}</div>` : "";
    const unv = (m.unverified_citations || []).length ? `<div class="err small">Unverified citations: ${esc(m.unverified_citations.join(", "))}</div>` : "";
    return `<div class="msg assistant">${tools ? `<div class="tools">${tools}</div>` : ""}${m.html || esc(m.text)}${err}${unv}</div>`;
  }
  async function loadChat() {
    if (chatLoaded || !messages) return;
    chatLoaded = true;
    try {
      const d = await (await fetch(api("/chat"))).json();
      (d.messages || []).forEach((m) => messages.insertAdjacentHTML("beforeend", msgHtml(m)));
      messages.scrollTop = messages.scrollHeight;
    } catch (e) {}
  }
  async function send(text) {
    text = (text || input.value).trim();
    if (!text) return;
    input.value = "";
    messages.insertAdjacentHTML("beforeend", msgHtml({role: "user", text}));
    const wait = document.createElement("div");
    wait.className = "msg assistant thinking";
    wait.textContent = "Reading the run";
    messages.appendChild(wait);
    messages.scrollTop = messages.scrollHeight;
    chatForm.querySelector("button.primary").disabled = true;
    try {
      const r = await fetch(api("/chat"), {method: "POST", headers: {"Content-Type": "application/json"}, body: JSON.stringify({message: text})});
      const d = await r.json();
      wait.remove();
      messages.insertAdjacentHTML("beforeend", d.error && !d.text ? `<div class="msg assistant"><div class="err">${esc(d.error)}</div></div>` : msgHtml(d));
    } catch (e) { wait.textContent = `Request failed: ${e}`; }
    chatForm.querySelector("button.primary").disabled = false;
    messages.scrollTop = messages.scrollHeight;
  }
  if (chatForm) {
    chatForm.addEventListener("submit", (e) => { e.preventDefault(); send(); });
    input.addEventListener("keydown", (e) => { if (e.key === "Enter" && !e.shiftKey) { e.preventDefault(); send(); } });
    $("#chat-clear").addEventListener("click", async () => {
      await fetch(api("/chat"), {method: "DELETE"});
      messages.querySelectorAll(".msg").forEach((m) => m.remove());
    });
    messages.addEventListener("click", (e) => {
      const b = e.target.closest(".suggest button");
      if (b) send(b.textContent);
    });
  }
})();
