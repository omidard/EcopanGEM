// EcopanGEM Flux Analysis Studio — FVA, Dynamic FBA & Multi-model tabs.
// (The Explore/Compare tab is handled by fba_ui.js.) All client-side.
import { runFBA, runPFBA, runFVA, runDFBA, exchangeReport, listExchanges } from './fba_engine.js';

const $ = (id) => document.getElementById(id);
const esc = (s) => String(s == null ? '' : s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const fmt = (x, n = 3) => (x == null || isNaN(x)) ? '—' : Number(x).toFixed(n);
let PRESETS = null;

const PHYLO_COLORS = { A: '#4e79a7', B1: '#59a14f', B2: '#e15759', C: '#f28e2b', D: '#76b7b2',
  E: '#edc948', F: '#b07aa1', G: '#ff9da7', clade: '#9c755f', Unknown: '#bab0ac', '': '#bab0ac' };
const phyloColor = (p) => PHYLO_COLORS[p] || '#8888aa';

// ── shared data helpers ───────────────────────────────────────────────────────
async function presets() { if (!PRESETS) PRESETS = await (await fetch('fba/media_presets.json')).json(); return PRESETS; }
function fillMedia(sel, descEl) {
  sel.innerHTML = '';
  for (const [k, v] of Object.entries(PRESETS)) { const o = document.createElement('option'); o.value = k; o.textContent = v.label; sel.appendChild(o); }
  const upd = () => { if (descEl) { const p = PRESETS[sel.value]; descEl.textContent = `${p.desc} · ${Object.keys(p.bounds).length} exchanges · ${p.source}`; } };
  sel.addEventListener('change', upd); upd();
}
const metaByFile = new Map();
function meta(gemFile) { if (!metaByFile.size && window.gemMetadata) window.gemMetadata.forEach(m => metaByFile.set(m.gem_file, m)); return metaByFile.get(gemFile) || {}; }

const modelCache = new Map();
async function loadModel(gemFile) {
  if (modelCache.has(gemFile)) return modelCache.get(gemFile);
  const bn = window.gemBatchMap && window.gemBatchMap[gemFile];
  if (!bn) throw new Error(`Unknown model "${gemFile}".`);
  const resp = await fetch(window.BATCH_URL_BASE + String(bn).padStart(2, '0') + '.zip');
  if (!resp.ok) throw new Error('Batch download failed.');
  const zip = await window.JSZip.loadAsync(await resp.blob());
  const entry = zip.file(gemFile); if (!entry) throw new Error('Model not in batch.');
  const model = JSON.parse(await entry.async('text')); modelCache.set(gemFile, model); return model;
}

// reusable model combobox. onPick(gemFile). If multi, input keeps focus.
function makeCombo(input, menu, onPick) {
  let items = [], active = -1;
  const render = () => {
    const q = input.value.trim().toLowerCase();
    const md = window.gemMetadata || [];
    items = (q ? md.filter(m => (m.gem_file || '').toLowerCase().includes(q) || (m.genome_name || '').toLowerCase().includes(q) || (m.strain || '').toLowerCase().includes(q)) : md).slice(0, 40);
    menu.innerHTML = items.length ? items.map((m, i) => `<div class="fba-combo-item${i === active ? ' active' : ''}" data-i="${i}"><div class="nm">${esc(m.genome_name || m.strain || m.gem_file)}</div><div class="id">${esc(m.gem_file)}</div><div class="meta">${m.phylogroup ? 'Phylogroup ' + esc(m.phylogroup) : ''}${m.MLST ? ' · ST' + esc(m.MLST) : ''}${m.isolation_source ? ' · ' + esc(m.isolation_source) : ''}</div></div>`).join('') : `<div class="fba-combo-empty">No models match.</div>`;
    menu.classList.add('show');
  };
  input.addEventListener('focus', render); input.addEventListener('input', () => { active = -1; render(); });
  input.addEventListener('keydown', (e) => {
    if (!menu.classList.contains('show')) return;
    if (e.key === 'ArrowDown') { active = Math.min(active + 1, items.length - 1); render(); e.preventDefault(); }
    else if (e.key === 'ArrowUp') { active = Math.max(active - 1, 0); render(); e.preventDefault(); }
    else if (e.key === 'Enter') { const it = items[active] || items[0]; if (it) { onPick(it.gem_file); e.preventDefault(); } }
    else if (e.key === 'Escape') menu.classList.remove('show');
  });
  menu.addEventListener('mousedown', (e) => { const it = e.target.closest('.fba-combo-item'); if (!it) return; e.preventDefault(); onPick(items[+it.dataset.i].gem_file); });
  input.addEventListener('blur', () => setTimeout(() => menu.classList.remove('show'), 150));
}

function modelCardHTML(gemFile, model) {
  const m = meta(gemFile);
  return `<div class="mc-name">${esc(m.genome_name || m.strain || gemFile)}</div><div class="mc-id">${esc(gemFile)}</div>
    <div class="mc-tags"><span class="fba-tag">${model.reactions.length} rxns</span><span class="fba-tag">${model.metabolites.length} mets</span>
    ${m.phylogroup ? `<span class="fba-tag">Phylogroup ${esc(m.phylogroup)}</span>` : ''}${m.MLST ? `<span class="fba-tag">ST ${esc(m.MLST)}</span>` : ''}${m.isolation_source ? `<span class="fba-tag">${esc(m.isolation_source)}</span>` : ''}</div>`;
}
function setStatus(id, msg, kind) { const e = $(id); e.textContent = msg; e.className = 'fba-status ' + (kind || ''); }
function prog(wrapId, barId, frac) { const w = $(wrapId); w.style.display = frac == null ? 'none' : 'block'; if (frac != null) $(barId).style.width = Math.round(frac * 100) + '%'; }
function saveCSV(csv, base) { const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([csv], { type: 'text/csv' })); a.download = base + '.csv'; a.click(); }
const bio = (id) => id.replace(/^EX_/, '').replace(/_e$/, '');

// ── Tabs ──────────────────────────────────────────────────────────────────────
function initTabs() {
  document.querySelectorAll('#studio-tabs .studio-tab').forEach(btn => btn.addEventListener('click', () => {
    const tab = btn.dataset.tab;
    document.querySelectorAll('#studio-tabs .studio-tab').forEach(b => b.classList.toggle('active', b === btn));
    document.querySelectorAll('.studio-panel').forEach(p => p.classList.toggle('active', p.id === 'panel-' + tab));
    // resize any Plotly plots in the now-visible panel
    setTimeout(() => document.querySelectorAll('#panel-' + tab + ' .js-plotly-plot').forEach(p => window.Plotly && window.Plotly.Plots.resize(p)), 30);
  }));
}

// ── FVA tab ───────────────────────────────────────────────────────────────────
const fvaState = { model: null, file: null };
function initFVA() {
  fillMedia($('fva-media'), $('fva-media-desc'));
  makeCombo($('fva-model-input'), $('fva-model-menu'), async (gemFile) => {
    $('fva-model-input').value = '';
    const mc = $('fva-modelcard'); mc.style.display = 'block'; mc.innerHTML = '<span class="fba-hint-inline">Loading…</span>';
    try { const model = await loadModel(gemFile); fvaState.model = model; fvaState.file = gemFile; mc.innerHTML = modelCardHTML(gemFile, model); }
    catch (e) { mc.innerHTML = `<span style="color:#c0392b">${esc(e.message)}</span>`; }
  });
  $('fva-frac').addEventListener('input', () => $('fva-frac-val').textContent = (+$('fva-frac').value).toFixed(2));
  $('fva-run').addEventListener('click', runFVAtab);
}
async function runFVAtab() {
  if (!fvaState.model) return setStatus('fva-status', 'Choose a model first.', 'err');
  const media = PRESETS[$('fva-media').value].bounds;
  const fraction = +$('fva-frac').value;
  const exIds = listExchanges(fvaState.model).map(e => e.id);
  $('fva-run').disabled = true; $('fva-results').style.display = 'none';
  setStatus('fva-status', `Running FVA on ${exIds.length} exchanges…`, 'busy'); prog('fva-prog', 'fva-prog-bar', 0);
  try {
    const t0 = performance.now();
    const res = await runFVA(fvaState.model, media, exIds, { fraction, onProgress: (d, t) => prog('fva-prog', 'fva-prog-bar', d / t) });
    prog('fva-prog', 'fva-prog-bar', null);
    if (!res.optimal) { setStatus('fva-status', 'Model does not grow on this medium.', 'err'); $('fva-run').disabled = false; return; }
    renderFVA(res, fvaState.model);
    setStatus('fva-status', `Done — ${exIds.length} exchanges in ${((performance.now() - t0) / 1000).toFixed(1)} s. Optimum growth ${fmt(res.z)} h⁻¹.`, 'ok');
  } catch (e) { setStatus('fva-status', 'Error: ' + e.message, 'err'); console.error(e); prog('fva-prog', 'fva-prog-bar', null); }
  finally { $('fva-run').disabled = false; }
}
function renderFVA(res, model) {
  const nameById = {}; model.reactions.forEach(r => nameById[r.id] = r.name || '');
  // keep exchanges with non-trivial range or nonzero, sort by span
  let rows = Object.entries(res.ranges).map(([id, r]) => ({ id, name: nameById[id], min: r.min, max: r.max, span: r.max - r.min }))
    .filter(r => Math.abs(r.max) > 1e-6 || Math.abs(r.min) > 1e-6);
  rows.sort((a, b) => b.span - a.span);
  const top = rows.slice(0, 30).reverse();
  $('fva-results').style.display = 'block';
  window.Plotly.newPlot('fva-plot', [{
    type: 'scatter', mode: 'markers', x: top.map(r => r.min), y: top.map(bioName), name: 'min',
    marker: { color: '#2c6fbb', size: 8, symbol: 'line-ns-open' }, hovertemplate: '%{y}<br>min %{x:.3f}<extra></extra>'
  }, {
    type: 'scatter', mode: 'markers', x: top.map(r => r.max), y: top.map(bioName), name: 'max',
    marker: { color: '#c0392b', size: 8 }, hovertemplate: '%{y}<br>max %{x:.3f}<extra></extra>'
  }, {
    type: 'scatter', mode: 'lines', x: top.flatMap(r => [r.min, r.max, null]), y: top.flatMap(r => [bioName(r), bioName(r), null]),
    line: { color: '#9db8d6', width: 6 }, hoverinfo: 'skip', showlegend: false
  }], {
    margin: { l: 110, r: 20, t: 10, b: 40 }, height: Math.max(420, top.length * 16),
    xaxis: { title: 'Flux range (mmol·gDW⁻¹·h⁻¹) — biomass ≥ ' + res.fraction + '× optimum', zeroline: true, zerolinecolor: '#ccc' },
    yaxis: { automargin: true }, legend: { orientation: 'h', y: 1.05 }, font: { size: 11 }
  }, { responsive: true, displaylogo: false });
  function bioName(r) { return bio(r.id); }
  // table
  const tr = rows.map(r => `<tr><td><code>${esc(r.id)}</code></td><td>${esc(r.name)}</td><td class="num">${fmt(r.min)}</td><td class="num">${fmt(r.max)}</td><td class="num">${fmt(r.span)}</td></tr>`).join('');
  $('fva-table').innerHTML = `<div class="fba-tablewrap" style="max-height:320px"><table class="fba-flux"><thead><tr><th>Exchange</th><th>Name</th><th>Min</th><th>Max</th><th>Span</th></tr></thead><tbody>${tr}</tbody></table></div>
    <button class="btn btn-sm btn-outline-secondary mt-2" id="fva-csv">⬇ Download FVA (CSV)</button>`;
  $('fva-csv').addEventListener('click', () => { let c = 'exchange_id,name,min,max,span\n'; rows.forEach(r => c += `${r.id},"${(r.name || '').replace(/"/g, '""')}",${r.min},${r.max},${r.span}\n`); saveCSV(c, `FVA_${fvaState.file.replace(/\.json.*/, '')}_${$('fva-media').value}`); });
}

// ── Dynamic FBA tab ────────────────────────────────────────────────────────────
const dfbaState = { model: null, file: null, series: null };
function initDFBA() {
  fillMedia($('dfba-media'), null);
  makeCombo($('dfba-model-input'), $('dfba-model-menu'), async (gemFile) => {
    $('dfba-model-input').value = '';
    const mc = $('dfba-modelcard'); mc.style.display = 'block'; mc.innerHTML = '<span class="fba-hint-inline">Loading…</span>';
    try {
      const model = await loadModel(gemFile); dfbaState.model = model; dfbaState.file = gemFile; mc.innerHTML = modelCardHTML(gemFile, model);
      // populate substrate options from carbon exchanges
      const sub = $('dfba-substrate'); const cur = sub.value;
      const exs = listExchanges(model).filter(e => /glc|glucose|fru|gal|lac|ac_e|succ|glyc|sucr|malt|xyl|arab/i.test(e.id));
      sub.innerHTML = exs.map(e => `<option value="${e.id}">${esc((e.name || e.id).replace(/ exchange$/i, ''))} (${e.id})</option>`).join('');
      if ([...sub.options].some(o => o.value === 'EX_glc__D_e')) sub.value = 'EX_glc__D_e'; else if (cur) sub.value = cur;
    } catch (e) { mc.innerHTML = `<span style="color:#c0392b">${esc(e.message)}</span>`; }
  });
  $('dfba-run').addEventListener('click', runDFBAtab);
}
async function runDFBAtab() {
  if (!dfbaState.model) return setStatus('dfba-status', 'Choose a model first.', 'err');
  const media = PRESETS[$('dfba-media').value].bounds;
  const substrateEx = $('dfba-substrate').value;
  const opts = {
    substrateEx, substrate0: +$('dfba-s0').value, biomass0: +$('dfba-x0').value,
    vmax: +$('dfba-vmax').value, km: +$('dfba-km').value, dt: +$('dfba-dt').value, tmax: +$('dfba-tmax').value,
    trackEx: ['EX_ac_e', 'EX_for_e', 'EX_etoh_e', 'EX_lac__D_e', 'EX_succ_e'].filter(e => dfbaState.model.reactions.some(r => r.id === e)),
    onProgress: (i, n) => prog('dfba-prog', 'dfba-prog-bar', i / n),
  };
  $('dfba-run').disabled = true; $('dfba-results').style.display = 'none';
  setStatus('dfba-status', 'Simulating batch culture…', 'busy'); prog('dfba-prog', 'dfba-prog-bar', 0);
  try {
    const t0 = performance.now();
    const s = await runDFBA(dfbaState.model, media, opts);
    prog('dfba-prog', 'dfba-prog-bar', null);
    dfbaState.series = s;
    renderDFBA(s, substrateEx, dfbaState.model);
    const finalX = s.biomass[s.biomass.length - 1];
    setStatus('dfba-status', `Done — ${s.t.length} steps in ${((performance.now() - t0) / 1000).toFixed(1)} s. Final biomass ${fmt(finalX)} gDW/L.`, 'ok');
  } catch (e) { setStatus('dfba-status', 'Error: ' + e.message, 'err'); console.error(e); prog('dfba-prog', 'dfba-prog-bar', null); }
  finally { $('dfba-run').disabled = false; }
}
function renderDFBA(s, substrateEx, model) {
  const nameById = {}; model.reactions.forEach(r => nameById[r.id] = r.name || '');
  $('dfba-results').style.display = 'block';
  const traces = [{ x: s.t, y: s.biomass, name: 'Biomass (gDW/L)', yaxis: 'y2', line: { color: '#1a7f4b', width: 3 } }];
  const palette = ['#2c6fbb', '#c0392b', '#e08a1e', '#7d3c98', '#16a085'];
  let ci = 0;
  for (const [ex, arr] of Object.entries(s.conc)) {
    const isSub = ex === substrateEx;
    traces.push({ x: s.t, y: arr, name: (isSub ? '▸ ' : '') + bio(ex) + ' (mM)', line: { color: isSub ? '#333' : palette[ci++ % palette.length], width: isSub ? 3 : 2, dash: isSub ? 'solid' : 'solid' } });
  }
  window.Plotly.newPlot('dfba-plot', traces, {
    margin: { l: 55, r: 55, t: 10, b: 45 }, height: 460,
    xaxis: { title: 'Time (h)' },
    yaxis: { title: 'Concentration (mM)', rangemode: 'tozero' },
    yaxis2: { title: 'Biomass (gDW/L)', overlaying: 'y', side: 'right', rangemode: 'tozero', showgrid: false },
    legend: { orientation: 'h', y: 1.08 }, font: { size: 11 }, hovermode: 'x unified',
  }, { responsive: true, displaylogo: false });
  $('dfba-note').textContent = `Substrate ▸ ${bio(substrateEx)} depletes via Michaelis–Menten uptake; biomass (right axis, green) grows until the substrate runs out; fermentation products accumulate.`;
  $('dfba-csv').onclick = () => {
    const keys = Object.keys(s.conc);
    let csv = 'time_h,biomass_gDW_L,' + keys.map(k => bio(k) + '_mM').join(',') + '\n';
    for (let i = 0; i < s.t.length; i++) csv += `${s.t[i]},${s.biomass[i]},` + keys.map(k => s.conc[k][i]).join(',') + '\n';
    saveCSV(csv, `dFBA_${dfbaState.file.replace(/\.json.*/, '')}_${$('dfba-media').value}`);
  };
}

// ── Multi-model tab ────────────────────────────────────────────────────────────
const mm = { selected: [], results: null };
function initMulti() {
  fillMedia($('mm-media'), $('mm-media-desc'));
  makeCombo($('mm-model-input'), $('mm-model-menu'), (gemFile) => { addModel(gemFile); $('mm-model-input').value = ''; });
  document.querySelectorAll('.mm-quick button[data-add]').forEach(b => b.addEventListener('click', () => quickAdd(b.dataset.add)));
  $('mm-clear').addEventListener('click', () => { mm.selected = []; renderChips(); });
  $('mm-run').addEventListener('click', runMulti);
}
function addModel(gemFile) { if (!mm.selected.includes(gemFile) && window.gemBatchMap[gemFile]) mm.selected.push(gemFile); renderChips(); }
function quickAdd(kind) {
  const md = window.gemMetadata || [];
  if (kind === 'perphylo') {
    const byP = {}; md.forEach(m => { const p = m.phylogroup || 'Unknown'; (byP[p] = byP[p] || []).push(m.gem_file); });
    Object.values(byP).forEach(list => { for (let i = 0; i < Math.min(4, list.length); i++) addModel(list[Math.floor(i * list.length / 4)]); });
  } else {
    const n = +kind, pool = md.map(m => m.gem_file);
    for (let k = 0; k < n && mm.selected.length < 120; k++) addModel(pool[Math.floor(seededRand(mm.selected.length + k) * pool.length)]);
  }
}
let _seed = 12345; function seededRand(i) { _seed = (1103515245 * (_seed + i) + 12345) & 0x7fffffff; return _seed / 0x7fffffff; }
function renderChips() {
  $('mm-count').textContent = `${mm.selected.length} models selected`;
  $('mm-chips').innerHTML = mm.selected.map(f => { const m = meta(f); return `<span class="mm-chip" title="${esc(f)}">${esc((m.genome_name || m.strain || f).slice(0, 26))} <span class="x" data-f="${esc(f)}">✕</span></span>`; }).join('');
  $('mm-chips').querySelectorAll('.x').forEach(x => x.addEventListener('click', () => { mm.selected = mm.selected.filter(f => f !== x.dataset.f); renderChips(); }));
}
async function runMulti() {
  if (mm.selected.length < 2) return setStatus('mm-status', 'Add at least 2 models.', 'err');
  const media = PRESETS[$('mm-media').value].bounds;
  const meth = document.querySelector('input[name="mm-method"]:checked').value;
  $('mm-run').disabled = true; $('mm-results').style.display = 'none';
  setStatus('mm-status', `Running ${meth.toUpperCase()} on ${mm.selected.length} models…`, 'busy'); prog('mm-prog', 'mm-prog-bar', 0);
  try {
    const t0 = performance.now();
    const rows = [];
    for (let i = 0; i < mm.selected.length; i++) {
      const f = mm.selected[i];
      try {
        const model = await loadModel(f);
        const fba = await runFBA(model, media);
        let res = fba;
        if (meth === 'pfba' && fba.optimal && fba.growth > 1e-9) res = await runPFBA(model, media, fba);
        const ex = {}; for (const [id, v] of Object.entries(res.fluxes)) if (id.startsWith('EX_') && Math.abs(v) > 1e-6) ex[id] = v;
        rows.push({ file: f, meta: meta(f), growth: fba.optimal ? fba.growth : 0, ex });
      } catch (e) { rows.push({ file: f, meta: meta(f), growth: 0, ex: {}, error: e.message }); }
      setStatus('mm-status', `Solved ${i + 1}/${mm.selected.length}…`, 'busy'); prog('mm-prog', 'mm-prog-bar', (i + 1) / mm.selected.length);
    }
    prog('mm-prog', 'mm-prog-bar', null);
    mm.results = rows;
    renderMulti(rows);
    setStatus('mm-status', `Done — ${rows.length} models in ${((performance.now() - t0) / 1000).toFixed(1)} s.`, 'ok');
  } catch (e) { setStatus('mm-status', 'Error: ' + e.message, 'err'); console.error(e); prog('mm-prog', 'mm-prog-bar', null); }
  finally { $('mm-run').disabled = false; }
}
function renderMulti(rows) {
  $('mm-results').style.display = 'block';
  const feas = rows.filter(r => r.growth > 1e-9);
  const meanG = feas.length ? feas.reduce((s, r) => s + r.growth, 0) / feas.length : 0;
  $('mm-kpis').innerHTML =
    `<div class="fba-kpi"><div class="v">${rows.length}</div><div class="l">Models run</div></div>
     <div class="fba-kpi"><div class="v">${feas.length}</div><div class="l">Feasible (grow)</div></div>
     <div class="fba-kpi"><div class="v">${fmt(meanG)}</div><div class="l">Mean growth (h⁻¹)</div></div>
     <div class="fba-kpi"><div class="v">${new Set(rows.map(r => r.meta.phylogroup || 'Unknown')).size}</div><div class="l">Phylogroups</div></div>`;

  growthByPhylo(rows);
  buildAxisSelectors(rows);
  drawScatter(rows);
  drawPCA(rows);
  drawHeatmap(rows);
  $('mm-csv').onclick = () => downloadMultiCSV(rows);
  $('mm-detail').innerHTML = '';
}
function growthByPhylo(rows) {
  const groups = {};
  rows.forEach(r => { const p = r.meta.phylogroup || 'Unknown'; (groups[p] = groups[p] || []).push(r); });
  const traces = Object.entries(groups).map(([p, rs]) => ({
    type: 'box', boxpoints: 'all', jitter: 0.5, pointpos: 0, name: p, y: rs.map(r => r.growth),
    text: rs.map(r => r.meta.genome_name || r.file), customdata: rs.map(r => r.file),
    marker: { color: phyloColor(p), size: 6 }, line: { color: phyloColor(p) },
    hovertemplate: '%{text}<br>growth %{y:.3f}<extra>' + p + '</extra>',
  }));
  window.Plotly.newPlot('mm-plot-growth', traces, { margin: { l: 45, r: 15, t: 10, b: 35 }, height: 340, yaxis: { title: 'Growth (h⁻¹)', rangemode: 'tozero' }, showlegend: false, font: { size: 11 } }, { responsive: true, displaylogo: false })
    .then(gd => gd.on('plotly_click', ev => showDetail(ev.points[0].customdata)));
}
function commonExchanges(rows, minCount = 2) {
  const count = {}; rows.forEach(r => Object.keys(r.ex).forEach(id => count[id] = (count[id] || 0) + 1));
  return Object.entries(count).filter(([, c]) => c >= minCount).map(([id]) => id).sort();
}
function buildAxisSelectors(rows) {
  const exs = commonExchanges(rows);
  const opts = [`<option value="__growth">Growth rate (h⁻¹)</option>`].concat(exs.map(id => `<option value="${id}">${esc(bio(id))} flux</option>`)).join('');
  $('mm-x').innerHTML = opts; $('mm-y').innerHTML = opts;
  $('mm-x').value = '__growth';
  $('mm-y').value = exs.includes('EX_ac_e') ? 'EX_ac_e' : (exs[0] || '__growth');
  $('mm-x').onchange = $('mm-y').onchange = () => drawScatter(rows);
}
function axisVal(r, key) { return key === '__growth' ? r.growth : (r.ex[key] || 0); }
function axisLabel(key) { return key === '__growth' ? 'Growth rate (h⁻¹)' : bio(key) + ' flux'; }
function drawScatter(rows) {
  const xk = $('mm-x').value, yk = $('mm-y').value;
  const groups = {};
  rows.forEach(r => { const p = r.meta.phylogroup || 'Unknown'; (groups[p] = groups[p] || []).push(r); });
  const traces = Object.entries(groups).map(([p, rs]) => ({
    type: 'scatter', mode: 'markers', name: p, x: rs.map(r => axisVal(r, xk)), y: rs.map(r => axisVal(r, yk)),
    text: rs.map(r => r.meta.genome_name || r.file), customdata: rs.map(r => r.file),
    marker: { color: phyloColor(p), size: 9, line: { color: '#fff', width: 0.5 } },
    hovertemplate: '%{text}<br>' + esc(axisLabel(xk)) + ' %{x:.3f}<br>' + esc(axisLabel(yk)) + ' %{y:.3f}<extra>' + p + '</extra>',
  }));
  window.Plotly.newPlot('mm-plot-scatter', traces, { margin: { l: 50, r: 15, t: 10, b: 45 }, height: 340, xaxis: { title: axisLabel(xk) }, yaxis: { title: axisLabel(yk) }, legend: { font: { size: 9 } }, font: { size: 11 } }, { responsive: true, displaylogo: false })
    .then(gd => gd.on('plotly_click', ev => showDetail(ev.points[0].customdata)));
}

// PCA via Gram-matrix (n×n) eigendecomposition (Jacobi)
function drawPCA(rows) {
  const feats = commonExchanges(rows, 2);
  if (feats.length < 2 || rows.length < 3) { $('mm-plot-pca').innerHTML = '<div class="fba-hint-inline" style="padding:1rem">Need ≥3 models and ≥2 shared exchanges for PCA.</div>'; return; }
  const n = rows.length;
  // matrix X (n×p), center columns
  let X = rows.map(r => feats.map(f => r.ex[f] || 0));
  const mean = feats.map((_, j) => X.reduce((s, row) => s + row[j], 0) / n);
  X = X.map(row => row.map((v, j) => v - mean[j]));
  // Gram G = X Xᵀ (n×n)
  const G = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, k) => { let s = 0; for (let j = 0; j < feats.length; j++) s += X[i][j] * X[k][j]; return s; }));
  const { values, vectors } = jacobiEigen(G);
  const order = values.map((v, i) => [v, i]).sort((a, b) => b[0] - a[0]);
  const tot = values.reduce((s, v) => s + Math.max(0, v), 0) || 1;
  const pc = (rank) => { const [lam, idx] = order[rank]; const sq = Math.sqrt(Math.max(0, lam)); return { scores: vectors.map(row => row[idx] * sq), ev: Math.max(0, lam) / tot }; };
  const p1 = pc(0), p2 = pc(1);
  const groups = {};
  rows.forEach((r, i) => { const p = r.meta.phylogroup || 'Unknown'; (groups[p] = groups[p] || []).push(i); });
  const traces = Object.entries(groups).map(([p, idxs]) => ({
    type: 'scatter', mode: 'markers', name: p, x: idxs.map(i => p1.scores[i]), y: idxs.map(i => p2.scores[i]),
    text: idxs.map(i => rows[i].meta.genome_name || rows[i].file), customdata: idxs.map(i => rows[i].file),
    marker: { color: phyloColor(p), size: 9, line: { color: '#fff', width: 0.5 } },
    hovertemplate: '%{text}<extra>' + p + '</extra>',
  }));
  window.Plotly.newPlot('mm-plot-pca', traces, { margin: { l: 50, r: 15, t: 10, b: 45 }, height: 340,
    xaxis: { title: `PC1 (${(p1.ev * 100).toFixed(0)}%)` }, yaxis: { title: `PC2 (${(p2.ev * 100).toFixed(0)}%)` }, legend: { font: { size: 9 } }, font: { size: 11 } }, { responsive: true, displaylogo: false })
    .then(gd => gd.on('plotly_click', ev => showDetail(ev.points[0].customdata)));
}
// Jacobi eigenvalue algorithm for symmetric matrices
function jacobiEigen(A) {
  const n = A.length; const a = A.map(r => r.slice());
  const V = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => i === j ? 1 : 0));
  for (let sweep = 0; sweep < 100; sweep++) {
    let off = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += a[i][j] * a[i][j];
    if (off < 1e-12) break;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
      if (Math.abs(a[p][q]) < 1e-14) continue;
      const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
      const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < n; k++) { const akp = a[k][p], akq = a[k][q]; a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq; }
      for (let k = 0; k < n; k++) { const apk = a[p][k], aqk = a[q][k]; a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk; }
      for (let k = 0; k < n; k++) { const vkp = V[k][p], vkq = V[k][q]; V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq; }
    }
  }
  return { values: a.map((r, i) => r[i]), vectors: V };
}
function drawHeatmap(rows) {
  const feats = commonExchanges(rows, 2);
  // most variable exchanges
  const varOf = (f) => { const vals = rows.map(r => r.ex[f] || 0); const m = vals.reduce((s, v) => s + v, 0) / vals.length; return vals.reduce((s, v) => s + (v - m) ** 2, 0); };
  const top = feats.map(f => [f, varOf(f)]).sort((a, b) => b[1] - a[1]).slice(0, 20).map(x => x[0]);
  const order = rows.map((r, i) => [r.growth, i]).sort((a, b) => b[0] - a[0]).map(x => x[1]);
  const z = order.map(i => top.map(f => rows[i].ex[f] || 0));
  window.Plotly.newPlot('mm-plot-heat', [{
    type: 'heatmap', z, x: top.map(bio), y: order.map(i => (rows[i].meta.genome_name || rows[i].file).slice(0, 22)),
    colorscale: 'RdBu', reversescale: true, zmid: 0, colorbar: { title: 'flux', thickness: 10 },
    hovertemplate: '%{y}<br>%{x}: %{z:.2f}<extra></extra>',
  }], { margin: { l: 150, r: 10, t: 10, b: 70 }, height: Math.max(300, order.length * 16 + 90), xaxis: { tickangle: -45, tickfont: { size: 9 } }, yaxis: { tickfont: { size: 9 } }, font: { size: 11 } }, { responsive: true, displaylogo: false });
}
function showDetail(file) {
  if (!file) return;
  const r = mm.results.find(x => x.file === file); if (!r) return;
  const m = r.meta; const rep = { sec: Object.entries(r.ex).filter(([, v]) => v > 1e-6).sort((a, b) => b[1] - a[1]).slice(0, 6) };
  $('mm-detail').innerHTML = `<div class="fba-note" style="background:var(--accent)">
    <strong>${esc(m.genome_name || m.strain || file)}</strong> <code>${esc(file)}</code>
    ${m.phylogroup ? '· Phylogroup ' + esc(m.phylogroup) : ''} · growth <strong>${fmt(r.growth)}</strong> h⁻¹
    · top secretions: ${rep.sec.map(([id, v]) => `${bio(id)} ${v.toFixed(1)}`).join(', ') || '—'}
    &nbsp;<a href="#" id="mm-open-explore">open in Explore →</a></div>`;
  $('mm-open-explore').onclick = (e) => { e.preventDefault(); if (window.fbaSetModel) { document.querySelector('#studio-tabs .studio-tab[data-tab="explore"]').click(); window.fbaSetModel('a', file); } };
}
function downloadMultiCSV(rows) {
  const feats = commonExchanges(rows, 1);
  let csv = 'gem_file,genome_name,phylogroup,growth_h-1,' + feats.map(bio).join(',') + '\n';
  rows.forEach(r => { csv += `${r.file},"${(r.meta.genome_name || '').replace(/"/g, '""')}",${r.meta.phylogroup || ''},${r.growth},` + feats.map(f => r.ex[f] || 0).join(',') + '\n'; });
  saveCSV(csv, `multimodel_${$('mm-media').value}`);
}

// ── init ──────────────────────────────────────────────────────────────────────
async function init() {
  await presets();
  initTabs(); initFVA(); initDFBA(); initMulti();
  // URL param handoff: ?tab=&model=&slot=&models=
  const q = new URLSearchParams(location.search);
  if (q.get('tab')) { const btn = document.querySelector(`#studio-tabs .studio-tab[data-tab="${q.get('tab')}"]`); if (btn) btn.click(); }
  if (q.get('models')) { q.get('models').split(',').filter(Boolean).forEach(addModel); document.querySelector('#studio-tabs .studio-tab[data-tab="multi"]').click(); }
}
if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
