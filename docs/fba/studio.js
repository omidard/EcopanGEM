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

// ── Statistics (validated vs scipy) ───────────────────────────────────────────
function erf(x) { const t = 1 / (1 + 0.3275911 * Math.abs(x)); const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x); return x >= 0 ? y : -y; }
function normalCdf(z) { return 0.5 * (1 + erf(z / Math.SQRT2)); }
function mannWhitneyU(a, b) {
  const n1 = a.length, n2 = b.length, N = n1 + n2;
  if (!n1 || !n2) return { U: NaN, p: 1 };
  const all = a.map(v => [v, 0]).concat(b.map(v => [v, 1])).sort((x, y) => x[0] - y[0]);
  const ranks = new Array(N); let i = 0; const tie = [];
  while (i < N) { let j = i; while (j < N - 1 && all[j + 1][0] === all[i][0]) j++; const r = (i + j + 2) / 2; for (let k = i; k <= j; k++) ranks[k] = r; if (j > i) tie.push(j - i + 1); i = j + 1; }
  let R1 = 0; for (let k = 0; k < N; k++) if (all[k][1] === 0) R1 += ranks[k];
  const U1 = R1 - n1 * (n1 + 1) / 2, U = Math.min(U1, n1 * n2 - U1), mu = n1 * n2 / 2;
  const tieTerm = tie.reduce((s, t) => s + (t * t * t - t), 0);
  const sigma = Math.sqrt(n1 * n2 / 12 * ((N + 1) - tieTerm / (N * (N - 1))));
  if (!(sigma > 0)) return { U, p: 1 };
  return { U, p: Math.min(1, 2 * normalCdf((U - mu + 0.5) / sigma)) };
}
function benjaminiHochberg(ps) {
  const m = ps.length, idx = ps.map((p, i) => [p, i]).sort((a, b) => a[0] - b[0]), q = new Array(m); let prev = 1;
  for (let k = m - 1; k >= 0; k--) { const [p, i] = idx[k]; const v = Math.min(prev, p * m / (k + 1)); q[i] = v; prev = v; }
  return q;
}

// ── Cohort comparison tab ──────────────────────────────────────────────────────
const COHORT_FIELDS = [['phylogroup', 'Phylogroup'], ['MLST', 'MLST (ST)'], ['pathovar', 'Pathotype'], ['isolation_source', 'Isolation source'], ['host_name', 'Host'], ['isolation_country', 'Country'], ['mash_cluster', 'MASH cluster'], ['serovar', 'Serovar'], ['oxygen_requirement', 'Oxygen requirement'], ['disease', 'Disease']];
const cohort = { a: { field: 'phylogroup', values: new Set() }, b: { field: 'phylogroup', values: new Set(), complement: false }, results: null };
const clean = (v) => { const s = String(v == null ? '' : v).trim(); return (s === '' || s === 'nan' || s === 'None' || s === '-1') ? '' : s; };

function initCohort() {
  fillMedia($('cohort-media'), null);
  ['a', 'b'].forEach(c => setupCohortBuilder(c));
  $('cohort-b-complement').addEventListener('change', () => { $('cohort-b-manual').style.display = $('cohort-b-complement').checked ? 'none' : 'block'; updateCohortCount('b'); });
  $('cohort-run').addEventListener('click', runCohort);
}
function builderEl(c) { return document.querySelector(`.cohort-build[data-cohort="${c}"]`); }
function setupCohortBuilder(c) {
  const root = c === 'b' ? $('cohort-b-manual') : builderEl('a');
  const fieldSel = root.querySelector('.cohort-field');
  fieldSel.innerHTML = COHORT_FIELDS.map(([k, l]) => `<option value="${k}">${l}</option>`).join('');
  fieldSel.value = cohort[c].field;
  fieldSel.addEventListener('change', () => { cohort[c].field = fieldSel.value; cohort[c].values = new Set(); renderCohortValues(c); });
  root.querySelector('.cohort-valsearch').addEventListener('input', () => renderCohortValues(c));
  renderCohortValues(c);
}
function fieldCounts(field) {
  const c = {}; (window.gemMetadata || []).forEach(m => { const v = clean(m[field]); if (v) c[v] = (c[v] || 0) + 1; });
  return Object.entries(c).sort((a, b) => b[1] - a[1]);
}
function renderCohortValues(c) {
  const root = c === 'b' ? $('cohort-b-manual') : builderEl('a');
  const box = root.querySelector('.cohort-values');
  const q = root.querySelector('.cohort-valsearch').value.trim().toLowerCase();
  let vals = fieldCounts(cohort[c].field);
  if (q) vals = vals.filter(([v]) => v.toLowerCase().includes(q));
  box.innerHTML = vals.slice(0, 300).map(([v, n]) =>
    `<label><input type="checkbox" value="${esc(v)}" ${cohort[c].values.has(v) ? 'checked' : ''}> ${esc(v)} <span class="cnt">${n}</span></label>`).join('') || '<div class="fba-hint-inline" style="padding:4px">no values</div>';
  box.querySelectorAll('input').forEach(cb => cb.addEventListener('change', () => { cb.checked ? cohort[c].values.add(cb.value) : cohort[c].values.delete(cb.value); updateCohortCount(c); }));
  updateCohortCount(c);
}
function cohortFiles(c) {
  const md = window.gemMetadata || [];
  if (c === 'b' && $('cohort-b-complement').checked) { const aset = new Set(cohortFiles('a')); return md.filter(m => !aset.has(m.gem_file)).map(m => m.gem_file); }
  const vals = cohort[c].values;
  if (!vals.size) return [];
  return md.filter(m => vals.has(clean(m[cohort[c].field]))).map(m => m.gem_file);
}
function updateCohortCount(c) { $('cohort-' + c + '-count').textContent = cohortFiles(c).length + ' GEMs'; }

async function runCohort() {
  const aFiles = cohortFiles('a'), bFiles = cohortFiles('b');
  if (aFiles.length < 3 || bFiles.length < 3) return setStatus('cohort-status', 'Each cohort needs ≥3 GEMs. Broaden your selection.', 'err');
  const cap = Math.max(3, Math.min(80, +$('cohort-cap').value || 25));
  const media = PRESETS[$('cohort-media').value].bounds;
  const meth = document.querySelector('input[name="cohort-method"]:checked').value;
  const sampA = sample(aFiles, cap), sampB = sample(bFiles, cap);
  const overlap = new Set(sampA).size + new Set(sampB).size - new Set([...sampA, ...sampB]).size;
  const jobs = [...sampA.map(f => ['A', f]), ...sampB.map(f => ['B', f])];
  $('cohort-run').disabled = true; $('cohort-results').style.display = 'none';
  setStatus('cohort-status', `Running ${meth.toUpperCase()} on ${jobs.length} models…`, 'busy'); prog('cohort-prog', 'cohort-prog-bar', 0);
  try {
    const t0 = performance.now();
    const rows = [];
    for (let i = 0; i < jobs.length; i++) {
      const [grp, f] = jobs[i];
      try { const model = await loadModel(f); const fba = await runFBA(model, media); let res = fba; if (meth === 'pfba' && fba.optimal && fba.growth > 1e-9) res = await runPFBA(model, media, fba); const ex = {}; for (const [id, v] of Object.entries(res.fluxes)) if (id.startsWith('EX_') && Math.abs(v) > 1e-6) ex[id] = v; rows.push({ grp, file: f, meta: meta(f), growth: fba.optimal ? fba.growth : 0, ex }); }
      catch (e) { rows.push({ grp, file: f, meta: meta(f), growth: 0, ex: {} }); }
      prog('cohort-prog', 'cohort-prog-bar', (i + 1) / jobs.length); setStatus('cohort-status', `Solved ${i + 1}/${jobs.length}…`, 'busy');
    }
    prog('cohort-prog', 'cohort-prog-bar', null);
    cohort.results = { rows, aFiles, bFiles, sampA, sampB, overlap, meth };
    renderCohort(cohort.results);
    setStatus('cohort-status', `Done — ${jobs.length} models in ${((performance.now() - t0) / 1000).toFixed(1)} s.` + (overlap ? ` (${overlap} strain(s) in both cohorts)` : ''), 'ok');
  } catch (e) { setStatus('cohort-status', 'Error: ' + e.message, 'err'); console.error(e); prog('cohort-prog', 'cohort-prog-bar', null); }
  finally { $('cohort-run').disabled = false; }
}
function sample(arr, n) { if (arr.length <= n) return arr.slice(); const a = arr.slice(), out = []; for (let k = 0; k < n; k++) { const j = k + Math.floor(seededRand(k) * (a.length - k)); [a[k], a[j]] = [a[j], a[k]]; out.push(a[k]); } return out; }
function labA() { return cohortLabel('a'); }
function labB() { return $('cohort-b-complement').checked ? 'B (complement)' : cohortLabel('b'); }
function cohortLabel(c) { const f = COHORT_FIELDS.find(x => x[0] === cohort[c].field)[1]; return `${f}=${[...cohort[c].values].join('/') || '—'}`; }

function renderCohort(R) {
  const A = R.rows.filter(r => r.grp === 'A'), B = R.rows.filter(r => r.grp === 'B');
  const gA = A.map(r => r.growth), gB = B.map(r => r.growth);
  const feasA = gA.filter(g => g > 1e-9), feasB = gB.filter(g => g > 1e-9);
  const meanA = feasA.length ? feasA.reduce((s, v) => s + v, 0) / feasA.length : 0;
  const meanB = feasB.length ? feasB.reduce((s, v) => s + v, 0) / feasB.length : 0;
  const gTest = mannWhitneyU(gA, gB);
  $('cohort-results').style.display = 'block';
  $('cohort-kpis').innerHTML =
    `<div class="fba-kpi"><div class="v" style="color:#2c6fbb">${A.length}</div><div class="l">Cohort A models</div></div>
     <div class="fba-kpi"><div class="v" style="color:#c0392b">${B.length}</div><div class="l">Cohort B models</div></div>
     <div class="fba-kpi"><div class="v">${fmt(meanA)} / ${fmt(meanB)}</div><div class="l">Mean growth A / B</div></div>
     <div class="fba-kpi"><div class="v" style="color:${gTest.p < 0.05 ? '#1a7f4b' : '#666'}">${gTest.p < 1e-4 ? gTest.p.toExponential(1) : gTest.p.toFixed(3)}</div><div class="l">Growth MWU p-value</div></div>`;

  // growth box
  window.Plotly.newPlot('cohort-plot-growth', [
    { type: 'box', boxpoints: 'all', jitter: 0.5, name: 'A: ' + labA(), y: gA, marker: { color: '#2c6fbb', size: 5 }, line: { color: '#2c6fbb' }, text: A.map(r => r.meta.genome_name || r.file) },
    { type: 'box', boxpoints: 'all', jitter: 0.5, name: 'B: ' + labB(), y: gB, marker: { color: '#c0392b', size: 5 }, line: { color: '#c0392b' }, text: B.map(r => r.meta.genome_name || r.file) },
  ], { margin: { l: 45, r: 10, t: 10, b: 30 }, height: 340, yaxis: { title: 'Growth (h⁻¹)', rangemode: 'tozero' }, showlegend: true, legend: { font: { size: 9 }, orientation: 'h', y: 1.12 }, font: { size: 11 } }, { responsive: true, displaylogo: false });

  // differential exchange fluxes
  const exIds = new Set(); R.rows.forEach(r => Object.keys(r.ex).forEach(id => exIds.add(id)));
  const diff = [];
  for (const id of exIds) {
    const va = A.map(r => r.ex[id] || 0), vb = B.map(r => r.ex[id] || 0);
    const nzA = va.filter(v => Math.abs(v) > 1e-9).length, nzB = vb.filter(v => Math.abs(v) > 1e-9).length;
    if (nzA + nzB < 2) continue;
    const mA = va.reduce((s, v) => s + v, 0) / va.length, mB = vb.reduce((s, v) => s + v, 0) / vb.length;
    const t = mannWhitneyU(va, vb);
    diff.push({ id, meanA: mA, meanB: mB, delta: mA - mB, p: t.p });
  }
  const qs = benjaminiHochberg(diff.map(d => d.p));
  diff.forEach((d, i) => d.q = qs[i]);
  diff.sort((a, b) => a.p - b.p);
  _cohortDiff = diff;
  renderVolcano(diff);
  renderCohortBars(diff);
  renderCohortPCA(R.rows);
  renderCohortTable(diff);
  $('cohort-detail').innerHTML = '';
  $('cohort-csv').onclick = () => { let c = 'exchange_id,metabolite,mean_flux_A,mean_flux_B,delta_A_minus_B,mwu_p,bh_q\n'; diff.forEach(d => { c += `${d.id},${bio(d.id)},${d.meanA},${d.meanB},${d.delta},${d.p},${d.q}\n`; }); saveCSV(c, `cohort_${$('cohort-media').value}`); };
}
function renderVolcano(diff) {
  const sig = diff.filter(d => d.q < 0.05), ns = diff.filter(d => !(d.q < 0.05));
  const mk = (arr, color) => ({ type: 'scatter', mode: 'markers', x: arr.map(d => d.delta), y: arr.map(d => -Math.log10(Math.max(d.p, 1e-300))), text: arr.map(d => bio(d.id)), customdata: arr.map(d => d.id), marker: { color, size: 8, line: { color: '#fff', width: 0.4 } }, hovertemplate: '%{text}<br>Δ(A−B) %{x:.3f}<br>-log10 p %{y:.2f}<extra></extra>' });
  window.Plotly.newPlot('cohort-plot-volcano', [
    Object.assign(mk(ns, '#b8c2cf'), { name: 'ns' }),
    Object.assign(mk(sig, '#7d3c98'), { name: 'FDR<0.05' }),
  ], { margin: { l: 45, r: 10, t: 10, b: 40 }, height: 340, xaxis: { title: 'Δ mean flux (A − B)', zeroline: true }, yaxis: { title: '−log₁₀ p' }, legend: { font: { size: 9 } }, font: { size: 11 } }, { responsive: true, displaylogo: false })
    .then(gd => gd.on('plotly_click', ev => showCohortDetail(ev.points[0].customdata)));
}
function renderCohortBars(diff) {
  const top = diff.slice().sort((a, b) => Math.abs(b.delta) - Math.abs(a.delta)).slice(0, 12).reverse();
  window.Plotly.newPlot('cohort-plot-bars', [
    { type: 'bar', orientation: 'h', name: 'A', y: top.map(d => bio(d.id)), x: top.map(d => d.meanA), marker: { color: '#2c6fbb' } },
    { type: 'bar', orientation: 'h', name: 'B', y: top.map(d => bio(d.id)), x: top.map(d => d.meanB), marker: { color: '#c0392b' } },
  ], { barmode: 'group', margin: { l: 70, r: 10, t: 10, b: 35 }, height: 340, xaxis: { title: 'Mean flux' }, legend: { font: { size: 9 } }, font: { size: 10 } }, { responsive: true, displaylogo: false });
}
function renderCohortPCA(rows) {
  const feats = commonExchanges(rows, 2);
  if (feats.length < 2 || rows.length < 4) { $('cohort-plot-pca').innerHTML = '<div class="fba-hint-inline" style="padding:1rem">Not enough data for PCA.</div>'; return; }
  const n = rows.length;
  let X = rows.map(r => feats.map(f => r.ex[f] || 0));
  const mean = feats.map((_, j) => X.reduce((s, row) => s + row[j], 0) / n);
  X = X.map(row => row.map((v, j) => v - mean[j]));
  const G = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, k) => { let s = 0; for (let j = 0; j < feats.length; j++) s += X[i][j] * X[k][j]; return s; }));
  const { values, vectors } = jacobiEigen(G);
  const order = values.map((v, i) => [v, i]).sort((a, b) => b[0] - a[0]);
  const tot = values.reduce((s, v) => s + Math.max(0, v), 0) || 1;
  const pc = (rank) => { const [lam, idx] = order[rank]; const sq = Math.sqrt(Math.max(0, lam)); return { s: vectors.map(row => row[idx] * sq), ev: Math.max(0, lam) / tot }; };
  const p1 = pc(0), p2 = pc(1);
  const grp = (g, color, name) => { const idxs = rows.map((r, i) => [r, i]).filter(([r]) => r.grp === g).map(([, i]) => i); return { type: 'scatter', mode: 'markers', name, x: idxs.map(i => p1.s[i]), y: idxs.map(i => p2.s[i]), text: idxs.map(i => rows[i].meta.genome_name || rows[i].file), customdata: idxs.map(i => rows[i].file), marker: { color, size: 9, line: { color: '#fff', width: 0.5 } }, hovertemplate: '%{text}<extra>' + name + '</extra>' }; };
  window.Plotly.newPlot('cohort-plot-pca', [grp('A', '#2c6fbb', 'Cohort A'), grp('B', '#c0392b', 'Cohort B')],
    { margin: { l: 45, r: 10, t: 10, b: 40 }, height: 340, xaxis: { title: `PC1 (${(p1.ev * 100).toFixed(0)}%)` }, yaxis: { title: `PC2 (${(p2.ev * 100).toFixed(0)}%)` }, legend: { font: { size: 9 } }, font: { size: 11 } }, { responsive: true, displaylogo: false });
}
function renderCohortTable(diff) {
  const rows = diff.slice(0, 40).map(d => `<tr><td><code>${esc(d.id)}</code></td><td>${esc(bio(d.id))}</td><td class="num">${fmt(d.meanA)}</td><td class="num">${fmt(d.meanB)}</td><td class="num" style="color:${d.delta >= 0 ? '#2c6fbb' : '#c0392b'}">${(d.delta >= 0 ? '+' : '') + fmt(d.delta)}</td><td class="num">${d.p < 1e-4 ? d.p.toExponential(1) : d.p.toFixed(4)}</td><td class="num" style="color:${d.q < 0.05 ? '#1a7f4b' : '#999'}">${d.q < 1e-4 ? d.q.toExponential(1) : d.q.toFixed(4)}</td></tr>`).join('');
  $('cohort-table').innerHTML = `<div class="fba-tablewrap" style="max-height:320px"><table class="fba-flux"><thead><tr><th>Exchange</th><th>Metabolite</th><th>Mean A</th><th>Mean B</th><th>Δ(A−B)</th><th>MWU p</th><th>BH q</th></tr></thead><tbody>${rows}</tbody></table></div>`;
}
function showCohortDetail(exId) {
  const d = cohort.results && cohortDiffLookup(exId); if (!d) return;
  $('cohort-detail').innerHTML = `<div class="fba-note" style="background:var(--accent)"><strong>${esc(bio(exId))}</strong> <code>${esc(exId)}</code> — mean flux A ${fmt(d.meanA)}, B ${fmt(d.meanB)}, Δ ${fmt(d.delta)} · MWU p ${d.p.toExponential(2)} · BH q ${d.q.toExponential(2)} ${d.q < 0.05 ? '<strong style="color:#1a7f4b">(significant)</strong>' : ''}</div>`;
}
let _cohortDiff = null;
function cohortDiffLookup(id) { return (_cohortDiff || []).find(d => d.id === id); }

// ── init ──────────────────────────────────────────────────────────────────────
async function init() {
  await presets();
  initTabs(); initFVA(); initDFBA(); initMulti(); initCohort();
  // URL param handoff: ?tab=&model=&slot=&models=
  const q = new URLSearchParams(location.search);
  if (q.get('tab')) { const btn = document.querySelector(`#studio-tabs .studio-tab[data-tab="${q.get('tab')}"]`); if (btn) btn.click(); }
  if (q.get('models')) { q.get('models').split(',').filter(Boolean).forEach(addModel); document.querySelector('#studio-tabs .studio-tab[data-tab="multi"]').click(); }
}
if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
