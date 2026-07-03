// EcopanGEM — Flux Balance Analysis UI controller (ES module).
// Wires the FBA panel: model + media pickers, Run button, results tables,
// Chart.js charts, and the Escher flux map. Depends on globals from the main
// page script: gemBatchMap, BATCH_URL_BASE, JSZip, Chart, escher.
import { runFBA, runPFBA, exchangeReport } from './fba_engine.js';

const $ = (id) => document.getElementById(id);
const BIOMASS_HINT = 'BIOMASS';
const state = { model: null, modelFile: null, media: null, presets: null, charts: {}, escher: null, lastResult: null };

// ── Load media presets ───────────────────────────────────────────────────────
async function loadPresets() {
  if (state.presets) return state.presets;
  state.presets = await (await fetch('fba/media_presets.json')).json();
  const sel = $('fba-media');
  sel.innerHTML = '';
  for (const [key, v] of Object.entries(state.presets)) {
    const o = document.createElement('option');
    o.value = key; o.textContent = v.label;
    sel.appendChild(o);
  }
  sel.addEventListener('change', updateMediaDesc);
  updateMediaDesc();
  return state.presets;
}
function updateMediaDesc() {
  const p = state.presets[$('fba-media').value];
  $('fba-media-desc').textContent = p ? `${p.desc}  ·  ${Object.keys(p.bounds).length} open exchanges  ·  ${p.source}` : '';
}

// ── Model picker (datalist populated from metadata) ──────────────────────────
function populateModels(meta) {
  const dl = $('fba-model-list');
  if (!dl || !meta) return;
  const frag = document.createDocumentFragment();
  meta.forEach((m) => {
    const o = document.createElement('option');
    o.value = m.gem_file;
    o.label = m.genome_name || '';
    frag.appendChild(o);
  });
  dl.innerHTML = '';
  dl.appendChild(frag);
  $('fba-model-hint').textContent = `${meta.length.toLocaleString()} models available — type a genome ID or pick a GEM from the table below.`;
}

// ── Load a model JSON from its batch zip (cached) ────────────────────────────
const modelCache = new Map();
async function loadModel(gemFile) {
  if (modelCache.has(gemFile)) return modelCache.get(gemFile);
  const batchNum = window.gemBatchMap && window.gemBatchMap[gemFile];
  if (!batchNum) throw new Error(`Unknown model "${gemFile}". Check the exact GEM file name (e.g. 562.12345.json.json).`);
  const url = window.BATCH_URL_BASE + String(batchNum).padStart(2, '0') + '.zip';
  const resp = await fetch(url);
  if (!resp.ok) throw new Error('Failed to download model batch.');
  const zip = await window.JSZip.loadAsync(await resp.blob());
  const entry = zip.file(gemFile);
  if (!entry) throw new Error('Model not found inside batch archive.');
  const model = JSON.parse(await entry.async('text'));
  modelCache.set(gemFile, model);
  return model;
}

// ── Run ──────────────────────────────────────────────────────────────────────
async function run() {
  const gemFile = $('fba-model').value.trim();
  const mediaKey = $('fba-media').value;
  const method = document.querySelector('input[name="fba-method"]:checked').value;
  if (!gemFile) { setStatus('Please choose a model first.', 'err'); return; }

  const btn = $('fba-run');
  btn.disabled = true;
  $('fba-results').style.display = 'none';
  setStatus('Loading model…', 'busy');
  try {
    const model = await loadModel(gemFile);
    state.model = model; state.modelFile = gemFile;
    const media = state.presets[mediaKey];

    setStatus(`Solving ${method.toUpperCase()} (in your browser)…`, 'busy');
    const t0 = performance.now();
    const fba = await runFBA(model, media.bounds);
    let result = fba;
    if (method === 'pfba' && fba.optimal && fba.growth > 1e-9) {
      result = await runPFBA(model, media.bounds, fba);
    }
    const solveMs = performance.now() - t0;

    // growth across all media (for the overview chart) — quick extra FBAs
    setStatus('Comparing growth across media…', 'busy');
    const across = {};
    for (const [k, v] of Object.entries(state.presets)) {
      across[k] = k === mediaKey ? fba.growth : (await runFBA(model, v.bounds)).growth;
    }

    state.lastResult = { gemFile, mediaKey, method, fba, result, across, solveMs };
    renderResults(state.lastResult);
    setStatus(`Done — ${method.toUpperCase()} solved in ${solveMs.toFixed(0)} ms on your machine.`, 'ok');
  } catch (e) {
    setStatus('Error: ' + e.message, 'err');
    console.error(e);
  } finally {
    btn.disabled = false;
  }
}

function setStatus(msg, kind) {
  const el = $('fba-status');
  el.textContent = msg;
  el.className = 'fba-status ' + (kind || '');
}

// ── Render results ────────────────────────────────────────────────────────────
function renderResults(R) {
  const { model, } = state;
  const fluxes = R.result.fluxes;
  const rep = exchangeReport(model, fluxes);
  const growth = R.result.growth || 0;
  const infeasible = !R.fba.optimal || growth <= 1e-9;

  $('fba-results').style.display = 'block';
  const label = state.presets[R.mediaKey].label;
  $('fba-result-title').innerHTML =
    `<code>${R.gemFile}</code> on <strong>${label}</strong> <span class="fba-badge">${R.method.toUpperCase()}</span>`;

  // KPI cards
  const totalFlux = R.result.pfba ? R.result.totalFlux : null;
  $('fba-kpis').innerHTML = [
    kpi(infeasible ? '0' : growth.toFixed(4), 'Growth rate (h⁻¹)', infeasible ? '#c0392b' : '#1a7f4b'),
    kpi(infeasible ? 'NO GROWTH' : 'FEASIBLE', 'Status', infeasible ? '#c0392b' : '#1a7f4b'),
    kpi(rep.uptake.length, 'Nutrients taken up'),
    kpi(rep.secretion.length, 'Products secreted'),
    kpi(totalFlux != null ? totalFlux.toFixed(0) : '—', 'Total flux Σ|v| (pFBA)'),
  ].join('');

  if (infeasible) {
    $('fba-charts').style.display = 'none';
    $('fba-map-wrap').style.display = 'none';
    $('fba-tables').innerHTML = `<div class="fba-note" style="color:#c0392b">This model cannot produce biomass on the selected medium (growth ≈ 0). Try a richer medium (e.g. Serum or Feces) or M9 + glucose.</div>`;
    return;
  }
  $('fba-charts').style.display = 'grid';
  $('fba-map-wrap').style.display = 'block';

  drawAcrossChart(R.across);
  drawExchangeChart('fba-chart-secretion', rep.secretion.slice(0, 10), '#c0392b', 'Secretion flux (mmol·gDW⁻¹·h⁻¹)');
  drawExchangeChart('fba-chart-uptake', rep.uptake.slice(0, 10).map(x => ({ ...x, flux: Math.abs(x.flux) })), '#2c6fbb', 'Uptake flux (mmol·gDW⁻¹·h⁻¹)');

  renderTables(rep);
  renderEscher(fluxes);
}

function kpi(val, lbl, color) {
  return `<div class="fba-kpi"><div class="v" ${color ? `style="color:${color}"` : ''}>${val}</div><div class="l">${lbl}</div></div>`;
}

// ── Charts (Chart.js) ─────────────────────────────────────────────────────────
function destroyChart(id) { if (state.charts[id]) { state.charts[id].destroy(); delete state.charts[id]; } }

function drawAcrossChart(across) {
  const labels = Object.keys(across).map(k => state.presets[k].label);
  const vals = Object.values(across).map(v => Math.max(0, v));
  const sel = state.lastResult.mediaKey;
  const colors = Object.keys(across).map(k => k === sel ? '#1a7f4b' : '#9db8d6');
  destroyChart('across');
  state.charts['across'] = new Chart($('fba-chart-across'), {
    type: 'bar',
    data: { labels, datasets: [{ data: vals, backgroundColor: colors, borderRadius: 4 }] },
    options: baseOpts('Predicted growth rate (h⁻¹)', false),
  });
}

function drawExchangeChart(canvasId, items, color, axisTitle) {
  const key = canvasId;
  const labels = items.map(x => bioLabel(x));
  const vals = items.map(x => x.flux);
  destroyChart(key);
  state.charts[key] = new Chart($(canvasId), {
    type: 'bar',
    data: { labels, datasets: [{ data: vals, backgroundColor: color, borderRadius: 4 }] },
    options: baseOpts(axisTitle, true),
  });
}

function bioLabel(x) {
  const met = x.id.replace(/^EX_/, '').replace(/_e$/, '');
  return `${met}`;
}

function baseOpts(axisTitle, horizontal) {
  return {
    indexAxis: horizontal ? 'y' : 'x',
    responsive: true, maintainAspectRatio: false,
    plugins: {
      legend: { display: false },
      tooltip: { callbacks: { label: (c) => ` ${c.parsed[horizontal ? 'x' : 'y'].toFixed(3)}` } },
    },
    scales: {
      [horizontal ? 'x' : 'y']: { title: { display: true, text: axisTitle, font: { size: 11 } }, ticks: { font: { size: 10 } } },
      [horizontal ? 'y' : 'x']: { ticks: { font: { size: 10 } } },
    },
  };
}

// ── Tables + CSV ──────────────────────────────────────────────────────────────
function renderTables(rep) {
  const rows = (arr, sign) => arr.map(x =>
    `<tr><td><code>${x.id}</code></td><td>${x.name || ''}</td><td class="num" style="color:${sign > 0 ? '#c0392b' : '#2c6fbb'}">${x.flux.toFixed(4)}</td></tr>`).join('');
  $('fba-tables').innerHTML = `
    <div class="fba-two">
      <div>
        <h6>Uptake (${rep.uptake.length})</h6>
        <div class="fba-tablewrap"><table class="fba-flux"><thead><tr><th>Exchange</th><th>Name</th><th>Flux</th></tr></thead>
          <tbody>${rows(rep.uptake, -1)}</tbody></table></div>
      </div>
      <div>
        <h6>Secretion / end products (${rep.secretion.length})</h6>
        <div class="fba-tablewrap"><table class="fba-flux"><thead><tr><th>Exchange</th><th>Name</th><th>Flux</th></tr></thead>
          <tbody>${rows(rep.secretion, 1)}</tbody></table></div>
      </div>
    </div>
    <button class="btn btn-sm btn-outline-secondary mt-2" id="fba-dl-csv">⬇ Download all fluxes (CSV)</button>`;
  $('fba-dl-csv').addEventListener('click', downloadFluxCSV);
}

function downloadFluxCSV() {
  const R = state.lastResult;
  const model = state.model;
  const nameById = {}; model.reactions.forEach(r => nameById[r.id] = r.name || '');
  let csv = 'reaction_id,name,flux\n';
  for (const [id, v] of Object.entries(R.result.fluxes)) {
    csv += `${id},"${(nameById[id] || '').replace(/"/g, '""')}",${v}\n`;
  }
  const blob = new Blob([csv], { type: 'text/csv' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = `FBA_${R.gemFile.replace(/\.json.*/, '')}_${R.mediaKey}_${R.method}.csv`;
  a.click();
}

// ── Escher flux map ───────────────────────────────────────────────────────────
let coreMapData = null;
async function renderEscher(fluxes) {
  const container = $('fba-map');
  if (!window.escher) { $('fba-map-note').textContent = 'Escher library unavailable.'; return; }
  try {
    if (!coreMapData) coreMapData = await (await fetch('fba/core_map.json')).json();
    container.innerHTML = '';
    const sel = window.escher.libs.d3_select('#fba-map');
    state.escher = window.escher.Builder(coreMapData, null, null, sel, {
      menu: 'zoom', scroll_behavior: 'zoom', fill_screen: false,
      never_ask_before_quit: true, reaction_data: fluxes,
      reaction_styles: ['color', 'size', 'text'],
      reaction_scale: [
        { type: 'min', color: '#d0d0d0', size: 6 },
        { type: 'value', value: 0, color: '#d0d0d0', size: 6 },
        { type: 'median', color: '#5b8ff9', size: 14 },
        { type: 'max', color: '#c0392b', size: 28 },
      ],
      tooltip_component: undefined,
      hide_secondary_metabolites: false,
    });
    $('fba-map-note').textContent = 'Central-carbon metabolism (e_coli_core map). Line colour & thickness = |flux|; hover a reaction for its value. Reactions absent from this core map are not shown.';
  } catch (e) {
    $('fba-map-note').textContent = 'Map render error: ' + e.message;
    console.error(e);
  }
}

// ── Init ──────────────────────────────────────────────────────────────────────
async function init() {
  await loadPresets();
  $('fba-run').addEventListener('click', run);
  if (window.gemMetadata) populateModels(window.gemMetadata);
  else window.addEventListener('gem-data-ready', () => populateModels(window.gemMetadata));
  // allow other code (detail modal) to preselect a model and run
  window.fbaSelectModel = (gemFile) => {
    $('fba-model').value = gemFile;
    document.getElementById('fba-section').scrollIntoView({ behavior: 'smooth' });
  };
}
if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
else init();
