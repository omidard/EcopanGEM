// EcopanGEM — Flux Balance Analysis UI controller (ES module).
// Single & comparative FBA/pFBA with editable media, knockouts, a rich model
// picker, charts (Chart.js) and Escher flux maps. All client-side.
// Depends on page globals: gemBatchMap, BATCH_URL_BASE, JSZip, Chart, escher.
import { runFBA, runPFBA, exchangeReport, listExchanges, listReactions } from './fba_engine.js';

const $ = (id) => document.getElementById(id);
const el = (tag, cls, html) => { const e = document.createElement(tag); if (cls) e.className = cls; if (html != null) e.innerHTML = html; return e; };
const fmt = (x, n = 4) => (x == null || isNaN(x)) ? '—' : Number(x).toFixed(n);

const S = {
  presets: null,
  mode: 'single',
  charts: {},
  conditions: {},   // slot -> condition state
  lastRun: null,
};

// ── Media presets ─────────────────────────────────────────────────────────────
async function loadPresets() {
  if (!S.presets) S.presets = await (await fetch('fba/media_presets.json')).json();
  return S.presets;
}

// ── Model loading (batch zip, cached) ────────────────────────────────────────
const modelCache = new Map();
async function loadModel(gemFile) {
  if (modelCache.has(gemFile)) return modelCache.get(gemFile);
  const batchNum = window.gemBatchMap && window.gemBatchMap[gemFile];
  if (!batchNum) throw new Error(`Unknown model "${gemFile}".`);
  const url = window.BATCH_URL_BASE + String(batchNum).padStart(2, '0') + '.zip';
  const resp = await fetch(url);
  if (!resp.ok) throw new Error('Failed to download model batch.');
  const zip = await window.JSZip.loadAsync(await resp.blob());
  const entry = zip.file(gemFile);
  if (!entry) throw new Error('Model not found inside batch.');
  const model = JSON.parse(await entry.async('text'));
  modelCache.set(gemFile, model);
  return model;
}
const metaByFile = new Map();
function indexMeta() {
  if (metaByFile.size || !window.gemMetadata) return;
  window.gemMetadata.forEach(m => metaByFile.set(m.gem_file, m));
}

// ── Condition card ────────────────────────────────────────────────────────────
function newCondition(slot) {
  return { slot, modelFile: null, model: null, meta: null, exchanges: null, reactions: null,
           mediaKey: 'M9_glucose_aerobic', mediaBounds: {}, mediaEdited: false, knockouts: new Set(), card: null };
}

function buildConditionCard(slot) {
  const cond = newCondition(slot);
  S.conditions[slot] = cond;
  const card = el('div', 'fba-cond'); card.dataset.slot = slot;
  card.innerHTML = `
    <div class="fba-cond-head"><span class="fba-cond-badge">${slot.toUpperCase()}</span>
      <span class="title">Condition ${slot.toUpperCase()}</span></div>
    <div class="fba-field">
      <label>Strain model (GEM)</label>
      <div class="fba-combo">
        <input class="form-control form-control-sm fba-combo-input" placeholder="Search genome name, strain or ID…" autocomplete="off">
        <div class="fba-combo-menu"></div>
      </div>
      <div class="fba-modelcard" style="display:none"></div>
    </div>
    <div class="fba-field">
      <label>Growth medium
        <button type="button" class="fba-linkbtn me-toggle" disabled>edit compounds</button>
        <button type="button" class="fba-linkbtn me-reset" disabled>reset</button></label>
      <select class="form-select form-select-sm fba-media-sel"></select>
      <div class="fba-media-desc"></div>
      <div class="fba-media-editor" style="display:none"></div>
    </div>
    <div class="fba-field">
      <label>Reaction knockouts <span class="fba-hint-inline">(optional — set flux to 0)</span></label>
      <div class="fba-ko">
        <input class="form-control form-control-sm fba-ko-input" placeholder="Search reaction ID or name to knock out…" autocomplete="off" disabled>
        <div class="fba-ko-menu"></div>
        <div class="fba-ko-chips"></div>
      </div>
    </div>`;
  cond.card = card;

  // media select
  const msel = card.querySelector('.fba-media-sel');
  for (const [k, v] of Object.entries(S.presets)) {
    const o = el('option'); o.value = k; o.textContent = v.label; msel.appendChild(o);
  }
  msel.value = cond.mediaKey;
  updateMediaDesc(cond);
  msel.addEventListener('change', () => {
    cond.mediaKey = msel.value;
    cond.mediaBounds = { ...S.presets[cond.mediaKey].bounds };
    cond.mediaEdited = false;
    updateMediaDesc(cond);
    if (card.querySelector('.fba-media-editor').style.display !== 'none') renderMediaEditor(cond);
  });

  wireCombo(cond);
  wireMediaEditor(cond);
  wireKO(cond);
  return card;
}

function updateMediaDesc(cond) {
  const p = S.presets[cond.mediaKey];
  const n = Object.keys(cond.mediaBounds).length || Object.keys(p.bounds).length;
  cond.card.querySelector('.fba-media-desc').textContent =
    `${p.desc} · ${n} open exchanges${cond.mediaEdited ? ' (edited)' : ''} · ${p.source}`;
}

// ── Model combobox ────────────────────────────────────────────────────────────
function wireCombo(cond) {
  const input = cond.card.querySelector('.fba-combo-input');
  const menu = cond.card.querySelector('.fba-combo-menu');
  let items = [], active = -1;

  const render = () => {
    const q = input.value.trim().toLowerCase();
    const meta = window.gemMetadata || [];
    items = (q ? meta.filter(m =>
        (m.gem_file || '').toLowerCase().includes(q) ||
        (m.genome_name || '').toLowerCase().includes(q) ||
        (m.strain || '').toLowerCase().includes(q))
      : meta).slice(0, 40);
    if (!items.length) { menu.innerHTML = `<div class="fba-combo-empty">No models match.</div>`; menu.classList.add('show'); return; }
    menu.innerHTML = items.map((m, i) => `
      <div class="fba-combo-item${i === active ? ' active' : ''}" data-i="${i}">
        <div class="nm">${esc(m.genome_name || m.strain || m.gem_file)}</div>
        <div class="id">${esc(m.gem_file)}</div>
        <div class="meta">${m.phylogroup ? 'Phylogroup ' + esc(m.phylogroup) : ''}${m.MLST ? ' · ST' + esc(m.MLST) : ''}${m.isolation_source ? ' · ' + esc(m.isolation_source) : ''}</div>
      </div>`).join('');
    menu.classList.add('show');
  };
  input.addEventListener('focus', render);
  input.addEventListener('input', () => { active = -1; render(); });
  input.addEventListener('keydown', (e) => {
    if (!menu.classList.contains('show')) return;
    if (e.key === 'ArrowDown') { active = Math.min(active + 1, items.length - 1); render(); e.preventDefault(); }
    else if (e.key === 'ArrowUp') { active = Math.max(active - 1, 0); render(); e.preventDefault(); }
    else if (e.key === 'Enter') { if (items[active] || items[0]) { pick((items[active] || items[0]).gem_file); e.preventDefault(); } }
    else if (e.key === 'Escape') menu.classList.remove('show');
  });
  menu.addEventListener('mousedown', (e) => {
    const it = e.target.closest('.fba-combo-item'); if (!it) return;
    e.preventDefault(); pick(items[+it.dataset.i].gem_file);
  });
  input.addEventListener('blur', () => setTimeout(() => menu.classList.remove('show'), 150));

  const pick = async (gemFile) => {
    menu.classList.remove('show');
    input.value = '';
    await selectModel(cond, gemFile);
  };
}

async function selectModel(cond, gemFile) {
  indexMeta();
  cond.modelFile = gemFile;
  cond.meta = metaByFile.get(gemFile) || {};
  const mc = cond.card.querySelector('.fba-modelcard');
  mc.style.display = 'block';
  mc.innerHTML = `<span class="fba-hint-inline">Loading model…</span>`;
  try {
    const model = await loadModel(gemFile);
    cond.model = model;
    cond.exchanges = listExchanges(model);
    cond.reactions = listReactions(model);
    cond.mediaBounds = { ...S.presets[cond.mediaKey].bounds };
    cond.mediaEdited = false;
    const m = cond.meta;
    mc.innerHTML = `
      <div class="mc-name">${esc(m.genome_name || m.strain || gemFile)}</div>
      <div class="mc-id">${esc(gemFile)}</div>
      <div class="mc-tags">
        <span class="fba-tag">${model.reactions.length} rxns</span>
        <span class="fba-tag">${model.metabolites.length} mets</span>
        <span class="fba-tag">${model.genes.length} genes</span>
        ${m.phylogroup ? `<span class="fba-tag">Phylogroup ${esc(m.phylogroup)}</span>` : ''}
        ${m.MLST ? `<span class="fba-tag">ST ${esc(m.MLST)}</span>` : ''}
        ${m.isolation_source ? `<span class="fba-tag">${esc(m.isolation_source)}</span>` : ''}
        ${m.isolation_country ? `<span class="fba-tag">${esc(m.isolation_country)}</span>` : ''}
      </div>`;
    // enable editors
    cond.card.querySelector('.me-toggle').disabled = false;
    cond.card.querySelector('.me-reset').disabled = false;
    cond.card.querySelector('.fba-ko-input').disabled = false;
    updateMediaDesc(cond);
  } catch (e) {
    mc.innerHTML = `<span style="color:#c0392b">${esc(e.message)}</span>`;
  }
}

// ── Media editor ──────────────────────────────────────────────────────────────
function wireMediaEditor(cond) {
  const toggle = cond.card.querySelector('.me-toggle');
  const reset = cond.card.querySelector('.me-reset');
  const box = cond.card.querySelector('.fba-media-editor');
  toggle.addEventListener('click', () => {
    if (box.style.display === 'none') { box.style.display = 'block'; renderMediaEditor(cond); toggle.textContent = 'hide compounds'; }
    else { box.style.display = 'none'; toggle.textContent = 'edit compounds'; }
  });
  reset.addEventListener('click', () => {
    cond.mediaBounds = { ...S.presets[cond.mediaKey].bounds };
    cond.mediaEdited = false; updateMediaDesc(cond);
    if (box.style.display !== 'none') renderMediaEditor(cond);
  });
}

function renderMediaEditor(cond) {
  const box = cond.card.querySelector('.fba-media-editor');
  const nameById = {}; (cond.exchanges || []).forEach(e => nameById[e.id] = e.name);
  const rows = Object.entries(cond.mediaBounds)
    .sort((a, b) => a[0].localeCompare(b[0]))
    .map(([id, lb]) => `
      <tr data-id="${id}">
        <td class="me-cmpd">${esc((nameById[id] || '').replace(/ exchange$/i, ''))}</td>
        <td class="me-id">${esc(id)}</td>
        <td><input type="number" step="0.1" class="me-rate" value="${lb}"></td>
        <td><span class="fba-me-rm" title="remove">✕</span></td>
      </tr>`).join('');
  box.innerHTML = `
    <div class="me-head"><b>Uptake bounds (mmol·gDW⁻¹·h⁻¹) — negative = uptake</b>
      <span class="fba-hint-inline">${Object.keys(cond.mediaBounds).length} compounds</span></div>
    <div class="fba-me-list"><table class="fba-me"><tbody>${rows}</tbody></table></div>
    <div class="fba-me-add">
      <input class="form-control form-control-sm me-add-input" placeholder="+ add compound (search exchange)…" autocomplete="off">
      <div class="fba-combo-menu me-add-menu"></div>
    </div>`;
  // rate edits
  box.querySelectorAll('.me-rate').forEach(inp => inp.addEventListener('change', () => {
    const id = inp.closest('tr').dataset.id;
    const v = parseFloat(inp.value);
    if (!isNaN(v)) { cond.mediaBounds[id] = v; cond.mediaEdited = true; updateMediaDesc(cond); }
  }));
  // remove
  box.querySelectorAll('.fba-me-rm').forEach(x => x.addEventListener('click', () => {
    const id = x.closest('tr').dataset.id;
    delete cond.mediaBounds[id]; cond.mediaEdited = true; renderMediaEditor(cond); updateMediaDesc(cond);
  }));
  // add-compound combobox
  const ai = box.querySelector('.me-add-input'), am = box.querySelector('.me-add-menu');
  let addItems = [];
  const renderAdd = () => {
    const q = ai.value.trim().toLowerCase();
    const avail = (cond.exchanges || []).filter(e => !(e.id in cond.mediaBounds));
    addItems = (q ? avail.filter(e => e.id.toLowerCase().includes(q) || (e.name || '').toLowerCase().includes(q)) : avail).slice(0, 30);
    am.innerHTML = addItems.length ? addItems.map((e, i) =>
      `<div class="fba-combo-item" data-i="${i}"><div class="nm">${esc((e.name || '').replace(/ exchange$/i, ''))}</div><div class="id">${esc(e.id)}</div></div>`).join('')
      : `<div class="fba-combo-empty">No more exchanges.</div>`;
    am.classList.add('show');
  };
  ai.addEventListener('focus', renderAdd); ai.addEventListener('input', renderAdd);
  ai.addEventListener('blur', () => setTimeout(() => am.classList.remove('show'), 150));
  am.addEventListener('mousedown', (e) => {
    const it = e.target.closest('.fba-combo-item'); if (!it) return; e.preventDefault();
    const ex = addItems[+it.dataset.i];
    cond.mediaBounds[ex.id] = -0.5; cond.mediaEdited = true; renderMediaEditor(cond); updateMediaDesc(cond);
  });
}

// ── Knockout editor ───────────────────────────────────────────────────────────
function wireKO(cond) {
  const input = cond.card.querySelector('.fba-ko-input');
  const menu = cond.card.querySelector('.fba-ko-menu');
  let items = [];
  const render = () => {
    const q = input.value.trim().toLowerCase();
    if (!q || !cond.reactions) { menu.classList.remove('show'); return; }
    items = cond.reactions.filter(r => !cond.knockouts.has(r.id) &&
      (r.id.toLowerCase().includes(q) || (r.name || '').toLowerCase().includes(q))).slice(0, 25);
    menu.innerHTML = items.length ? items.map((r, i) =>
      `<div class="fba-ko-item" data-i="${i}"><span class="id">${esc(r.id)}</span> <span class="nm">${esc(r.name || '')}</span></div>`).join('')
      : `<div class="fba-combo-empty">No match.</div>`;
    menu.classList.add('show');
  };
  input.addEventListener('input', render); input.addEventListener('focus', render);
  input.addEventListener('blur', () => setTimeout(() => menu.classList.remove('show'), 150));
  menu.addEventListener('mousedown', (e) => {
    const it = e.target.closest('.fba-ko-item'); if (!it) return; e.preventDefault();
    cond.knockouts.add(items[+it.dataset.i].id); input.value = ''; menu.classList.remove('show'); renderChips(cond);
  });
}
function renderChips(cond) {
  const wrap = cond.card.querySelector('.fba-ko-chips');
  wrap.innerHTML = [...cond.knockouts].map(id =>
    `<span class="fba-ko-chip">${esc(id)} <span class="x" data-id="${esc(id)}">✕</span></span>`).join('');
  wrap.querySelectorAll('.x').forEach(x => x.addEventListener('click', () => { cond.knockouts.delete(x.dataset.id); renderChips(cond); }));
}

// ── Run ───────────────────────────────────────────────────────────────────────
function method() { return document.querySelector('input[name="fba-method"]:checked').value; }

async function runAnalysis(cond, meth) {
  const opts = { knockouts: [...cond.knockouts] };
  const fba = await runFBA(cond.model, cond.mediaBounds, opts);
  let result = fba;
  if (meth === 'pfba' && fba.optimal && fba.growth > 1e-9) result = await runPFBA(cond.model, cond.mediaBounds, fba, opts);
  return { cond, fba, result, meth };
}

async function run() {
  const meth = method();
  const btn = $('fba-run'); btn.disabled = true;
  $('fba-results').style.display = 'none';
  try {
    if (S.mode === 'single') {
      const a = S.conditions.a;
      if (!a.model) return setStatus('Choose a model for Condition A first.', 'err');
      setStatus(`Solving ${meth.toUpperCase()} in your browser…`, 'busy');
      const t0 = performance.now();
      const RA = await runAnalysis(a, meth);
      // growth across media (with current knockouts)
      const across = {};
      for (const [k, v] of Object.entries(S.presets)) {
        across[k] = (k === a.mediaKey && !a.mediaEdited) ? RA.fba.growth
          : (await runFBA(a.model, v.bounds, { knockouts: [...a.knockouts] })).growth;
      }
      RA.across = across;
      S.lastRun = { mode: 'single', RA };
      renderSingle(RA);
      setStatus(`Done — solved in ${(performance.now() - t0).toFixed(0)} ms on your machine.`, 'ok');
    } else {
      const a = S.conditions.a, b = S.conditions.b;
      if (!a.model || !b.model) return setStatus('Choose a model for both Condition A and B.', 'err');
      setStatus(`Solving both conditions (${meth.toUpperCase()})…`, 'busy');
      const t0 = performance.now();
      const RA = await runAnalysis(a, meth);
      const RB = await runAnalysis(b, meth);
      S.lastRun = { mode: 'compare', RA, RB };
      renderCompare(RA, RB);
      setStatus(`Done — both solved in ${(performance.now() - t0).toFixed(0)} ms.`, 'ok');
    }
  } catch (e) { setStatus('Error: ' + e.message, 'err'); console.error(e); }
  finally { btn.disabled = false; }
}

function setStatus(msg, kind) { const e = $('fba-status'); e.textContent = msg; e.className = 'fba-status ' + (kind || ''); }
function kpi(v, l, color) { return `<div class="fba-kpi"><div class="v"${color ? ` style="color:${color}"` : ''}>${v}</div><div class="l">${l}</div></div>`; }
function condLabel(cond) { return `${S.presets[cond.mediaKey].label}${cond.mediaEdited ? ' *' : ''}${cond.knockouts.size ? ` · KO:${[...cond.knockouts].join(',')}` : ''}`; }

// ── Single results ────────────────────────────────────────────────────────────
function renderSingle(R) {
  const cond = R.cond, fluxes = R.result.fluxes, rep = exchangeReport(cond.model, fluxes);
  const growth = R.result.growth || 0, infeasible = !R.fba.optimal || growth <= 1e-9;
  const box = $('fba-results'); box.style.display = 'block';
  box.innerHTML = `
    <hr>
    <div style="font-size:0.95rem;margin-bottom:0.4rem"><code>${esc(cond.modelFile)}</code> on <strong>${esc(S.presets[cond.mediaKey].label)}</strong>
      <span class="fba-badge">${R.meth.toUpperCase()}</span> ${cond.knockouts.size ? `<span class="fba-badge" style="background:#c0392b">${cond.knockouts.size} KO</span>` : ''}</div>
    <div class="fba-kpis">
      ${kpi(infeasible ? '0' : fmt(growth), 'Growth rate (h⁻¹)', infeasible ? '#c0392b' : '#1a7f4b')}
      ${kpi(infeasible ? 'NO GROWTH' : 'FEASIBLE', 'Status', infeasible ? '#c0392b' : '#1a7f4b')}
      ${kpi(rep.uptake.length, 'Nutrients taken up')}
      ${kpi(rep.secretion.length, 'Products secreted')}
      ${kpi(R.result.pfba ? fmt(R.result.totalFlux, 0) : '—', 'Total flux Σ|v| (pFBA)')}
    </div>`;
  if (infeasible) {
    box.insertAdjacentHTML('beforeend', `<div class="fba-note" style="color:#c0392b">No biomass on this medium/knockout set (growth ≈ 0). Try a richer medium or remove knockouts.</div>`);
    return;
  }
  box.insertAdjacentHTML('beforeend', `
    <div class="fba-charts">
      <div class="fba-chart-card"><h6>Growth across media${cond.knockouts.size ? ' (with knockouts)' : ''}</h6><div class="fba-chart-box"><canvas id="ch-across"></canvas></div></div>
      <div class="fba-chart-card"><h6>Top secreted end products</h6><div class="fba-chart-box"><canvas id="ch-sec"></canvas></div></div>
      <div class="fba-chart-card"><h6>Top nutrient uptakes</h6><div class="fba-chart-box"><canvas id="ch-up"></canvas></div></div>
    </div>
    <div class="fba-map-wrap"><h6 style="font-size:0.85rem;font-weight:700;color:#444;margin-bottom:0.4rem">Metabolic flux map</h6>
      <div id="ch-map" style="width:100%;height:560px;border:1px solid #e5e8ec;border-radius:8px;background:#fff;overflow:hidden"></div>
      <div class="fba-map-note" id="ch-map-note"></div></div>
    <div id="ch-tables" style="margin-top:1.2rem"></div>`);

  barChart('ch-across', Object.keys(R.across).map(k => S.presets[k].label),
    Object.values(R.across).map(v => Math.max(0, v)),
    Object.keys(R.across).map(k => k === cond.mediaKey ? '#1a7f4b' : '#9db8d6'), 'Growth (h⁻¹)', false);
  barChart('ch-sec', rep.secretion.slice(0, 10).map(bio), rep.secretion.slice(0, 10).map(x => x.flux), '#c0392b', 'Flux', true);
  barChart('ch-up', rep.uptake.slice(0, 10).map(bio), rep.uptake.slice(0, 10).map(x => Math.abs(x.flux)), '#2c6fbb', 'Flux', true);
  renderEscher('ch-map', 'ch-map-note', fluxes);
  singleTables(cond, R, rep);
}

function singleTables(cond, R, rep) {
  const rows = (arr, s) => arr.map(x => `<tr><td><code>${esc(x.id)}</code></td><td>${esc(x.name || '')}</td><td class="num" style="color:${s > 0 ? '#c0392b' : '#2c6fbb'}">${fmt(x.flux)}</td></tr>`).join('');
  $('ch-tables').innerHTML = `
    <div class="fba-two">
      <div><h6>Uptake (${rep.uptake.length})</h6><div class="fba-tablewrap"><table class="fba-flux"><thead><tr><th>Exchange</th><th>Name</th><th>Flux</th></tr></thead><tbody>${rows(rep.uptake, -1)}</tbody></table></div></div>
      <div><h6>Secretion / end products (${rep.secretion.length})</h6><div class="fba-tablewrap"><table class="fba-flux"><thead><tr><th>Exchange</th><th>Name</th><th>Flux</th></tr></thead><tbody>${rows(rep.secretion, 1)}</tbody></table></div></div>
    </div>
    <button class="btn btn-sm btn-outline-secondary mt-2" id="ch-csv">⬇ Download all fluxes (CSV)</button>`;
  $('ch-csv').addEventListener('click', () => downloadFluxCSV(cond.model, R.result.fluxes, `FBA_${cond.modelFile.replace(/\.json.*/, '')}_${cond.mediaKey}_${R.meth}`));
}

// ── Compare results ───────────────────────────────────────────────────────────
function renderCompare(RA, RB) {
  const gA = RA.result.growth || 0, gB = RB.result.growth || 0;
  const fA = RA.result.fluxes, fB = RB.result.fluxes;
  const infA = !RA.fba.optimal || gA <= 1e-9, infB = !RB.fba.optimal || gB <= 1e-9;
  const delta = gB - gA, pct = gA > 1e-9 ? (delta / gA * 100) : null;
  const box = $('fba-results'); box.style.display = 'block';
  box.innerHTML = `
    <hr>
    <div class="fba-cmp-heads">
      <div><span class="fba-cond-badge">A</span> <code>${esc(RA.cond.modelFile)}</code> — ${esc(condLabel(RA.cond))}</div>
      <div><span class="fba-cond-badge" style="background:#c0392b">B</span> <code>${esc(RB.cond.modelFile)}</code> — ${esc(condLabel(RB.cond))}</div>
    </div>
    <div class="fba-kpis">
      ${kpi(infA ? '0' : fmt(gA), 'Growth A (h⁻¹)', '#2c6fbb')}
      ${kpi(infB ? '0' : fmt(gB), 'Growth B (h⁻¹)', '#c0392b')}
      ${kpi((delta >= 0 ? '+' : '') + fmt(delta), 'Δ Growth (B−A)', delta >= 0 ? '#1a7f4b' : '#c0392b')}
      ${kpi(pct == null ? '—' : (pct >= 0 ? '+' : '') + pct.toFixed(1) + '%', 'Δ Growth %')}
    </div>
    <div class="fba-charts" style="grid-template-columns:1fr 1fr">
      <div class="fba-chart-card"><h6>Growth: A vs B</h6><div class="fba-chart-box"><canvas id="cc-growth"></canvas></div></div>
      <div class="fba-chart-card"><h6>Secreted end products: A vs B</h6><div class="fba-chart-box"><canvas id="cc-sec"></canvas></div></div>
    </div>
    <div class="fba-two" style="margin-top:0.4rem">
      <div class="fba-map-wrap"><h6 style="font-size:0.82rem;font-weight:700;color:#2c6fbb">Flux map — Condition A</h6>
        <div id="cc-map-a" style="width:100%;height:460px;border:1px solid #e5e8ec;border-radius:8px;background:#fff;overflow:hidden"></div></div>
      <div class="fba-map-wrap"><h6 style="font-size:0.82rem;font-weight:700;color:#c0392b">Flux map — Condition B</h6>
        <div id="cc-map-b" style="width:100%;height:460px;border:1px solid #e5e8ec;border-radius:8px;background:#fff;overflow:hidden"></div></div>
    </div>
    <div class="fba-map-note">Each map is coloured by its own condition's |flux|. Reactions outside the e_coli_core central-metabolism map are not shown.</div>
    <div id="cc-diff" style="margin-top:1.2rem"></div>`;

  barChart('cc-growth', ['Condition A', 'Condition B'], [Math.max(0, gA), Math.max(0, gB)], ['#2c6fbb', '#c0392b'], 'Growth (h⁻¹)', false);

  // secretion comparison over union of top secretions
  const repA = exchangeReport(RA.cond.model, fA), repB = exchangeReport(RB.cond.model, fB);
  const topIds = [...new Set([...repA.secretion.slice(0, 8), ...repB.secretion.slice(0, 8)].map(x => x.id))];
  const mapA = Object.fromEntries(repA.secretion.map(x => [x.id, x.flux]));
  const mapB = Object.fromEntries(repB.secretion.map(x => [x.id, x.flux]));
  groupedChart('cc-sec', topIds.map(id => id.replace(/^EX_/, '').replace(/_e$/, '')),
    topIds.map(id => mapA[id] || 0), topIds.map(id => mapB[id] || 0), 'Secretion flux');

  renderEscher('cc-map-a', null, fA);
  renderEscher('cc-map-b', null, fB);
  diffTable(RA, RB, fA, fB);
}

function diffTable(RA, RB, fA, fB) {
  const ids = new Set([...Object.keys(fA), ...Object.keys(fB)]);
  const nameById = {}; RA.cond.model.reactions.forEach(r => nameById[r.id] = r.name || ''); RB.cond.model.reactions.forEach(r => { if (!nameById[r.id]) nameById[r.id] = r.name || ''; });
  const rows = [...ids].map(id => ({ id, a: fA[id] || 0, b: fB[id] || 0, d: (fB[id] || 0) - (fA[id] || 0), name: nameById[id] || '' }))
    .filter(r => Math.abs(r.d) > 1e-6).sort((x, y) => Math.abs(y.d) - Math.abs(x.d));
  const top = rows.slice(0, 60);
  $('cc-diff').innerHTML = `
    <h6>Largest flux differences (B − A) <span class="fba-hint-inline">${rows.length} reactions differ; top ${top.length} shown</span></h6>
    <div class="fba-tablewrap" style="max-height:360px"><table class="fba-flux"><thead><tr><th>Reaction</th><th>Name</th><th>Flux A</th><th>Flux B</th><th>Δ (B−A)</th></tr></thead>
      <tbody>${top.map(r => `<tr><td><code>${esc(r.id)}</code></td><td>${esc(r.name)}</td><td class="num">${fmt(r.a, 3)}</td><td class="num">${fmt(r.b, 3)}</td><td class="num" style="color:${r.d >= 0 ? '#1a7f4b' : '#c0392b'}">${(r.d >= 0 ? '+' : '') + fmt(r.d, 3)}</td></tr>`).join('')}</tbody></table></div>
    <button class="btn btn-sm btn-outline-secondary mt-2" id="cc-csv">⬇ Download flux comparison (CSV)</button>`;
  $('cc-csv').addEventListener('click', () => {
    let csv = 'reaction_id,name,flux_A,flux_B,delta_B_minus_A\n';
    rows.forEach(r => { csv += `${r.id},"${(r.name || '').replace(/"/g, '""')}",${r.a},${r.b},${r.d}\n`; });
    saveCSV(csv, `FBA_compare_${RA.cond.modelFile.replace(/\.json.*/, '')}_vs_${RB.cond.modelFile.replace(/\.json.*/, '')}`);
  });
}

// ── Charts ────────────────────────────────────────────────────────────────────
function destroyChart(id) { if (S.charts[id]) { S.charts[id].destroy(); delete S.charts[id]; } }
function barChart(canvasId, labels, data, colors, axisTitle, horizontal) {
  destroyChart(canvasId);
  S.charts[canvasId] = new Chart($(canvasId), {
    type: 'bar', data: { labels, datasets: [{ data, backgroundColor: colors, borderRadius: 4 }] },
    options: {
      indexAxis: horizontal ? 'y' : 'x', responsive: true, maintainAspectRatio: false,
      plugins: { legend: { display: false }, tooltip: { callbacks: { label: c => ' ' + c.parsed[horizontal ? 'x' : 'y'].toFixed(3) } } },
      scales: { [horizontal ? 'x' : 'y']: { title: { display: true, text: axisTitle, font: { size: 11 } }, ticks: { font: { size: 10 } } }, [horizontal ? 'y' : 'x']: { ticks: { font: { size: 10 } } } },
    },
  });
}
function groupedChart(canvasId, labels, dataA, dataB, axisTitle) {
  destroyChart(canvasId);
  S.charts[canvasId] = new Chart($(canvasId), {
    type: 'bar',
    data: { labels, datasets: [{ label: 'A', data: dataA, backgroundColor: '#2c6fbb', borderRadius: 3 }, { label: 'B', data: dataB, backgroundColor: '#c0392b', borderRadius: 3 }] },
    options: { indexAxis: 'y', responsive: true, maintainAspectRatio: false,
      plugins: { legend: { display: true, labels: { font: { size: 10 }, boxWidth: 12 } } },
      scales: { x: { title: { display: true, text: axisTitle, font: { size: 11 } }, ticks: { font: { size: 10 } } }, y: { ticks: { font: { size: 10 } } } } },
  });
}
function bio(x) { return x.id.replace(/^EX_/, '').replace(/_e$/, ''); }

// ── Escher ────────────────────────────────────────────────────────────────────
let coreMap = null;
async function renderEscher(containerId, noteId, fluxes) {
  const note = noteId ? $(noteId) : null;
  if (!window.escher) { if (note) note.textContent = 'Escher unavailable.'; return; }
  try {
    if (!coreMap) coreMap = await (await fetch('fba/core_map.json')).json();
    $(containerId).innerHTML = '';
    window.escher.Builder(coreMap, null, null, window.escher.libs.d3_select('#' + containerId), {
      menu: 'zoom', scroll_behavior: 'zoom', fill_screen: false, never_ask_before_quit: true,
      reaction_data: fluxes, reaction_styles: ['color', 'size', 'text'],
      reaction_scale: [
        { type: 'min', color: '#d0d0d0', size: 6 }, { type: 'value', value: 0, color: '#d0d0d0', size: 6 },
        { type: 'median', color: '#5b8ff9', size: 14 }, { type: 'max', color: '#c0392b', size: 28 }],
    });
    if (note) note.textContent = 'Central-carbon metabolism (e_coli_core). Colour & thickness = |flux|; hover a reaction for its value. Reactions absent from this core map are not shown.';
  } catch (e) { if (note) note.textContent = 'Map error: ' + e.message; console.error(e); }
}

// ── CSV ───────────────────────────────────────────────────────────────────────
function downloadFluxCSV(model, fluxes, base) {
  const nameById = {}; model.reactions.forEach(r => nameById[r.id] = r.name || '');
  let csv = 'reaction_id,name,flux\n';
  for (const [id, v] of Object.entries(fluxes)) csv += `${id},"${(nameById[id] || '').replace(/"/g, '""')}",${v}\n`;
  saveCSV(csv, base);
}
function saveCSV(csv, base) {
  const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([csv], { type: 'text/csv' }));
  a.download = base + '.csv'; a.click();
}
function esc(s) { return String(s == null ? '' : s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c])); }

// ── Mode toggle & init ────────────────────────────────────────────────────────
function setMode(mode) {
  S.mode = mode;
  document.querySelectorAll('#fba-modetoggle button').forEach(b => b.classList.toggle('active', b.dataset.mode === mode));
  const wrap = $('fba-conditions');
  wrap.classList.toggle('compare', mode === 'compare');
  S.conditions.b.card.style.display = mode === 'compare' ? '' : 'none';
}

async function init() {
  await loadPresets();
  const wrap = $('fba-conditions');
  wrap.appendChild(buildConditionCard('a'));
  wrap.appendChild(buildConditionCard('b'));
  S.conditions.b.card.style.display = 'none';
  document.querySelectorAll('#fba-modetoggle button').forEach(b => b.addEventListener('click', () => setMode(b.dataset.mode)));
  $('fba-run').addEventListener('click', run);

  // handoff from the GEM detail modal / table
  window.fbaSetModel = (slot, gemFile) => {
    slot = (slot === 'b') ? 'b' : 'a';
    if (slot === 'b') setMode('compare');
    selectModel(S.conditions[slot], gemFile);
    $('fba-section').scrollIntoView({ behavior: 'smooth' });
  };
  window.fbaSelectModel = (gemFile) => window.fbaSetModel('a', gemFile); // back-compat
}
if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
else init();
