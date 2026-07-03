// EcopanGEM in-browser FBA / pFBA engine (ES module).
// Solves genome-scale COBRA-JSON models entirely client-side with glpk.js (WASM).
// Validated against COBRApy 0.27 (identical growth & pFBA sum|v|).
import GLPK from '../vendor/glpk.esm.js';

let _glpk = null;
export async function getGLPK() {
  if (!_glpk) _glpk = await GLPK();
  return _glpk;
}

// Build the LP: max c·v  s.t.  S·v = 0,  lb <= v <= ub.
// mediaBounds (optional): {EX_id: lower_bound}. When given, every exchange not
// listed is closed (lb=0), reproducing the "defined medium" convention.
// knockouts (optional): Set/array of reaction ids forced to lb=ub=0.
function buildLP(glpk, model, mediaBounds, knockouts) {
  const koSet = knockouts ? (knockouts instanceof Set ? knockouts : new Set(knockouts)) : null;
  const metRows = {};
  for (const met of model.metabolites)
    metRows[met.id] = { name: met.id, vars: [], bnds: { type: glpk.GLP_FX, ub: 0, lb: 0 } };

  const objVars = [];
  const bounds = [];
  for (const rxn of model.reactions) {
    let lb = rxn.lower_bound, ub = rxn.upper_bound;
    if (mediaBounds && rxn.id.startsWith('EX_'))
      lb = Object.prototype.hasOwnProperty.call(mediaBounds, rxn.id) ? mediaBounds[rxn.id] : 0;
    if (koSet && koSet.has(rxn.id)) { lb = 0; ub = 0; }
    let type;
    if (lb === ub) type = glpk.GLP_FX;
    else if (lb <= -1e30 && ub >= 1e30) type = glpk.GLP_FR;
    else if (lb <= -1e30) type = glpk.GLP_UP;
    else if (ub >= 1e30) type = glpk.GLP_LO;
    else type = glpk.GLP_DB;
    bounds.push({ name: rxn.id, type, ub, lb });
    for (const [met, coef] of Object.entries(rxn.metabolites))
      if (metRows[met]) metRows[met].vars.push({ name: rxn.id, coef });
    const oc = rxn.objective_coefficient || 0;
    if (oc !== 0) objVars.push({ name: rxn.id, coef: oc });
  }
  return {
    name: 'fba',
    objective: { direction: glpk.GLP_MAX, name: 'growth', vars: objVars },
    subjectTo: Object.values(metRows),
    bounds,
  };
}

export async function runFBA(model, mediaBounds, opts) {
  const glpk = await getGLPK();
  const lp = buildLP(glpk, model, mediaBounds, opts && opts.knockouts);
  const res = await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true });
  const r = res.result;
  return {
    status: r.status, optimal: r.status === glpk.GLP_OPT,
    growth: r.z, objectiveId: (lp.objective.vars[0] || {}).name, fluxes: r.vars || {},
  };
}

// Parsimonious FBA: fix biomass at the FBA optimum, then minimize total flux
// sum|v|. |v_i| is linearized with aux vars a_i >= v_i, a_i >= -v_i (min sum a_i).
export async function runPFBA(model, mediaBounds, fbaResult, opts) {
  const glpk = await getGLPK();
  const fba = fbaResult || (await runFBA(model, mediaBounds, opts));
  if (!fba.optimal || !(fba.growth > 1e-9)) return { ...fba, pfba: false };

  const lp = buildLP(glpk, model, mediaBounds, opts && opts.knockouts);
  const objId = fba.objectiveId;
  const gTarget = fba.growth * (1 - 1e-6);
  for (const b of lp.bounds)
    if (b.name === objId) { b.lb = gTarget; b.ub = Math.max(b.ub, fba.growth); b.type = glpk.GLP_DB; }

  const absVars = [], extraCons = [];
  for (const rxn of model.reactions) {
    const a = 'a_' + rxn.id;
    absVars.push({ name: a, coef: 1 });
    lp.bounds.push({ name: a, type: glpk.GLP_LO, lb: 0, ub: 1e30 });
    extraCons.push({ name: 'ap_' + rxn.id, vars: [{ name: a, coef: 1 }, { name: rxn.id, coef: -1 }], bnds: { type: glpk.GLP_LO, lb: 0, ub: 0 } });
    extraCons.push({ name: 'an_' + rxn.id, vars: [{ name: a, coef: 1 }, { name: rxn.id, coef: 1 }], bnds: { type: glpk.GLP_LO, lb: 0, ub: 0 } });
  }
  lp.subjectTo = lp.subjectTo.concat(extraCons);
  lp.objective = { direction: glpk.GLP_MIN, name: 'total_flux', vars: absVars };

  const res = await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true });
  const r = res.result;
  const all = r.vars || {};
  const fluxes = {};
  for (const rxn of model.reactions) fluxes[rxn.id] = all[rxn.id] || 0;
  return {
    status: r.status, optimal: r.status === glpk.GLP_OPT,
    growth: fluxes[objId], objectiveId: objId, totalFlux: r.z, fluxes, pfba: true,
  };
}

// Flux Variability Analysis: min & max flux for each reaction in reactionIds,
// subject to biomass >= fraction * optimum. opts: {fraction, knockouts, onProgress}.
export async function runFVA(model, mediaBounds, reactionIds, opts = {}) {
  const glpk = await getGLPK();
  const fraction = opts.fraction != null ? opts.fraction : 1.0;
  const lp = buildLP(glpk, model, mediaBounds, opts.knockouts);
  const objId = (lp.objective.vars[0] || {}).name;
  const fba = (await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true })).result;
  if (fba.status !== glpk.GLP_OPT || !(fba.z > 1e-9)) return { optimal: false, z: fba.z || 0, ranges: {} };
  // fix biomass >= fraction * optimum
  for (const b of lp.bounds) if (b.name === objId) { b.lb = fraction * fba.z; b.ub = Math.max(b.ub, fba.z); b.type = glpk.GLP_DB; }
  const ranges = {};
  let done = 0;
  for (const rid of reactionIds) {
    lp.objective = { direction: glpk.GLP_MIN, name: 'fva', vars: [{ name: rid, coef: 1 }] };
    const mn = (await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true })).result.z;
    lp.objective.direction = glpk.GLP_MAX;
    const mx = (await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true })).result.z;
    ranges[rid] = { min: mn, max: mx };
    if (opts.onProgress) opts.onProgress(++done, reactionIds.length);
  }
  return { optimal: true, z: fba.z, fraction, ranges };
}

// Dynamic (batch) FBA by static optimization. Substrate uptake follows
// Michaelis-Menten; biomass and tracked metabolite concentrations integrate by
// forward Euler. opts: {substrateEx, substrate0, biomass0, vmax, km, dt, tmax, trackEx, knockouts, onProgress}.
export async function runDFBA(model, mediaBounds, opts = {}) {
  const o = Object.assign({ substrateEx: 'EX_glc__D_e', substrate0: 10, biomass0: 0.01,
    vmax: 10, km: 0.5, dt: 0.1, tmax: 15, trackEx: [] }, opts);
  const media = { ...mediaBounds };
  const conc = { [o.substrateEx]: o.substrate0 };
  for (const e of o.trackEx) if (!(e in conc)) conc[e] = 0;
  const series = { t: [], biomass: [], conc: {}, substrateEx: o.substrateEx };
  for (const e of Object.keys(conc)) series.conc[e] = [];
  let X = o.biomass0, t = 0;
  const nsteps = Math.ceil(o.tmax / o.dt);
  for (let i = 0; i <= nsteps; i++) {
    const S = Math.max(0, conc[o.substrateEx]);
    let vUp = o.vmax * S / (o.km + S);
    if (X > 0 && o.dt > 0) vUp = Math.min(vUp, S / (X * o.dt)); // never consume more than present
    media[o.substrateEx] = -vUp;
    const fba = await runFBA(model, media, { knockouts: o.knockouts });
    const mu = fba.optimal ? fba.growth : 0;
    series.t.push(+t.toFixed(4)); series.biomass.push(X);
    for (const e of Object.keys(conc)) series.conc[e].push(Math.max(0, conc[e]));
    for (const e of Object.keys(conc)) {
      const flux = fba.fluxes[e] || 0;
      conc[e] += flux * X * o.dt; if (conc[e] < 0) conc[e] = 0;
    }
    X += mu * X * o.dt; t += o.dt;
    if (opts.onProgress) opts.onProgress(i, nsteps);
    if (mu <= 1e-9 && S <= 1e-9) break;
  }
  return series;
}

// Production envelope: trade-off between biomass and a target product. Fix the
// product flux at a grid of values from 0 to its maximum and record the max &
// min biomass at each. opts: {points, knockouts}.
export async function productionEnvelope(model, media, productId, opts = {}) {
  const glpk = await getGLPK();
  const points = opts.points || 20;
  const lpP = buildLP(glpk, model, media, opts.knockouts);
  lpP.objective = { direction: glpk.GLP_MAX, name: 'p', vars: [{ name: productId, coef: 1 }] };
  const prodMax = (await glpk.solve(lpP, { msglev: glpk.GLP_MSG_OFF, presol: true })).result.z || 0;
  const out = [];
  for (let i = 0; i < points; i++) {
    const v = prodMax * i / (points - 1);
    const lp = buildLP(glpk, model, media, opts.knockouts);
    for (const b of lp.bounds) if (b.name === productId) { b.lb = v; b.ub = v; b.type = glpk.GLP_FX; }
    const mx = (await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true })).result;
    lp.objective = { ...lp.objective, direction: glpk.GLP_MIN };
    const mn = (await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true })).result;
    out.push({ product: v, growthMax: mx.status === glpk.GLP_OPT ? mx.z : 0, growthMin: mn.status === glpk.GLP_OPT ? mn.z : 0 });
    if (opts.onProgress) opts.onProgress(i + 1, points);
  }
  return { productId, prodMax, points: out };
}

// Phenotype phase plane: max biomass over a grid of two uptake capacities.
// opts: {n, xMax, yMax, knockouts, onProgress}. Returns {xs, ys, Z} (Z[j][i]).
export async function phasePlane(model, media, xId, yId, opts = {}) {
  const n = opts.n || 20;
  const xMax = opts.xMax != null ? opts.xMax : Math.abs(media[xId] != null ? media[xId] : 20);
  const yMax = opts.yMax != null ? opts.yMax : Math.abs(media[yId] != null ? media[yId] : 20);
  const xs = [], ys = [], Z = [];
  for (let i = 0; i < n; i++) xs.push(xMax * i / (n - 1));
  for (let j = 0; j < n; j++) ys.push(yMax * j / (n - 1));
  let done = 0;
  for (let j = 0; j < n; j++) {
    const row = [];
    for (let i = 0; i < n; i++) {
      const m2 = { ...media, [xId]: -xs[i], [yId]: -ys[j] };
      const fba = await runFBA(model, m2, { knockouts: opts.knockouts });
      row.push(fba.optimal ? fba.growth : 0);
      if (opts.onProgress) opts.onProgress(++done, n * n);
    }
    Z.push(row);
  }
  return { xId, yId, xs, ys, Z };
}

// Single-reaction deletion (essentiality) over a list of reactions.
// Returns [{id, name, subsystem, growth, ratio}]. opts: {knockouts, onProgress}.
export async function essentialityScan(model, mediaBounds, reactionIds, opts = {}) {
  const wt = await runFBA(model, mediaBounds, { knockouts: opts.knockouts });
  const wtG = wt.optimal ? wt.growth : 0;
  const nameById = {}, subById = {};
  model.reactions.forEach(r => { nameById[r.id] = r.name || ''; subById[r.id] = r.subsystem || ''; });
  const base = opts.knockouts ? [...opts.knockouts] : [];
  const out = [];
  let done = 0;
  for (const rid of reactionIds) {
    const fba = await runFBA(model, mediaBounds, { knockouts: base.concat(rid) });
    const g = fba.optimal ? fba.growth : 0;
    out.push({ id: rid, name: nameById[rid], subsystem: subById[rid], growth: g, ratio: wtG > 1e-9 ? g / wtG : 0 });
    if (opts.onProgress) opts.onProgress(++done, reactionIds.length);
  }
  return { wtGrowth: wtG, results: out };
}

// All exchange reactions in a model: {id, name, met, defaultLb, defaultUb}. For the media editor.
export function listExchanges(model) {
  return model.reactions
    .filter(r => r.id.startsWith('EX_'))
    .map(r => ({ id: r.id, name: r.name || r.id, met: Object.keys(r.metabolites || {})[0] || '',
                 defaultLb: r.lower_bound, defaultUb: r.upper_bound }))
    .sort((a, b) => a.id.localeCompare(b.id));
}

// All reactions (id + name) for the knockout search.
export function listReactions(model) {
  return model.reactions.map(r => ({ id: r.id, name: r.name || '', subsystem: r.subsystem || '' }));
}

// Split exchange fluxes into uptake (negative) and secretion (positive), sorted by magnitude.
export function exchangeReport(model, fluxes, tol = 1e-6) {
  const nameById = {};
  for (const r of model.reactions) nameById[r.id] = r.name || r.id;
  const uptake = [], secretion = [];
  for (const [id, v] of Object.entries(fluxes)) {
    if (!id.startsWith('EX_')) continue;
    if (v < -tol) uptake.push({ id, name: nameById[id], flux: v });
    else if (v > tol) secretion.push({ id, name: nameById[id], flux: v });
  }
  uptake.sort((a, b) => a.flux - b.flux);
  secretion.sort((a, b) => b.flux - a.flux);
  return { uptake, secretion };
}
