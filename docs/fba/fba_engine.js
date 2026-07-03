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
