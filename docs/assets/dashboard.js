/* EcopanGEM landing dashboard — model-size distribution, genomes per phylogroup, isolation world map.
   Globals: Chart (chart.js) + d3 (v7). Data: assets/dashboard.json + assets/world.geojson.
   Interactive: click a phylogroup bar or a map country to filter the table below. */
(function () {
  'use strict';
  const PRIMARY = '#2563EB';
  const INK = '#0F1B2D', MUTED = '#64748B', GRID = '#E3E8EF';
  const FONT = "'Fira Sans', system-ui, -apple-system, 'Segoe UI', sans-serif";
  const PHYLO_COLORS = { 'A':'#e41a1c','B1':'#377eb8','B2':'#4daf4a','C':'#984ea3','D':'#ff7f00',
    'E':'#a65628','F':'#f781bf','G':'#999999','U':'#666666','cryptic':'#333333' };
  const phyColor = n => PHYLO_COLORS[n] || '#8aa0b6';

  Promise.all([
    fetch('assets/dashboard.json').then(r => r.json()),
    fetch('assets/world.geojson').then(r => r.json())
  ]).then(([db, world]) => {
    if (window.Chart) { Chart.defaults.font.family = FONT; Chart.defaults.color = MUTED; Chart.defaults.font.size = 12; }
    renderSize(db); renderPhylo(db); renderMap(db, world);
  }).catch(e => { console.error('dashboard load failed', e); const s = document.getElementById('dashboard-section'); if (s) s.style.display = 'none'; });

  // ---- 1. Model size distribution (metric toggle) ----
  let sizeChart = null;
  function renderSize(db) {
    const ctx = document.getElementById('chart-size'); if (!ctx || !window.Chart) return;
    const LBL = { n_reactions: 'Reactions', n_genes: 'Genes', n_metabolites: 'Metabolites' };
    function build(metric) {
      const h = db.size[metric];
      return { labels: h.labels, datasets: [{ label: 'E. coli GEMs', data: h.series.EcopanGEM,
        backgroundColor: PRIMARY + 'cc', borderColor: PRIMARY, borderWidth: 1, borderRadius: 3,
        categoryPercentage: 0.98, barPercentage: 0.96 }] };
    }
    function draw(metric) {
      const data = build(metric);
      if (sizeChart) { sizeChart.data = data; sizeChart.options.scales.x.title.text = LBL[metric] + ' per model'; sizeChart.update(); return; }
      sizeChart = new Chart(ctx, { type: 'bar', data, options: {
        responsive: true, maintainAspectRatio: false,
        plugins: { legend: { display: false }, tooltip: { callbacks: { title: it => LBL[metric] + ': ' + it[0].label, label: c => `${c.parsed.y.toLocaleString()} models` } } },
        scales: { x: { title: { display: true, text: 'Reactions per model', color: MUTED, font: { size: 11 } }, grid: { display: false }, ticks: { maxRotation: 60, minRotation: 60, autoSkip: true, maxTicksLimit: 14, font: { size: 10 } } },
                  y: { title: { display: true, text: 'Models', color: MUTED, font: { size: 11 } }, grid: { color: GRID }, ticks: { precision: 0 } } } } });
    }
    draw('n_reactions');
    document.querySelectorAll('#size-toggle .metric-btn').forEach(btn => btn.addEventListener('click', () => {
      document.querySelectorAll('#size-toggle .metric-btn').forEach(b => b.classList.remove('active'));
      btn.classList.add('active'); draw(btn.dataset.metric);
    }));
  }

  // ---- 2. Genomes per phylogroup (click to filter) ----
  function renderPhylo(db) {
    const ctx = document.getElementById('chart-phylo'); if (!ctx || !window.Chart) return;
    const rows = db.phylogroups;
    const chart = new Chart(ctx, { type: 'bar',
      data: { labels: rows.map(r => r.name), datasets: [{ label: 'Genomes', data: rows.map(r => r.count),
        backgroundColor: rows.map(r => phyColor(r.name) + 'dd'), borderColor: rows.map(r => phyColor(r.name)), borderWidth: 1,
        borderRadius: 3, barPercentage: 0.82, categoryPercentage: 0.82 }] },
      options: { indexAxis: 'y', responsive: true, maintainAspectRatio: false,
        onClick: (ev, els) => { if (!els.length) return; const name = rows[els[0].index].name; filterTable('filter-phylo', name); },
        plugins: { legend: { display: false }, tooltip: { callbacks: { label: c => `${c.parsed.x.toLocaleString()} genomes · click to filter` } } },
        scales: { x: { title: { display: true, text: 'Genome-scale models', color: MUTED, font: { size: 11 } }, grid: { color: GRID }, ticks: { precision: 0 } },
                  y: { grid: { display: false }, ticks: { font: { size: 11, weight: '600' }, color: INK, autoSkip: false } } } } });
    ctx.style.cursor = 'pointer';
  }

  // ---- 3. Isolation world map (click to filter) ----
  function renderMap(db, world) {
    const host = document.getElementById('map'); if (!host || !window.d3) return;
    const geo = db.geo, max = db.geo_max || 1;
    const color = c => c > 0 ? d3.interpolateBlues(0.18 + 0.82 * Math.sqrt(c / max)) : '#EEF2F7';
    let tip = document.getElementById('map-tip'); if (!tip) { tip = document.createElement('div'); tip.id = 'map-tip'; tip.className = 'map-tip'; document.body.appendChild(tip); }
    function paint() {
      const w = host.clientWidth || 900, h = Math.max(320, Math.round(w * 0.5));
      host.innerHTML = '';
      const svg = d3.select(host).append('svg').attr('width', w).attr('height', h).attr('viewBox', `0 0 ${w} ${h}`).style('display', 'block');
      const fc = { type: 'FeatureCollection', features: world.features.filter(f => f.properties.iso !== 'ATA') };
      const proj = d3.geoNaturalEarth1().fitSize([w, h - 8], fc);
      const path = d3.geoPath(proj);
      svg.append('g').selectAll('path').data(fc.features).join('path')
        .attr('d', path).attr('fill', d => color(geo[d.properties.iso] || 0))
        .attr('stroke', '#fff').attr('stroke-width', 0.4)
        .style('cursor', d => geo[d.properties.iso] ? 'pointer' : 'default')
        .on('mousemove', function (ev, d) {
          const n = geo[d.properties.iso] || 0;
          tip.innerHTML = `<strong>${d.properties.name}</strong><br>${n ? n.toLocaleString() + ' model' + (n === 1 ? '' : 's') : 'no isolates'}`;
          tip.style.opacity = 1; tip.style.left = (ev.pageX + 14) + 'px'; tip.style.top = (ev.pageY - 10) + 'px';
          d3.select(this).attr('stroke', '#0F1B2D').attr('stroke-width', 1.2).raise();
        })
        .on('mouseleave', function () { tip.style.opacity = 0; d3.select(this).attr('stroke', '#fff').attr('stroke-width', 0.4); })
        .on('click', function (ev, d) { const n = geo[d.properties.iso] || 0; if (!n) return; tip.style.opacity = 0; filterTable('filter-country', d.properties.name); });
    }
    paint();
    let t; window.addEventListener('resize', () => { clearTimeout(t); t = setTimeout(paint, 200); });
    const leg = document.getElementById('map-legend');
    if (leg) { const stops = []; for (let i = 0; i <= 6; i++) stops.push(color(Math.round(max * (i / 6) * (i / 6)))); leg.querySelector('.grad').style.background = `linear-gradient(90deg, ${['#EEF2F7'].concat(stops.slice(1)).join(',')})`; leg.querySelector('.max').textContent = max.toLocaleString(); }
  }

  // ---- filter the DataTables table below (matches an <option> value, dispatches change) ----
  function filterTable(selectId, value) {
    const sel = document.getElementById(selectId);
    if (sel && [...sel.options].some(o => o.value === value)) {
      sel.value = value; sel.dispatchEvent(new Event('change'));
      const tw = document.querySelector('.table-wrap'); if (tw) tw.scrollIntoView({ behavior: 'smooth', block: 'start' });
    }
  }
})();
