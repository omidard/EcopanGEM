/* ── GEM Browser: theme ───────────────────────────────────────────────────────
   Small on purpose. The charts do not need to know about the theme: a chart box
   is white paper in both light and dark (--canvas), exactly as in Flux Studio, so
   a figure looks the same wherever it is read and the same wherever it is printed.
   Only the page around it changes.

   The initial theme is applied by an inline snippet in <head>, before first paint,
   so the page never flashes the wrong one. This file only wires the toggle.
   ────────────────────────────────────────────────────────────────────────────── */
(function () {
  'use strict';
  var root = document.documentElement;
  var KEY = 'gem-browser-theme';

  function apply(mode) {
    root.setAttribute('data-theme', mode);
    document.querySelectorAll('.theme-toggle').forEach(function (b) {
      b.setAttribute('aria-label', 'Switch to ' + (mode === 'dark' ? 'light' : 'dark') + ' theme');
      b.setAttribute('aria-pressed', String(mode === 'dark'));
    });
  }

  function wire() {
    apply(root.getAttribute('data-theme') || 'light');
    document.querySelectorAll('.theme-toggle').forEach(function (btn) {
      btn.addEventListener('click', function () {
        var next = root.getAttribute('data-theme') === 'dark' ? 'light' : 'dark';
        apply(next);
        try { localStorage.setItem(KEY, next); } catch (e) { /* private mode */ }
      });
    });
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', wire);
  else wire();
})();
