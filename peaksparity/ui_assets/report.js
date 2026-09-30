/* peaks-parity results viewer.  Self-contained: no external libraries, no network access.
   Data come from two <script type="application/json"> elements: #pp-data (tables, summary,
   parameters) and #pp-spectra (peak lists, parsed lazily on first use). */
(function () {
  'use strict';
  // Show script errors on the page instead of failing silently.
  window.addEventListener('error', e => {
    const d = document.createElement('div');
    d.className = 'pp-error';
    d.textContent = 'Viewer error: ' + (e.message || e) + (e.lineno ? ' (line ' + e.lineno + ')' : '');
    document.body.insertBefore(d, document.body.firstChild);
  });

  // ------------------------------------------------------------------ data
  const D = JSON.parse(document.getElementById('pp-data').textContent);
  let SPEC = null;
  function spectra() {
    if (SPEC === null) {
      const el = document.getElementById('pp-spectra');
      SPEC = el ? JSON.parse(el.textContent) : {};
    }
    return SPEC;
  }
  const T = D.tables || {};

  // ------------------------------------------------------------------ masses
  const PROTON = 1.007276467, H2O = 18.010565, NH3 = 17.026549;
  const AA = { G: 57.021464, A: 71.037114, S: 87.032028, P: 97.052764, V: 99.068414, T: 101.047679,
    C: 103.009185, L: 113.084064, I: 113.084064, J: 113.084064, N: 114.042927, D: 115.026943,
    Q: 128.058578, K: 128.094963, E: 129.042593, M: 131.040485, H: 137.058912, F: 147.068414,
    R: 156.101111, Y: 163.063329, W: 186.079313, U: 150.953636, O: 237.147727, X: 0 };

  // ------------------------------------------------------------------ helpers
  function esc(s) {
    return String(s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  }
  function h(tag, attrs, ...kids) {
    const e = document.createElement(tag);
    if (attrs) for (const k in attrs) {
      const v = attrs[k];
      if (v === null || v === undefined || v === false) continue;
      if (k === 'class') e.className = v;
      else if (k === 'html') e.innerHTML = v;
      else if (k === 'text') e.textContent = v;
      else if (k.startsWith('on')) e.addEventListener(k.slice(2), v);
      else e.setAttribute(k, v === true ? '' : v);
    }
    for (const k of kids.flat()) if (k !== null && k !== undefined && k !== false)
      e.appendChild(typeof k === 'string' || typeof k === 'number' ? document.createTextNode(String(k)) : k);
    return e;
  }
  function sci(v) {
    if (v === 0) return '0';
    let e = Math.floor(Math.log10(Math.abs(v)));
    let m = Math.round(v / Math.pow(10, e) * 100) / 100;
    if (Math.abs(m) >= 10) { m /= 10; e += 1; }
    return String(+m.toFixed(2)) + 'E' + e;
  }
  function fmt(v, kind) {
    if (v === null || v === undefined) return '';
    if (typeof v === 'boolean') return v ? 'true' : 'false';
    if (typeof v !== 'number') return String(v);
    if (!isFinite(v)) return String(v);
    if (kind && kind.startsWith('fixed')) return v.toFixed(+kind.slice(5));
    if (kind === 'sci') return sci(v);
    if (kind === 'int') return String(Math.round(v));
    if (kind === 'dec2') return String(Math.round(v * 100) / 100).replace(/^(-?)0\./, '$1.');
    return String(v);
  }
  const nf = new Intl.NumberFormat('en-US');
  function num(v) { return v === null || v === undefined ? '–' : (typeof v === 'number' ? nf.format(v) : String(v)); }
  function pct(x) { return x === null || x === undefined ? 'n/a' : (100 * x).toFixed(2) + '%'; }
  function range(n) { const a = new Array(n); for (let i = 0; i < n; i++) a[i] = i; return a; }
  function col(t, name) { return t ? t.columns.indexOf(name) : -1; }
  function colLike(t, re) { if (!t) return -1; for (let i = 0; i < t.columns.length; i++) if (re.test(t.columns[i])) return i; return -1; }
  const SRC_EXT = /(\.(raw|mzml|mzxml|mgf|d|wiff|ms2|gz|bz2|zip))+$/i;
  function normSource(s) { const b = String(s || '').replace(/\\/g, '/').split('/').pop(); return b.replace(SRC_EXT, ''); }
  function stripped(p) { const i = D.pepinfo[p]; return i ? i[0] : String(p).replace(/\([^)]*\)|\[[^\]]*\]/g, '').replace(/[^A-Z]/g, ''); }
  function debounce(fn, ms) { let t = null; return function (...a) { clearTimeout(t); t = setTimeout(() => fn.apply(this, a), ms); }; }
  function pepHtml(p) {
    return esc(p).replace(/(\([+-]?[\d.]+\)|\[[^\]]*\])/g, '<span class="pep-mod">$1</span>');
  }

  // ------------------------------------------------------------------ indexes
  const IDX = {};
  (function buildIndexes() {
    const ps = T.psms;
    IDX.psmByPep = new Map();
    if (ps) {
      const ip = col(ps, 'Peptide');
      ps.rows.forEach((r, i) => { const k = r[ip]; if (!IDX.psmByPep.has(k)) IDX.psmByPep.set(k, []); IDX.psmByPep.get(k).push(i); });
    }
    const pp = T.peptides;
    IDX.pepRow = new Map();
    if (pp) { const ip = col(pp, 'Peptide'); pp.rows.forEach((r, i) => { if (!IDX.pepRow.has(r[ip])) IDX.pepRow.set(r[ip], i); }); }
    const pr = T.proteins;
    IDX.protByAcc = new Map(); IDX.protByGroup = new Map();
    if (pr) {
      const ia = col(pr, 'Accession'), ig = col(pr, 'Protein Group');
      pr.rows.forEach((r, i) => {
        IDX.protByAcc.set(r[ia], i);
        const g = ig >= 0 ? r[ig] : i;
        if (!IDX.protByGroup.has(g)) IDX.protByGroup.set(g, []);
        IDX.protByGroup.get(g).push(i);
      });
    }
    const x = T.protein_peptides;
    IDX.ppByAcc = new Map();
    if (x) { const ia = col(x, 'Accession'); x.rows.forEach((r, i) => { if (!IDX.ppByAcc.has(r[ia])) IDX.ppByAcc.set(r[ia], []); IDX.ppByAcc.get(r[ia]).push(i); }); }
  })();

  function psmObj(i) {
    const t = T.psms, r = t.rows[i];
    const g = n => { const j = col(t, n); return j >= 0 ? r[j] : null; };
    return { row: i, peptide: g('Peptide'), z: g('z'), mz: g('m/z'), scan: g('Scan'), source: g('Source File'),
      score: g('-10LgP'), ppm: g('ppm'), rt: g('RT'), rawrt: g('raw RT'), mass: g('Mass'), acc: g('Accession'),
      feature: g('Feature Id'), ptm: g('PTM') };
  }
  function specKey(src, scan) { return normSource(src) + '|' + scan; }
  function hasSpectrum(src, scan) { return Object.prototype.hasOwnProperty.call(spectra(), specKey(src, scan)); }
  /** PSM row index of the best PSM of a peptide (the peptide table's representative scan). */
  function bestPsm(pep) {
    const list = IDX.psmByPep.get(pep) || [];
    if (!list.length) return -1;
    const ps = T.psms, isc = col(ps, 'Scan'), isr = col(ps, 'Source File'), is = col(ps, '-10LgP');
    const pr = IDX.pepRow.get(pep);
    if (pr !== undefined) {
      const pt = T.peptides, r = pt.rows[pr], sc = r[col(pt, 'Scan')], src = normSource(r[col(pt, 'Source File')]);
      const hit = list.filter(i => ps.rows[i][isc] === sc && normSource(ps.rows[i][isr]) === src);
      if (hit.length) return hit.reduce((a, b) => (ps.rows[b][is] > ps.rows[a][is] ? b : a));
    }
    let best = list[0];
    for (const i of list) {
      const a = hasSpectrum(ps.rows[i][isr], ps.rows[i][isc]), b = hasSpectrum(ps.rows[best][isr], ps.rows[best][isc]);
      if ((a && !b) || (a === b && ps.rows[i][is] > ps.rows[best][is])) best = i;
    }
    return best;
  }

  // ================================================================== DataTable
  const TEXT_COLS = ['Peptide', 'Accession', 'Gene', 'Description', 'peptide', 'modseq', 'accessions',
    'top_accessions', 'accession', 'stripped', 'category', 'Extended Peptide', 'PTM'];
  const NUMERIC = new Set(['int', 'jfloat', 'sci', 'dec2', 'fixed1', 'fixed2', 'fixed3', 'fixed4', 'fixed5', 'fixed6', 'num']);

  class DataTable {
    constructor(host, o) {
      this.o = o; this.t = o.table; this.host = host;
      this.base = o.rowsIdx ? o.rowsIdx.slice() : range(this.t.rows.length);
      this.pageSize = o.pageSize || 50; this.page = 0; this.text = ''; this.minScore = null; this.selected = -1;
      this.scoreCol = col(this.t, o.scoreCol === undefined ? '-10LgP' : o.scoreCol);
      this.textIdx = (o.textCols || TEXT_COLS).map(n => col(this.t, n)).filter(i => i >= 0);
      this.pepIdx = col(this.t, 'Peptide') >= 0 ? col(this.t, 'Peptide') : col(this.t, 'peptide');
      const hide = new Set(o.hidden || []);
      this.visible = this.t.columns.map(c => !hide.has(c));
      this.sortCol = o.sort ? col(this.t, o.sort[0]) : -1; this.sortDesc = o.sort ? !!o.sort[1] : false;
      this.hay = null; this.build(); this.refresh();
    }
    isNum(j) { const k = this.t.kinds[j]; if (NUMERIC.has(k)) return true; if (k !== 'str') return false;
      for (const i of this.base.slice(0, 50)) { const v = this.t.rows[i][j]; if (v !== null && v !== '') return typeof v === 'number'; } return false; }
    build() {
      const o = this.o, self = this;
      this.host.innerHTML = '';
      const root = h('div', { class: 'dt' });
      const bar = h('div', { class: 'dt-bar' });
      this.q = h('input', { type: 'search', placeholder: o.placeholder || 'Filter: sequence, accession, gene…',
        oninput: debounce(() => { self.text = self.q.value.trim().toLowerCase(); self.page = 0; self.refresh(); }, 180) });
      bar.appendChild(this.q);
      if (this.scoreCol >= 0) {
        this.ms = h('input', { type: 'number', step: '1', placeholder: 'min',
          oninput: debounce(() => { const v = parseFloat(self.ms.value); self.minScore = isNaN(v) ? null : v; self.page = 0; self.refresh(); }, 180) });
        bar.appendChild(h('label', null, this.t.columns[this.scoreCol], ' ≥ ', this.ms));
      }
      if (o.controls) o.controls(bar, this);
      const ps = h('select', { onchange: () => { self.pageSize = +ps.value; self.page = 0; self.render(); } },
        [25, 50, 100, 250, 1000].map(n => h('option', { value: n, selected: n === self.pageSize }, n + ' / page')));
      bar.appendChild(ps);
      const colsBox = h('div');
      this.t.columns.forEach((c, j) => colsBox.appendChild(h('label', null,
        h('input', { type: 'checkbox', checked: this.visible[j], onchange: e => { self.visible[j] = e.target.checked; self.render(); } }), ' ', c)));
      bar.appendChild(h('details', { class: 'cols' }, h('summary', { class: 'btn' }, 'Columns'), colsBox));
      bar.appendChild(h('button', { class: 'btn', title: 'Download the filtered rows as CSV', onclick: () => self.download() }, 'CSV'));
      this.count = h('span', { class: 'count' });
      bar.appendChild(this.count);
      root.appendChild(bar);
      this.table = h('table');
      this.thead = h('thead'); this.tbody = h('tbody');
      this.table.appendChild(this.thead); this.table.appendChild(this.tbody);
      this.wrap = h('div', { class: 'dt-wrap' }, this.table);
      if (o.maxHeight) this.wrap.style.maxHeight = o.maxHeight;
      root.appendChild(this.wrap);
      this.pager = h('div', { class: 'dt-pager' });
      root.appendChild(this.pager);
      this.host.appendChild(root);
      this.thead.addEventListener('click', e => {
        const th = e.target.closest('th'); if (!th) return;
        const j = +th.dataset.j;
        if (this.sortCol === j) this.sortDesc = !this.sortDesc; else { this.sortCol = j; this.sortDesc = this.isNum(j); }
        this.sort(); this.page = 0; this.render();
      });
      this.tbody.addEventListener('click', e => {
        const tr = e.target.closest('tr'); if (!tr || tr.dataset.i === undefined) return;
        this.selected = +tr.dataset.i;
        for (const x of this.tbody.querySelectorAll('tr.sel')) x.classList.remove('sel');
        tr.classList.add('sel');
        if (o.onRowClick) o.onRowClick(this.selected, this);
      });
    }
    haystack(i) {
      if (!this.hay) this.hay = new Array(this.t.rows.length);
      let s = this.hay[i];
      if (s === undefined) {
        const r = this.t.rows[i];
        s = this.textIdx.map(j => (r[j] === null ? '' : String(r[j]))).join('\u0001');
        if (this.pepIdx >= 0 && r[this.pepIdx]) s += '\u0001' + stripped(r[this.pepIdx]);
        s = s.toLowerCase(); this.hay[i] = s;
      }
      return s;
    }
    setBase(idx) { this.base = idx.slice(); this.page = 0; this.refresh(); }
    setText(txt) { this.q.value = txt; this.text = txt.trim().toLowerCase(); this.page = 0; this.refresh(); }
    refresh() {
      const txt = this.text, ms = this.minScore, sc = this.scoreCol, rows = this.t.rows, ex = this.o.filter;
      const terms = txt ? txt.split(/\s+/) : [];
      this.view = this.base.filter(i => {
        if (ms !== null && sc >= 0) { const v = rows[i][sc]; if (v === null || v < ms) return false; }
        if (ex && !ex(i)) return false;
        if (terms.length) { const hs = this.haystack(i); for (const t of terms) if (hs.indexOf(t) < 0) return false; }
        return true;
      });
      this.sort(); this.render();
    }
    sort() {
      const j = this.sortCol; if (j < 0) return;
      const rows = this.t.rows, d = this.sortDesc ? -1 : 1;
      const numeric = this.isNum(j);
      this.view.sort((a, b) => {
        let x = rows[a][j], y = rows[b][j];
        if (x === null || x === '') return (y === null || y === '') ? a - b : 1;
        if (y === null || y === '') return -1;
        if (numeric) { if (typeof x !== 'number') x = parseFloat(x); if (typeof y !== 'number') y = parseFloat(y); }
        else if (typeof x === 'boolean') { x = x ? 1 : 0; y = y ? 1 : 0; }
        else { x = String(x); y = String(y); }
        return x < y ? -d : x > y ? d : a - b;
      });
    }
    render() {
      const t = this.t, o = this.o, vis = this.visible;
      const cols = t.columns.map((c, j) => j).filter(j => vis[j]);
      const numCls = j => (this.isNumCache || (this.isNumCache = {}))[j] !== undefined ? this.isNumCache[j] : (this.isNumCache[j] = this.isNum(j));
      this.thead.innerHTML = '<tr>' + cols.map(j => `<th data-j="${j}" class="${numCls(j) ? 'num' : ''}" title="Sort by ${esc(t.columns[j])}">${esc(t.columns[j])}${j === this.sortCol ? `<span class="arr">${this.sortDesc ? '▼' : '▲'}</span>` : ''}</th>`).join('') + '</tr>';
      const n = this.view.length, np = Math.max(1, Math.ceil(n / this.pageSize));
      if (this.page >= np) this.page = np - 1;
      const a = this.page * this.pageSize, b = Math.min(n, a + this.pageSize);
      const out = [];
      for (let k = a; k < b; k++) {
        const i = this.view[k], r = t.rows[i];
        const cls = (i === this.selected ? 'sel ' : '') + (o.rowClass ? o.rowClass(i) : '');
        out.push(`<tr data-i="${i}" class="${cls}">`);
        for (const j of cols) {
          const v = r[j], name = t.columns[j];
          let s;
          if (o.cell && o.cell[name]) s = o.cell[name](v, i);
          else if (name === 'In ref' && typeof v === 'boolean') s = v ? '<span class="flag-yes">yes</span>' : '<span class="flag-no">no</span>';
          else if (name === 'Peptide' || name === 'peptide' || name === 'modseq') s = v ? pepHtml(v) : '';
          else s = esc(fmt(v, t.kinds[j]));
          const title = typeof v === 'string' && v.length > 40 ? ` title="${esc(v)}"` : '';
          out.push(`<td class="${numCls(j) ? 'num' : ''}"${title}>${s}</td>`);
        }
        out.push('</tr>');
      }
      this.tbody.innerHTML = out.join('') || `<tr><td colspan="${cols.length}" class="muted">No rows</td></tr>`;
      this.count.textContent = n === this.base.length ? `${num(n)} rows` : `${num(n)} of ${num(this.base.length)} rows`;
      const self = this;
      this.pager.innerHTML = '';
      const go = p => { self.page = Math.max(0, Math.min(np - 1, p)); self.render(); self.wrap.scrollTop = 0; };
      const pin = h('input', { type: 'number', min: 1, max: np, value: this.page + 1, onchange: () => go((+pin.value || 1) - 1) });
      this.pager.append(
        h('button', { class: 'btn', disabled: this.page === 0, onclick: () => go(0) }, '«'),
        h('button', { class: 'btn', disabled: this.page === 0, onclick: () => go(self.page - 1) }, '‹ Prev'),
        h('span', null, 'Page ', pin, ` of ${num(np)}`),
        h('button', { class: 'btn', disabled: this.page >= np - 1, onclick: () => go(self.page + 1) }, 'Next ›'),
        h('button', { class: 'btn', disabled: this.page >= np - 1, onclick: () => go(np - 1) }, '»'),
        h('span', { class: 'muted' }, n ? `rows ${num(a + 1)}–${num(b)}` : ''));
    }
    selectRow(i) {
      this.selected = i;
      const k = this.view.indexOf(i);
      if (k >= 0) this.page = Math.floor(k / this.pageSize);
      this.render();
      const tr = this.tbody.querySelector(`tr[data-i="${i}"]`);
      if (tr && tr.scrollIntoView) tr.scrollIntoView({ block: 'nearest' });
    }
    download() {
      const t = this.t, cols = t.columns.map((c, j) => j).filter(j => this.visible[j]);
      const q = s => { s = s === null || s === undefined ? '' : String(s); return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; };
      const lines = [cols.map(j => q(t.columns[j])).join(',')];
      for (const i of this.view) lines.push(cols.map(j => q(fmt(t.rows[i][j], t.kinds[j]))).join(','));
      const blob = new Blob([lines.join('\n') + '\n'], { type: 'text/csv' });
      const a = h('a', { href: URL.createObjectURL(blob), download: (this.o.name || 'table') + '.csv' });
      document.body.appendChild(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
    }
  }

  // ================================================================== charts (SVG)
  function niceTicks(lo, hi, n) {
    if (!(hi > lo)) { hi = lo + 1; }
    const span = hi - lo, step0 = span / Math.max(1, n), mag = Math.pow(10, Math.floor(Math.log10(step0)));
    const err = step0 / mag, step = mag * (err >= 7.5 ? 10 : err >= 3.5 ? 5 : err >= 1.5 ? 2 : 1);
    const out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(+v.toFixed(10));
    return out;
  }
  function tickLabel(v) { const a = Math.abs(v); return a >= 1e5 ? sci(v) : a >= 100 || v === Math.round(v) ? String(Math.round(v * 100) / 100) : String(+v.toFixed(3)); }
  function axes(W, H, m, xlo, xhi, ylo, yhi, xlabel, ylabel) {
    const X = v => m.l + (v - xlo) / (xhi - xlo || 1) * (W - m.l - m.r);
    const Y = v => H - m.b - (v - ylo) / (yhi - ylo || 1) * (H - m.t - m.b);
    let s = '<g class="grid">';
    for (const v of niceTicks(ylo, yhi, 5)) s += `<line x1="${m.l}" x2="${W - m.r}" y1="${Y(v)}" y2="${Y(v)}"/>`;
    s += '</g><g class="axis">';
    s += `<line x1="${m.l}" x2="${W - m.r}" y1="${H - m.b}" y2="${H - m.b}"/><line x1="${m.l}" x2="${m.l}" y1="${m.t}" y2="${H - m.b}"/>`;
    for (const v of niceTicks(xlo, xhi, 8)) s += `<line x1="${X(v)}" x2="${X(v)}" y1="${H - m.b}" y2="${H - m.b + 4}"/><text class="lbl" x="${X(v)}" y="${H - m.b + 16}" text-anchor="middle">${tickLabel(v)}</text>`;
    for (const v of niceTicks(ylo, yhi, 5)) s += `<text class="lbl" x="${m.l - 6}" y="${Y(v) + 4}" text-anchor="end">${tickLabel(v)}</text>`;
    s += '</g>';
    if (xlabel) s += `<text class="lbl" x="${(m.l + W - m.r) / 2}" y="${H - 4}" text-anchor="middle">${esc(xlabel)}</text>`;
    if (ylabel) s += `<text class="lbl" transform="translate(12 ${(m.t + H - m.b) / 2}) rotate(-90)" text-anchor="middle">${esc(ylabel)}</text>`;
    return { s, X, Y };
  }
  function quantile(sorted, q) { if (!sorted.length) return NaN; const p = (sorted.length - 1) * q, i = Math.floor(p); return sorted[i] + (sorted[Math.min(i + 1, sorted.length - 1)] - sorted[i]) * (p - i); }
  function histogramSVG(values, o) {
    o = o || {};
    const v = values.filter(x => typeof x === 'number' && isFinite(x)).sort((a, b) => a - b);
    const W = 560, H = 250, m = { l: 50, r: 12, t: 12, b: 38 };
    if (!v.length) return `<p class="muted">No data.</p>`;
    let lo = o.clip ? quantile(v, o.clip) : v[0], hi = o.clip ? quantile(v, 1 - o.clip) : v[v.length - 1];
    if (o.symmetric) { const a = Math.max(Math.abs(lo), Math.abs(hi)); lo = -a; hi = a; }
    if (hi <= lo) { hi = lo + 1; }
    const nb = o.bins || 40, w = (hi - lo) / nb, cnt = new Array(nb).fill(0);
    let outside = 0;
    for (const x of v) { if (x < lo || x > hi) { outside++; continue; } cnt[Math.min(nb - 1, Math.floor((x - lo) / w))]++; }
    const ymax = Math.max(...cnt) * 1.08;
    const ax = axes(W, H, m, lo, hi, 0, ymax, o.xlabel, o.ylabel || 'count');
    let s = `<svg class="chart" viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(o.xlabel || 'histogram')}">` + ax.s;
    cnt.forEach((c, i) => { if (c) { const x0 = ax.X(lo + i * w), x1 = ax.X(lo + (i + 1) * w); s += `<rect class="bar" x="${x0 + 0.5}" y="${ax.Y(c)}" width="${Math.max(0.5, x1 - x0 - 1)}" height="${ax.Y(0) - ax.Y(c)}"><title>${tickLabel(lo + i * w)} – ${tickLabel(lo + (i + 1) * w)}: ${c}</title></rect>`; } });
    for (const L of o.lines || []) if (L.x >= lo && L.x <= hi) s += `<line class="thr" x1="${ax.X(L.x)}" x2="${ax.X(L.x)}" y1="${m.t}" y2="${H - m.b}"/><text class="lbl" x="${ax.X(L.x) + 4}" y="${m.t + 10}">${esc(L.label)}</text>`;
    s += '</svg>';
    const med = quantile(v, 0.5);
    s += `<div class="muted small">n = ${num(v.length)}; median ${med.toFixed(2)}${outside ? `; ${num(outside)} outside the plotted range` : ''}</div>`;
    return s;
  }
  function scatterSVG(xs, ys, o) {
    o = o || {};
    const W = 560, H = 250, m = { l: 50, r: 12, t: 12, b: 38 };
    const pts = [];
    for (let i = 0; i < xs.length; i++) if (typeof xs[i] === 'number' && typeof ys[i] === 'number' && isFinite(xs[i]) && isFinite(ys[i])) pts.push([xs[i], ys[i]]);
    if (!pts.length) return `<p class="muted">No data.</p>`;
    const sx = pts.map(p => p[0]).sort((a, b) => a - b), sy = pts.map(p => p[1]).sort((a, b) => a - b);
    let ylo = o.clip ? quantile(sy, o.clip) : sy[0], yhi = o.clip ? quantile(sy, 1 - o.clip) : sy[sy.length - 1];
    if (o.symmetric) { const a = Math.max(Math.abs(ylo), Math.abs(yhi), 1e-9); ylo = -a; yhi = a; }
    const xlo = sx[0], xhi = sx[sx.length - 1];
    const ax = axes(W, H, m, xlo, xhi, ylo, yhi, o.xlabel, o.ylabel);
    let d = '';
    for (const [x, y] of pts) { if (y < ylo || y > yhi) continue; d += `M${ax.X(x).toFixed(1)} ${ax.Y(y).toFixed(1)}h0`; }
    let s = `<svg class="chart" viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(o.ylabel || '')} vs ${esc(o.xlabel || '')}">` + ax.s;
    if (ylo < 0 && yhi > 0) s += `<line x1="${m.l}" x2="${W - m.r}" y1="${ax.Y(0)}" y2="${ax.Y(0)}" stroke="currentColor" opacity=".3"/>`;
    s += `<path d="${d}" stroke="var(--accent)" stroke-opacity=".45" stroke-width="3" stroke-linecap="round" fill="none"/></svg>`;
    return s;
  }

  // ================================================================== layout
  const app = document.getElementById('app');
  app.innerHTML = '';
  const meta = D.meta || {};
  const TABS = [['summary', 'Summary'], ['proteins', 'Proteins'], ['peptides', 'Peptides'], ['psms', 'PSMs']];
  const nav = h('nav', { class: 'tabs' });
  const header = h('header', { class: 'top' },
    h('div', null, h('div', { class: 'title' }, meta.title || 'Results'),
      h('div', { class: 'sub' }, [meta.sample ? `sample ${meta.sample}` : null, meta.generated ? `generated ${meta.generated}` : null, meta.tool].filter(Boolean).join(' · '))),
    nav);
  const mainEl = h('main');
  app.append(header, mainEl);
  const sections = {};
  for (const [id, label] of TABS) {
    sections[id] = h('section', { class: 'tab', id: 'tab-' + id });
    mainEl.appendChild(sections[id]);
    nav.appendChild(h('button', { 'data-tab': id, onclick: () => showTab(id, true) }, label));
  }
  const built = {};
  const builders = {};
  function showTab(id, push) {
    if (!sections[id]) id = 'summary';
    for (const k in sections) sections[k].classList.toggle('on', k === id);
    for (const b of nav.children) b.classList.toggle('on', b.dataset.tab === id);
    if (!built[id]) { built[id] = true; builders[id](sections[id]); }
    if (push) setHash({ tab: id });
  }
  function setHash(o) { try { history.replaceState(null, '', '#' + Object.keys(o).map(k => k + '=' + encodeURIComponent(o[k])).join('&')); } catch (e) { /* file:// in some browsers */ } }

  // modal (spectrum from the Proteins tab / concordance lists)
  const modalBox = h('div', { class: 'box' });
  const modal = h('div', { class: 'modal', onclick: e => { if (e.target === modal) closeModal(); } }, modalBox);
  document.body.appendChild(modal);
  function openModal(fill) {
    modalBox.innerHTML = '';
    modalBox.appendChild(h('button', { class: 'btn close', onclick: closeModal, title: 'Close (Esc)' }, 'Close'));
    const host = h('div'); modalBox.appendChild(host); fill(host);
    modal.classList.add('on');
  }
  function closeModal() { modal.classList.remove('on'); modalBox.innerHTML = ''; }
  document.addEventListener('keydown', e => { if (e.key === 'Escape') closeModal(); });
  function openSpectrumModal(psm) { openModal(host => new SpectrumView(host, psm)); }
  function openPeptideModal(pep) {
    const i = bestPsm(pep);
    if (i >= 0) openSpectrumModal(psmObj(i));
    else openModal(host => { host.appendChild(h('p', null, 'No PSM for ', h('span', { class: 'mono' }, pep), ' in this export.')); });
  }

  // ================================================================== Summary
  builders.summary = function (sec) {
    const c = D.counts || {}, st = D.statistics || {}, p = D.params || {};
    const tile = (k, v, x) => h('div', { class: 'tile' }, h('div', { class: 'v' }, num(v)), h('div', { class: 'k' }, k), x ? h('div', { class: 'x' }, x) : null);
    const tiles = h('div', { class: 'tiles' },
      tile('PSMs', c.psms !== undefined ? c.psms : c.psm_keys, `${num(c.psm_rows)} rows · ${num(c.psm_keys)} scan/peptide/z keys` + (st.psms !== undefined ? ` · summary ${num(st.psms)}` : '')),
      tile('Peptides', c.peptides, `${num(c.peptide_rows)} rows (I/L separate)` + (st.peptides !== undefined ? ` · summary ${num(st.peptides)}` : '')),
      tile('Sequences', c.sequences, 'I/L-collapsed, unmodified' + (st.sequences !== undefined ? ` · summary ${num(st.sequences)}` : '')),
      tile('Protein groups', c.protein_groups, `${num(c.proteins_all)} proteins · ${num(c.proteins_top)} top` + (st.protein_groups !== undefined ? ` · summary ${num(st.protein_groups)}` : '')),
      tile('Identified scans', c.identified_scans, st.ms2_scans !== undefined ? `of ${num(st.ms2_scans)} MS2 scans` : null),
      tile('Spectra embedded', (D.spectra_info || {}).embedded, `${(((D.spectra_info || {}).bytes || 0) / 1e6).toFixed(1)} MB · top ${(D.spectra_info || {}).top_n} peaks`));
    sec.appendChild(h('div', { class: 'card' }, h('h2', null, 'Identifications'), tiles,
      h('div', { class: 'muted small', style: 'margin-top:8px' }, 'Counts are computed from the exported CSV files; "summary" values come from db.summary.html when present. Peptides and sequences are counted I/L-collapsed, as in the PEAKS summary.')));

    // filters + parameters
    const kv = rows => h('table', { class: 'kv' }, rows.filter(r => r[1] !== undefined && r[1] !== null && r[1] !== '').map(r => h('tr', null, h('td', null, r[0]), h('td', null, String(r[1])))));
    const unit = (v, u) => (v === undefined || v === null ? undefined : `${v} ${u || ''}`.trim());
    const filt = [
      ['Peptide -10lgP ≥', p.peptide_min_10lgp], ['Protein -10lgP ≥', p.protein_min_10lgp],
      ['Protein unique peptides ≥', p.protein_min_unique_peptides], ['De novo ALC (%) ≥', p.denovo_min_alc_pct],
      ['FDR PSM', unit(p.fdr_psm_pct, '%')], ['FDR peptide sequences', unit(p.fdr_peptide_pct, '%')],
      ['FDR protein groups', unit(p.fdr_protein_group_pct, '%')], ['De novo only spectra', p.denovo_only_spectra]];
    const ps = T.psms, pepT = T.peptides, prT = T.proteins;
    const minOf = (t, name) => { const j = col(t, name); if (j < 0) return undefined; let m = Infinity; for (const r of t.rows) if (typeof r[j] === 'number' && r[j] < m) m = r[j]; return m === Infinity ? undefined : m.toFixed(2); };
    const observed = [['Lowest PSM -10lgP', minOf(ps, '-10LgP')], ['Lowest peptide -10lgP', minOf(pepT, '-10LgP')], ['Lowest protein -10lgP', minOf(prT, '-10LgP')]];
    const fcard = h('div', { class: 'card' }, h('h2', null, 'Filters and FDR'),
      filt.some(r => r[1] !== undefined) ? kv(filt) : h('p', { class: 'muted' }, 'Not in the export (no db.summary.html).'),
      h('h3', null, 'Observed in the tables'), kv(observed));
    const raw = D.raw_params || {};
    let prows = Object.keys(raw).map(k => [k, Array.isArray(raw[k]) ? raw[k].join('; ') : raw[k]]);
    if (!prows.length) {
      const skip = new Set(['statistics', 'missed_cleavages', 'peptide_min_10lgp', 'protein_min_10lgp', 'protein_min_unique_peptides', 'denovo_min_alc_pct', 'fdr_psm_pct', 'fdr_peptide_pct', 'fdr_protein_group_pct', 'denovo_only_spectra']);
      prows = Object.keys(p).filter(k => !skip.has(k)).map(k => [k, typeof p[k] === 'object' ? JSON.stringify(p[k]) : p[k]]);
    }
    const ft = D.frag_tol || {};
    prows.push(['Fragment tolerance used by this viewer', `${ft.value} ${ft.unit} (${ft.source})`]);
    const pcard = h('div', { class: 'card' }, h('h2', null, 'Search parameters'), kv(prows),
      h('h3', null, 'Inputs'), kv([['Export', meta.export_dir], ['FASTA', meta.fasta || 'none'], ['mzML', (meta.mzml || []).join(', ') || 'none'], ['Reference', meta.ref]]));
    sec.appendChild(h('div', { class: 'grid2', style: 'margin-bottom:14px' }, fcard, pcard));

    // plots
    const colVals = (t, name) => { const j = col(t, name); return j < 0 ? [] : t.rows.map(r => r[j]); };
    const thr = typeof p.peptide_min_10lgp === 'number' ? [{ x: p.peptide_min_10lgp, label: `peptide threshold ${p.peptide_min_10lgp}` }] : [];
    const plots = h('div', { class: 'grid2', style: 'margin-bottom:14px' },
      h('div', { class: 'card', html: '<h2>PSM score distribution</h2>' + histogramSVG(colVals(ps, '-10LgP'), { xlabel: 'PSM -10lgP', bins: 50, lines: thr }) }),
      h('div', { class: 'card', html: '<h2>Peptide score distribution</h2>' + histogramSVG(colVals(pepT, '-10LgP'), { xlabel: 'peptide -10lgP', bins: 50, lines: thr }) }),
      h('div', { class: 'card', html: '<h2>Precursor mass error</h2>' + histogramSVG(colVals(ps, 'ppm'), { xlabel: 'precursor error (ppm)', bins: 60, clip: 0.002 }) }),
      h('div', { class: 'card', html: '<h2>Precursor mass error vs m/z</h2>' + scatterSVG(colVals(ps, 'm/z'), colVals(ps, 'ppm'), { xlabel: 'precursor m/z', ylabel: 'error (ppm)', clip: 0.002 }) }));
    sec.appendChild(plots);

    // missed cleavages + PTM profile
    const mc = D.missed_cleavages || {};
    const keys = ['0', '1', '2', '3', '4+'];
    const mrows = [['computed from peptides', mc.computed]];
    if (mc.summary) for (const s in mc.summary) mrows.push([`db.summary.html (${s})`, mc.summary[s]]);
    const mtab = h('table', { class: 'simple' }, h('thead', null, h('tr', null, h('th', null, 'source'), keys.map(k => h('th', null, k)), h('th', null, 'total'))),
      h('tbody', null, mrows.filter(r => r[1]).map(([n, m]) => h('tr', null, h('td', null, n), keys.map(k => h('td', null, num(m[k]))), h('td', null, num(keys.reduce((a, k) => a + (m[k] || 0), 0)))))));
    const enz = mc.enzyme || {};
    const mcard = h('div', { class: 'card' }, h('h2', null, 'Missed cleavages'), mtab,
      h('div', { class: 'muted small', style: 'margin-top:6px' }, `Internal sites after ${enz.cleave_after || '–'}${enz.not_before ? ', not before ' + enz.not_before : ''}; one count per peptide row.`));
    const pt = D.ptm_profile || { columns: [], rows: [] };
    const ptab = h('table', { class: 'simple' }, h('thead', null, h('tr', null, pt.columns.map(c => h('th', null, c)))),
      h('tbody', null, pt.rows.map(r => h('tr', null, r.map((v, j) => h('td', null, fmt(v, pt.kinds[j])))))));
    const ptcard = h('div', { class: 'card' }, h('h2', null, 'PTM profile'), pt.rows.length ? ptab : h('p', { class: 'muted' }, 'No modifications.'),
      h('div', { class: 'muted small', style: 'margin-top:6px' }, 'Source: ' + (pt.source || '')));
    sec.appendChild(h('div', { class: 'grid2', style: 'margin-bottom:14px' }, mcard, ptcard));

    if (D.concordance) sec.appendChild(concordancePanel(D.concordance));
    if ((D.notes || []).length) sec.appendChild(h('div', { class: 'card' }, h('h2', null, 'Notes'), h('ul', { class: 'notes' }, D.notes.map(n => h('li', null, n)))));
  };

  function concordancePanel(C) {
    const m = C.metrics || {}, L = m.levels || {};
    const card = h('div', { class: 'card', id: 'concordance' }, h('h2', null, 'Concordance with the reference'));
    card.appendChild(h('div', { class: 'muted small' }, 'Reference: ', h('span', { class: 'mono' }, m.ref || ''), h('br'), 'This export: ', h('span', { class: 'mono' }, m.test || '')));
    const t3 = h('div', { class: 'tiles', style: 'margin:10px 0' });
    for (const [k, lab] of [['psm', 'PSM'], ['peptide', 'Peptide'], ['protein_group', 'Protein group']]) {
      const x = L[k]; if (!x) continue;
      t3.appendChild(h('div', { class: 'tile' }, h('div', { class: 'k' }, lab), h('div', { class: 'v' }, `${pct(x.recall)}`), h('div', { class: 'x' }, `recall · precision ${pct(x.precision)}`), h('div', { class: 'x' }, `${num(x.shared)} shared · ${num(x.ref_only)} ref-only · ${num(x.test_only)} test-only`)));
    }
    card.appendChild(t3);
    const order = ['psm', 'scan', 'peptide', 'peptide_sequence', 'peptide_IL_collapsed', 'protein_group', 'protein_top_accession', 'protein_accession', 'feature'];
    const main3 = new Set(['psm', 'peptide', 'protein_group']);
    card.appendChild(h('table', { class: 'simple' },
      h('thead', null, h('tr', null, ['level', 'ref', 'test', 'shared', 'ref-only', 'test-only', 'recall', 'precision'].map(x => h('th', null, x)))),
      h('tbody', null, order.filter(k => L[k]).map(k => h('tr', { class: main3.has(k) ? 'hl' : '' }, h('td', null, k), [L[k].n_ref, L[k].n_test, L[k].shared, L[k].ref_only, L[k].test_only].map(v => h('td', null, num(v))), h('td', null, pct(L[k].recall)), h('td', null, pct(L[k].precision)))))));
    card.appendChild(h('div', { class: 'muted small', style: 'margin-top:4px' }, 'PSM key: source file + raw scan + exact modified sequence + charge (I ≠ L). Protein groups match when they share an accession. Levels in bold are the headline levels.'));
    const rb = m.psm_ref_only_breakdown || {}, tb = m.psm_test_only_breakdown || {};
    if (Object.keys(rb).length) {
      card.appendChild(h('h3', null, 'Why PSMs differ'));
      card.appendChild(h('table', { class: 'simple' }, h('thead', null, h('tr', null, h('th', null, 'category'), h('th', null, 'ref-only'), h('th', null, 'test-only'))),
        h('tbody', null, Object.keys(rb).map(k => h('tr', null, h('td', null, k.replace(/_/g, ' ')), h('td', null, num(rb[k])), h('td', null, num(tb[k])))))));
    }
    // lists
    const lists = [['psm_ref_only', 'PSMs: reference only'], ['psm_test_only', 'PSMs: this export only'], ['peptide_ref_only', 'Peptides: reference only'],
      ['peptide_test_only', 'Peptides: this export only'], ['protein_group_ref_only', 'Protein groups: reference only'], ['protein_group_test_only', 'Protein groups: this export only'],
      ['protein_top_ref_only', 'Top accessions: reference only'], ['protein_top_test_only', 'Top accessions: this export only']].filter(x => C.tables[x[0]]);
    if (lists.length) {
      card.appendChild(h('h3', null, 'Reference-only and test-only items'));
      const host = h('div');
      const sel = h('select', { onchange: () => show(sel.value) }, lists.map(([k, lab]) => h('option', { value: k }, `${lab} (${num(C.tables[k].rows.length)})`)));
      card.appendChild(h('div', { class: 'ctl' }, h('label', null, 'List ', sel), h('span', { class: 'muted small' }, 'Click a PSM row to view its spectrum (when embedded) annotated with that peptide.')));
      card.appendChild(host);
      const show = k => {
        const t = C.tables[k];
        new DataTable(host, { table: t, name: k, scoreCol: 'score', pageSize: 25, maxHeight: '420px', sort: col(t, 'score') >= 0 ? ['score', true] : null,
          onRowClick: i => {
            const r = t.rows[i], g = n => r[col(t, n)];
            if (k.startsWith('psm_')) openSpectrumModal({ peptide: g('peptide'), z: g('z'), scan: g('scan'), source: g('source'), score: g('score'), origin: k.startsWith('psm_ref') ? 'reference' : 'this export', note: g('category') ? 'category: ' + g('category') + (g('other_side_at_scan') ? ' · other side at this scan: ' + g('other_side_at_scan') : '') : '' });
            else if (k === 'peptide_test_only') openPeptideModal(g('peptide'));
            else if (k.startsWith('protein') && k.endsWith('test_only')) { const a = (g('top_accessions') || g('accession') || '').split(';')[0]; if (IDX.protByAcc.has(a)) { showTab('proteins', true); PROT.select(a); } }
          } });
      };
      show(lists[0][0]);
    }
    if ((m.warnings || []).length) card.appendChild(h('ul', { class: 'notes' }, m.warnings.map(w => h('li', null, w))));
    return card;
  }

  // ================================================================== Proteins
  const PROT = {};
  builders.proteins = function (sec) {
    const t = T.proteins;
    if (!t) { sec.appendChild(h('div', { class: 'card' }, 'This export has no db.proteins.csv.')); return; }
    const itop = col(t, 'Top');
    let topOnly = false;
    const card = h('div', { class: 'card' }, h('h2', null, 'Protein groups'));
    const host = h('div'); card.appendChild(host);
    const detail = h('div', { class: 'detail' });
    sec.append(card, detail);
    const hidden = t.columns.filter(c => /^Coverage\(%\) /.test(c) || c === 'Average Mass' || c === 'PTM');
    const dt = new DataTable(host, { table: t, name: 'proteins', sort: ['-10LgP', true], hidden,
      filter: i => !topOnly || itop < 0 || t.rows[i][itop] === true,
      controls: bar => { if (itop >= 0) bar.appendChild(h('label', null, h('input', { type: 'checkbox', onchange: e => { topOnly = e.target.checked; dt.page = 0; dt.refresh(); } }), ' top proteins only')); },
      rowClass: i => (itop >= 0 && t.rows[i][itop] !== true ? 'dim' : ''),
      onRowClick: i => showProtein(i) });
    function showProtein(i) {
      const r = t.rows[i], acc = r[col(t, 'Accession')];
      setHash({ tab: 'proteins', protein: acc });
      detail.innerHTML = '';
      const g = n => { const j = col(t, n); return j >= 0 ? r[j] : null; };
      const pcard = h('div', { class: 'card' },
        h('h2', null, acc, g('Gene') ? ` · ${g('Gene')}` : ''),
        h('div', { class: 'muted' }, g('Description') || ''),
        h('div', { class: 'ctl' }, h('span', null, `group ${g('Protein Group')}`), h('span', null, `-10lgP ${fmt(g('-10LgP'), t.kinds[col(t, '-10LgP')])}`),
          h('span', null, `coverage ${fmt(g('Coverage(%)'), 'fixed2')}%`), h('span', null, `${g('#Peptides')} peptides (${g('#Unique')} unique)`),
          g('PTM') ? h('span', null, 'PTM: ' + g('PTM')) : null));
      detail.appendChild(pcard);
      // coverage
      const seq = (D.sequences || {})[acc];
      const pp = T.protein_peptides, rows = IDX.ppByAcc.get(acc) || [];
      const peps = [];
      if (pp) {
        const ip = col(pp, 'Peptide'), is = col(pp, 'Start'), iu = col(pp, 'Unique'), isc = col(pp, '-10LgP');
        for (const k of rows) {
          const rr = pp.rows[k], pep = rr[ip];
          let start = is >= 0 && typeof rr[is] === 'number' ? rr[is] - 1 : -1;
          if (start < 0 && seq) start = seq.replace(/I/g, 'L').indexOf(stripped(pep).replace(/I/g, 'L'));
          peps.push({ pep, start, len: stripped(pep).length, unique: iu >= 0 ? rr[iu] === true : true, score: isc >= 0 ? rr[isc] : null });
        }
      }
      pcard.appendChild(h('h3', null, 'Sequence coverage'));
      if (seq) {
        pcard.appendChild(h('div', { class: 'legend' }, h('span', null, h('i', { style: 'background:var(--cov)' }), 'covered residue'), h('span', null, h('i', { style: 'background:var(--mod)' }), 'modified residue'),
          h('span', null, h('i', { style: 'background:var(--bar-unique)' }), 'unique peptide'), h('span', null, h('i', { style: 'background:var(--bar-shared)' }), 'shared peptide'), h('span', null, 'click a bar for its spectrum')));
        pcard.appendChild(h('div', { class: 'cov-wrap', html: coverageSVG(seq, peps) }));
        pcard.querySelector('.cov-wrap').addEventListener('click', e => { const b = e.target.closest('[data-pep]'); if (b) openPeptideModal(b.getAttribute('data-pep')); });
      } else pcard.appendChild(h('p', { class: 'muted' }, 'Protein sequence not available (run with --fasta).'));
      // members of the group
      const grp = IDX.protByGroup.get(g('Protein Group')) || [];
      if (grp.length > 1) {
        pcard.appendChild(h('h3', null, `Proteins in group ${g('Protein Group')} (${grp.length})`));
        pcard.appendChild(h('table', { class: 'simple' }, h('thead', null, h('tr', null, ['Accession', 'Top', 'Gene', '-10lgP', 'Coverage(%)', '#Peptides', '#Unique'].map(x => h('th', null, x)))),
          h('tbody', null, grp.map(k => { const q = t.rows[k], G = n => { const j = col(t, n); return j >= 0 ? q[j] : ''; };
            return h('tr', { style: 'cursor:pointer', onclick: () => { dt.selectRow(k); showProtein(k); } }, h('td', null, G('Accession')), h('td', null, String(G('Top'))), h('td', null, G('Gene')), h('td', null, fmt(G('-10LgP'), 'fixed2')), h('td', null, fmt(G('Coverage(%)'), 'fixed2')), h('td', null, num(G('#Peptides'))), h('td', null, num(G('#Unique')))); }))));
      }
      // supporting peptides
      const scard = h('div', { class: 'card' }, h('h2', null, `Supporting peptides (${rows.length})`), h('div', { class: 'muted small' }, 'Click a row to view the annotated spectrum of the peptide\'s best PSM.'));
      const shost = h('div'); scard.appendChild(shost);
      detail.appendChild(scard);
      if (pp && rows.length) new DataTable(shost, { table: pp, rowsIdx: rows, name: 'protein-peptides-' + acc, pageSize: 25, sort: ['Start', false], maxHeight: '420px',
        hidden: ['Protein Group', 'Accession', 'Gene', 'Top', 'Found By'], onRowClick: k => openPeptideModal(pp.rows[k][col(pp, 'Peptide')]) });
      else shost.appendChild(h('p', { class: 'muted' }, 'No protein-peptide rows for this protein.'));
    }
    PROT.select = acc => { const i = IDX.protByAcc.get(acc); if (i === undefined) return; dt.setText(''); dt.selectRow(i); showProtein(i); };
    PROT.dt = dt;
  };

  function coverageSVG(seq, peps) {
    const PER = 60, CW = 10.2, LX = 56, n = seq.length;
    const cov = new Uint16Array(n), mod = new Array(n).fill(null);
    for (const p of peps) {
      if (p.start < 0) continue;
      for (let k = p.start; k < Math.min(n, p.start + p.len); k++) cov[k]++;
      const info = D.pepinfo[p.pep];
      if (info) for (const [pos, d] of info[1]) {
        const k = p.start + Math.min(Math.max(pos, 1), p.len) - 1;
        if (k >= 0 && k < n) (mod[k] = mod[k] || new Set()).add((d >= 0 ? '+' : '') + d.toFixed(2));
      }
    }
    let covered = 0; for (let k = 0; k < n; k++) if (cov[k]) covered++;
    let y = 4, s = '';
    const nl = Math.ceil(n / PER);
    for (let L = 0; L < nl; L++) {
      const a = L * PER, b = Math.min(n, a + PER);
      y += 16;
      s += `<text class="pos" x="${LX - 8}" y="${y}" text-anchor="end">${a + 1}</text>`;
      let txt = '';
      for (let k = a; k < b; k++) {
        const x = LX + (k - a) * CW + Math.floor((k - a) / 10) * 5;
        if (cov[k]) s += `<rect class="cb" x="${x - 0.5}" y="${y - 12}" width="${CW + 0.5}" height="16"/>`;
        const cls = mod[k] ? 'res m' : cov[k] ? 'res c' : 'res';
        txt += `<text class="${cls}" x="${x + CW / 2}" y="${y}" text-anchor="middle">${seq[k]}${mod[k] ? `<title>${k + 1}${seq[k]}: ${[...mod[k]].join(', ')}</title>` : ''}</text>`;
      }
      s += txt;
      // peptide bars overlapping this line, greedy packing into rows
      const lanes = [];
      const segs = peps.filter(p => p.start >= 0 && p.start < b && p.start + p.len > a).sort((p, q) => p.start - q.start || q.len - p.len);
      y += 6;
      for (const p of segs) {
        const s0 = Math.max(a, p.start), s1 = Math.min(b, p.start + p.len);
        let lane = lanes.findIndex(end => end <= s0);
        if (lane < 0) { lane = lanes.length; lanes.push(0); }
        lanes[lane] = s1;
        const X = k => LX + (k - a) * CW + Math.floor((k - a) / 10) * 5;
        const yy = y + lane * 7;
        s += `<rect class="pb ${p.unique ? 'u' : 's'}" data-pep="${esc(p.pep)}" x="${X(s0)}" y="${yy}" width="${X(s1 - 1) + CW - X(s0)}" height="5" rx="1.5"><title>${esc(p.pep)} (${p.start + 1}-${p.start + p.len})${p.score !== null ? ' -10lgP ' + p.score : ''}</title></rect>`;
        const info = D.pepinfo[p.pep];
        if (info) for (const [pos] of info[1]) { const k = p.start + Math.min(Math.max(pos, 1), p.len) - 1; if (k >= s0 && k < s1) s += `<circle class="mk" cx="${X(k) + CW / 2}" cy="${yy + 2.5}" r="2.6"/>`; }
      }
      y += lanes.length * 7 + 4;
    }
    const W = LX + PER * CW + 6 * 5 + 10;
    const unplaced = peps.filter(p => p.start < 0).length;
    return `<div class="muted small">${n} residues · ${covered} covered (${(100 * covered / n).toFixed(1)}%) · ${peps.length} peptide rows${unplaced ? ` · ${unplaced} not located` : ''}</div>` +
      `<svg class="cov" viewBox="0 0 ${W} ${y + 4}" width="${W}" height="${y + 4}">${s}</svg>`;
  }

  // ================================================================== Peptides
  builders.peptides = function (sec) {
    const t = T.peptides;
    if (!t) { sec.appendChild(h('div', { class: 'card' }, 'This export has no db.peptides.csv.')); return; }
    const card = h('div', { class: 'card' }, h('h2', null, 'Peptides'));
    const host = h('div'); card.appendChild(host);
    const detail = h('div', { class: 'detail' });
    sec.append(card, detail);
    const dt = new DataTable(host, { table: t, name: 'peptides', sort: ['-10LgP', true], onRowClick: i => showPeptide(i) });
    function showPeptide(i) {
      const pep = t.rows[i][col(t, 'Peptide')];
      setHash({ tab: 'peptides', peptide: pep });
      detail.innerHTML = '';
      const sc = h('div', { class: 'card' }); detail.appendChild(sc);
      const bi = bestPsm(pep);
      if (bi >= 0) new SpectrumView(sc, psmObj(bi)); else sc.appendChild(h('p', { class: 'muted' }, 'No PSM rows for this peptide.'));
      const list = IDX.psmByPep.get(pep) || [];
      const pc = h('div', { class: 'card' }, h('h2', null, `PSMs of ${pep} (${list.length})`), h('div', { class: 'muted small' }, 'Click a PSM to show its spectrum above.'));
      const ph = h('div'); pc.appendChild(ph); detail.appendChild(pc);
      if (list.length) new DataTable(ph, { table: T.psms, rowsIdx: list, name: 'psms-' + pep, pageSize: 25, sort: ['-10LgP', true], maxHeight: '360px',
        onRowClick: k => { sc.innerHTML = ''; new SpectrumView(sc, psmObj(k)); sc.scrollIntoView({ block: 'start' }); } });
      const acc = t.rows[i][col(t, 'Accession')];
      if (acc) {
        const accs = String(acc).split(/[:;]/).filter(Boolean);
        pc.appendChild(h('h3', null, 'Proteins'));
        pc.appendChild(h('div', null, accs.map(a => IDX.protByAcc.has(a) ? h('a', { href: '#', style: 'margin-right:12px', onclick: e => { e.preventDefault(); showTab('proteins', true); PROT.select(a); } }, a) : h('span', { style: 'margin-right:12px' }, a))));
      }
    }
    builders.peptides.select = pep => { const i = IDX.pepRow.get(pep); if (i === undefined) return; dt.selectRow(i); showPeptide(i); };
  };

  // ================================================================== PSMs
  builders.psms = function (sec) {
    const t = T.psms;
    if (!t) { sec.appendChild(h('div', { class: 'card' }, 'This export has no db.psms.csv.')); return; }
    const card = h('div', { class: 'card' }, h('h2', null, 'Peptide-spectrum matches'));
    const host = h('div'); card.appendChild(host);
    const detail = h('div', { class: 'detail' });
    sec.append(card, detail);
    let withSpec = false;
    const isc = col(t, 'Scan'), isr = col(t, 'Source File');
    const dt = new DataTable(host, { table: t, name: 'psms', sort: ['-10LgP', true],
      hidden: ['Positional Confidence', 'tag(>=0.0%)'],
      filter: i => !withSpec || hasSpectrum(t.rows[i][isr], t.rows[i][isc]),
      controls: bar => bar.appendChild(h('label', { title: 'Only PSMs whose MS2 peaks are embedded in this report' }, h('input', { type: 'checkbox', onchange: e => { withSpec = e.target.checked; dt.page = 0; dt.refresh(); } }), ' with spectrum')),
      onRowClick: i => showPsm(i) });
    function showPsm(i) {
      setHash({ tab: 'psms', psm: i });
      detail.innerHTML = '';
      const sc = h('div', { class: 'card' }); detail.appendChild(sc);
      new SpectrumView(sc, psmObj(i));
    }
    builders.psms.select = i => { if (i < 0 || i >= t.rows.length) return; dt.selectRow(i); showPsm(i); };
  };

  // ================================================================== spectrum
  function theoreticalIons(info, z, opt) {
    const s = info[0], n = s.length, md = new Float64Array(n + 2);
    for (const [p, d] of info[1]) md[Math.min(Math.max(p, 0), n + 1)] += d;
    const r = new Float64Array(n);
    for (let k = 0; k < n; k++) r[k] = (AA[s[k]] || 0) + md[k + 1];
    const pre = new Float64Array(n + 1); pre[0] = md[0];
    for (let k = 0; k < n; k++) pre[k + 1] = pre[k] + r[k];
    const M = pre[n] + md[n + 1] + H2O;
    const charges = [1]; if (opt.two) charges.push(2);
    const ions = [];
    const has = (a, b, set) => { for (let k = a; k < b; k++) if (set.indexOf(s[k]) >= 0) return true; return false; };
    for (let i = 1; i < n; i++) {
      const bN = pre[i], yN = M - pre[n - i];
      const segs = [['b', i, bN, 0, i], ['y', i, yN, n - i, n]];
      for (const [ser, k, neutral, a, b] of segs) {
        if ((ser === 'b' && !opt.b) || (ser === 'y' && !opt.y)) continue;
        const losses = [['', 0]];
        if (opt.loss && has(a, b, 'STED')) losses.push(['-H2O', H2O]);
        if (opt.loss && has(a, b, 'RKNQ')) losses.push(['-NH3', NH3]);
        for (const c of charges) for (const [ln, lm] of losses) ions.push({ ser, k, c, loss: ln, mz: (neutral - lm + c * PROTON) / c });
      }
    }
    return { ions, M, mzPrec: z ? (M + z * PROTON) / z : null };
  }
  function matchIons(ions, sp, tol, unit) {
    const mz = sp.mz, I = sp.i;
    for (const ion of ions) {
      const w = unit === 'ppm' ? ion.mz * tol * 1e-6 : tol;
      let lo = 0, hi = mz.length;
      while (lo < hi) { const m = (lo + hi) >> 1; if (mz[m] < ion.mz - w) lo = m + 1; else hi = m; }
      let best = -1;
      for (let k = lo; k < mz.length && mz[k] <= ion.mz + w; k++) if (best < 0 || I[k] > I[best]) best = k;
      ion.peak = best;
      if (best >= 0) { ion.err = mz[best] - ion.mz; ion.ppm = ion.err / ion.mz * 1e6; }
    }
  }
  function ionLabel(ion, plain) {
    const sup = plain ? (ion.c > 1 ? ion.c + '+' : '') : (ion.c > 1 ? (ion.c === 2 ? '²⁺' : ion.c + '+') : '');
    const loss = plain ? ion.loss : (ion.loss === '-H2O' ? '°' : ion.loss === '-NH3' ? '*' : '');
    return ion.ser + ion.k + (plain ? (sup ? ' ' + sup : '') + ion.loss : sup + loss);
  }

  class SpectrumView {
    constructor(host, psm) {
      this.host = host; this.psm = psm;
      this.info = D.pepinfo[psm.peptide] || null;
      this.sp = psm.scan !== null && psm.scan !== undefined ? spectra()[specKey(psm.source, psm.scan)] || null : null;
      const z = +psm.z || (this.sp && this.sp.pz) || 2;
      const ft = D.frag_tol || { value: 0.5, unit: 'Da' };
      this.opt = { b: true, y: true, two: z >= 3, loss: false, tol: ft.value, unit: ft.unit };
      this.z = z; this.xr = null;
      this.build();
    }
    build() {
      const p = this.psm, self = this;
      this.host.innerHTML = '';
      const head = h('div', { class: 'spec-head' },
        h('span', { class: 'pep', html: pepHtml(p.peptide || '') }),
        h('span', null, `z = ${this.z}+`), p.mz !== null && p.mz !== undefined ? h('span', null, `m/z ${fmt(p.mz, 'fixed5')}`) : null,
        h('span', null, `scan ${p.scan}`), h('span', { class: 'muted' }, p.source || ''),
        p.score !== null && p.score !== undefined ? h('span', null, `-10lgP ${fmt(p.score, 'fixed2')}`) : null,
        p.ppm !== null && p.ppm !== undefined ? h('span', null, `${fmt(p.ppm, 'fixed1')} ppm`) : null,
        p.rt !== null && p.rt !== undefined ? h('span', null, `RT ${fmt(p.rt, 'fixed2')} min`) : null,
        p.origin ? h('span', { class: 'pill' }, 'peptide from ' + p.origin) : null);
      this.host.appendChild(h('h2', null, 'Annotated MS2 spectrum'));
      this.host.appendChild(head);
      if (p.note) this.host.appendChild(h('div', { class: 'muted small' }, p.note));
      if (!this.info) { this.host.appendChild(h('div', { class: 'notice' }, 'Cannot parse this peptide sequence.')); return; }
      const cb = (k, lab, title) => h('label', { title: title || null }, h('input', { type: 'checkbox', checked: this.opt[k], onchange: e => { self.opt[k] = e.target.checked; self.update(); } }), ' ', lab);
      const tolIn = h('input', { type: 'number', step: 'any', min: '0', value: this.opt.tol, onchange: () => { const v = parseFloat(tolIn.value); if (v > 0) { self.opt.tol = v; self.update(); } } });
      const unitSel = h('select', { onchange: () => { self.opt.unit = unitSel.value; self.update(); } }, ['Da', 'ppm'].map(u => h('option', { value: u, selected: u === this.opt.unit }, u)));
      this.host.appendChild(h('div', { class: 'ctl' }, cb('b', 'b ions'), cb('y', 'y ions'), cb('two', '2+ fragments', 'default on for precursor charge ≥ 3'), cb('loss', 'neutral losses (H₂O °, NH₃ *)'),
        h('label', null, 'tolerance ', tolIn, ' ', unitSel), h('button', { class: 'btn', onclick: () => { self.xr = null; self.draw(); } }, 'Reset zoom'),
        h('span', { class: 'legend' }, h('span', null, h('i', { style: 'background:var(--ion-b)' }), 'b'), h('span', null, h('i', { style: 'background:var(--ion-y)' }), 'y'), h('span', null, h('i', { style: 'background:var(--ion-loss)' }), 'loss'), h('span', null, h('i', { style: 'background:var(--peak)' }), 'unassigned'))));
      if (!this.sp) this.host.appendChild(h('div', { class: 'notice' }, 'The MS2 peaks of this scan are not embedded in the report (only the best PSM of each peptide is embedded, and only when --mzml is given). Theoretical fragments are listed below.'));
      this.ladder = h('div'); this.plot = h('div', { class: 'spec-plot' }); this.status = h('div', { class: 'spec-status' }); this.errp = h('div');
      this.stats = h('div', { class: 'muted small' }); this.table = h('div');
      this.host.append(this.ladder, this.plot, this.status, this.errp, this.stats, h('h3', null, 'Fragment ions'), this.table);
      this.update();
    }
    update() {
      const th = theoreticalIons(this.info, this.z, this.opt);
      this.ions = th.ions; this.M = th.M;
      if (this.sp) matchIons(this.ions, this.sp, this.opt.tol, this.opt.unit); else for (const i of this.ions) i.peak = -1;
      this.byPeak = new Map();
      const pri = i => (i.loss ? 10 : 0) + (i.c - 1) * 3 + (i.ser === 'y' ? 0 : 1);
      for (const ion of this.ions) if (ion.peak >= 0) {
        const cur = this.byPeak.get(ion.peak);
        if (!cur) this.byPeak.set(ion.peak, [ion]); else { cur.push(ion); cur.sort((a, b) => pri(a) - pri(b)); }
      }
      this.draw(); this.drawLadder(); this.drawTable(); this.drawStats();
    }
    drawStats() {
      const main = this.ions.filter(i => !i.loss && i.c === 1), mm = main.filter(i => i.peak >= 0).length;
      let tot = 0, ex = 0;
      if (this.sp) { for (const v of this.sp.i) tot += v; for (const k of this.byPeak.keys()) ex += this.sp.i[k]; }
      const sp = this.sp;
      this.stats.textContent = `Theoretical [M] ${this.M.toFixed(4)} Da · 1+ b/y matched ${mm}/${main.length} · all ions matched ${this.ions.filter(i => i.peak >= 0).length}/${this.ions.length}` +
        (sp ? ` · ${(100 * ex / (tot || 1)).toFixed(1)}% of the embedded intensity explained · ${sp.mz.length} peaks shown (top-N) · base peak ${sci(sp.bp)}${sp.pmz ? ' · scan precursor ' + sp.pmz + (sp.pz ? ' (' + sp.pz + '+)' : '') : ''}${sp.rt !== undefined ? ' · RT ' + sp.rt.toFixed(2) + ' min' : ''}` : '');
    }
    draw() {
      const sp = this.sp;
      if (!sp) { this.plot.innerHTML = ''; this.errp.innerHTML = ''; return; }
      const W = 1000, H = 360, m = { l: 52, r: 14, t: 34, b: 36 };
      const mzs = sp.mz;
      let lo, hi;
      if (this.xr) [lo, hi] = this.xr;
      else { const mn = mzs.length ? mzs[0] : 0, mx = mzs.length ? mzs[mzs.length - 1] : 1000; lo = Math.max(0, Math.floor((mn - 25) / 50) * 50); hi = Math.ceil((mx + 25) / 50) * 50; }
      let ymax = 0;
      for (let k = 0; k < mzs.length; k++) if (mzs[k] >= lo && mzs[k] <= hi && sp.i[k] > ymax) ymax = sp.i[k];
      ymax = (ymax || 100) * 1.12;
      const ax = axes(W, H, m, lo, hi, 0, ymax, 'm/z', 'relative intensity (%)');
      this.X = ax.X; this.lo = lo; this.hi = hi; this.W = W; this.m = m;
      let s = `<svg class="chart" viewBox="0 0 ${W} ${H}">` + ax.s;
      const lab = [];
      let peaks = '';
      for (let k = 0; k < mzs.length; k++) {
        if (mzs[k] < lo || mzs[k] > hi) continue;
        const x = ax.X(mzs[k]).toFixed(1), y = ax.Y(sp.i[k]).toFixed(1);
        const ann = this.byPeak.get(k);
        const cls = ann ? (ann[0].loss ? 'l' : ann[0].ser) : '';
        peaks += `<line class="pk ${cls}" data-k="${k}" x1="${x}" x2="${x}" y1="${ax.Y(0)}" y2="${y}"/>`;
        if (ann) lab.push([k, +x, +y, ann[0]]);
      }
      // labels: most intense first, skip if too close to one already placed
      lab.sort((a, b) => sp.i[b[0]] - sp.i[a[0]]);
      const placed = [];
      let labels = '';
      for (const [k, x, y, ion] of lab) {
        if (placed.some(px => Math.abs(px - x) < 11)) continue;
        placed.push(x);
        const cls = ion.loss ? 'l' : ion.ser;
        labels += `<text class="il ${cls}" transform="translate(${x + 3.5} ${Math.max(m.t + 2, y - 4)}) rotate(-90)">${esc(ionLabel(ion))}</text>`;
      }
      if (this.psm.mz && this.psm.mz >= lo && this.psm.mz <= hi) s += `<text class="lbl" x="${ax.X(this.psm.mz)}" y="${m.t - 20}" text-anchor="middle">▼ precursor</text>`;
      s += peaks + labels + `<rect class="zoom" x="0" y="${m.t}" width="0" height="${H - m.t - m.b}"/></svg>`;
      this.plot.innerHTML = s;
      this.bindPlot(this.plot.querySelector('svg'), ax, H);
      this.drawErr();
    }
    bindPlot(svg, ax, H) {
      const self = this, sp = this.sp, zr = svg.querySelector('.zoom');
      const toX = e => { const b = svg.getBoundingClientRect(); return (e.clientX - b.left) * self.W / b.width; };
      const toMz = x => self.lo + (x - self.m.l) / (self.W - self.m.l - self.m.r) * (self.hi - self.lo);
      let x0 = null;
      svg.addEventListener('mousedown', e => { x0 = toX(e); zr.setAttribute('x', x0); zr.setAttribute('width', 0); e.preventDefault(); });
      svg.addEventListener('mousemove', e => {
        const x = toX(e);
        if (x0 !== null) { zr.setAttribute('x', Math.min(x0, x)); zr.setAttribute('width', Math.abs(x - x0)); return; }
        // nearest peak within 6 px
        let best = -1, bd = 6;
        for (let k = 0; k < sp.mz.length; k++) { if (sp.mz[k] < self.lo || sp.mz[k] > self.hi) continue; const d = Math.abs(ax.X(sp.mz[k]) - x); if (d < bd) { bd = d; best = k; } }
        for (const el of svg.querySelectorAll('.pk.hot')) el.classList.remove('hot');
        if (best >= 0) {
          const el = svg.querySelector(`.pk[data-k="${best}"]`); if (el) el.classList.add('hot');
          const ann = self.byPeak.get(best);
          self.status.textContent = `m/z ${sp.mz[best].toFixed(4)}  ${sp.i[best].toFixed(1)}%` + (ann ? '  ' + ann.map(i => `${ionLabel(i, true)} (${i.mz.toFixed(4)}, ${i.err >= 0 ? '+' : ''}${i.err.toFixed(4)} Da / ${i.ppm.toFixed(1)} ppm)`).join(', ') : '  unassigned');
        } else self.status.textContent = 'drag to zoom · double-click to reset';
      });
      const end = e => {
        if (x0 === null) return;
        const x1 = toX(e), a = Math.min(x0, x1), b = Math.max(x0, x1); x0 = null;
        zr.setAttribute('width', 0);
        if (b - a > 6) { self.xr = [toMz(a), toMz(b)]; self.draw(); }
      };
      svg.addEventListener('mouseup', end);
      svg.addEventListener('mouseleave', e => { if (x0 !== null) end(e); });
      svg.addEventListener('dblclick', () => { self.xr = null; self.draw(); });
    }
    drawErr() {
      const W = 1000, H = 130, m = { l: 52, r: 14, t: 10, b: 30 };
      const ppm = this.opt.unit === 'ppm', tol = this.opt.tol;
      const ax = axes(W, H, m, this.lo, this.hi, -tol, tol, 'm/z', ppm ? 'error (ppm)' : 'error (Da)');
      let s = `<svg class="chart" viewBox="0 0 ${W} ${H}">` + ax.s + `<line x1="${m.l}" x2="${W - m.r}" y1="${ax.Y(0)}" y2="${ax.Y(0)}" stroke="currentColor" opacity=".3"/>`;
      for (const i of this.ions) if (i.peak >= 0 && i.mz >= this.lo && i.mz <= this.hi) {
        const cls = i.loss ? 'l' : i.ser;
        s += `<circle class="ep ${cls}" cx="${ax.X(i.mz).toFixed(1)}" cy="${ax.Y(ppm ? i.ppm : i.err).toFixed(1)}" r="3"><title>${esc(ionLabel(i, true))}: ${i.err.toFixed(4)} Da, ${i.ppm.toFixed(1)} ppm</title></circle>`;
      }
      this.errp.innerHTML = s + '</svg>';
    }
    drawLadder() {
      const s = this.info[0], n = s.length, CW = 26, W = Math.max(360, n * CW + 40), H = 64, x0 = 20;
      const mods = new Map(); for (const [p, d] of this.info[1]) { const q = Math.min(Math.max(p, 1), n); mods.set(q, (mods.get(q) || 0) + d); }
      const got = { b: new Set(), y: new Set() };
      for (const i of this.ions) if (i.peak >= 0) got[i.ser].add(i.k);
      let g = `<svg class="chart ladder" viewBox="0 0 ${W} ${H}" style="max-width:${W}px">`;
      for (let k = 0; k < n; k++) {
        const x = x0 + k * CW + CW / 2, md = mods.get(k + 1);
        g += `<text x="${x}" y="38" text-anchor="middle" class="${md !== undefined ? 'mod' : ''}">${s[k]}${md !== undefined ? `<title>${s[k]}${k + 1} ${md >= 0 ? '+' : ''}${md.toFixed(4)}</title>` : ''}</text>`;
        if (md !== undefined) g += `<text x="${x}" y="54" text-anchor="middle" class="mod" style="font-size:9px">${md >= 0 ? '+' : ''}${md.toFixed(2)}</text>`;
      }
      for (let k = 1; k < n; k++) {
        const xb = x0 + k * CW;
        if (got.b.has(k)) g += `<path class="tb" d="M${xb} 26 V44 H${xb - 8}"/>`;
        if (got.y.has(n - k)) g += `<path class="ty" d="M${xb} 44 V26 H${xb + 8}"/>`;
      }
      g += `<text class="lbl" x="${x0}" y="12">b ions below-left, y ions above-right of each cleavage</text></svg>`;
      this.ladder.innerHTML = g;
    }
    drawTable() {
      const s = this.info[0], n = s.length, o = this.opt;
      const kinds = [];
      const charges = o.two ? [1, 2] : [1];
      const losses = o.loss ? ['', '-H2O', '-NH3'] : [''];
      for (const c of charges) for (const l of losses) kinds.push({ c, l });
      const find = new Map();
      for (const i of this.ions) find.set(`${i.ser}|${i.k}|${i.c}|${i.loss}`, i);
      const hd = (ser, c, l) => `${ser}${c > 1 ? ' ' + c + '+' : ''}${l}`;
      const cell = (ser, k, c, l) => {
        const i = find.get(`${ser}|${k}|${c}|${l}`);
        if (!i) return '<td></td>';
        const cls = i.peak >= 0 ? (l ? 'ml' : ser === 'b' ? 'mb' : 'my') : '';
        const t = i.peak >= 0 ? ` title="observed ${this.sp.mz[i.peak].toFixed(4)} (${this.sp.i[i.peak].toFixed(1)}%), error ${i.err >= 0 ? '+' : ''}${i.err.toFixed(4)} Da / ${i.ppm.toFixed(1)} ppm"` : '';
        return `<td class="${cls}"${t}>${i.mz.toFixed(4)}</td>`;
      };
      const bk = o.b ? kinds : [], yk = o.y ? kinds : [];
      const mods = new Map(); for (const [p, d] of this.info[1]) { const q = Math.min(Math.max(p, 1), n); mods.set(q, (mods.get(q) || 0) + d); }
      let html = '<table class="frag"><thead><tr><th>#</th>' + bk.map(k => `<th>${hd('b', k.c, k.l)}</th>`).join('') + '<th>AA</th>' + yk.map(k => `<th>${hd('y', k.c, k.l)}</th>`).join('') + '<th>#</th><th>err b (Da)</th><th>err y (Da)</th></tr></thead><tbody>';
      for (let k = 1; k <= n; k++) {
        const md = mods.get(k), yk_ = n - k + 1;
        const eb = find.get(`b|${k}|1|`), ey = find.get(`y|${yk_}|1|`);
        html += `<tr><td>${k}</td>` + bk.map(q => cell('b', k, q.c, q.l)).join('') +
          `<td class="aa">${s[k - 1]}${md !== undefined ? `<sub class="pep-mod">${md >= 0 ? '+' : ''}${md.toFixed(2)}</sub>` : ''}</td>` +
          yk.map(q => cell('y', yk_, q.c, q.l)).join('') + `<td>${yk_}</td>` +
          `<td>${eb && eb.peak >= 0 ? (eb.err >= 0 ? '+' : '') + eb.err.toFixed(4) : ''}</td><td>${ey && ey.peak >= 0 ? (ey.err >= 0 ? '+' : '') + ey.err.toFixed(4) : ''}</td></tr>`;
      }
      html += '</tbody></table>';
      this.table.innerHTML = `<div class="frag-wrap">${html}</div>`;
    }
  }

  // ================================================================== start
  function parseHash() {
    const o = {};
    for (const kv of (location.hash || '').replace(/^#/, '').split('&')) { if (!kv) continue; const i = kv.indexOf('='); o[decodeURIComponent(i < 0 ? kv : kv.slice(0, i))] = i < 0 ? '' : decodeURIComponent(kv.slice(i + 1)); }
    return o;
  }
  const H0 = parseHash();
  showTab(H0.tab || 'summary', false);
  if (H0.tab === 'proteins' && H0.protein !== undefined && PROT.select) {
    PROT.select(H0.protein === 'top' && T.proteins ? T.proteins.rows[PROT.dt.view[0]][col(T.proteins, 'Accession')] : H0.protein);
  }
  if (H0.tab === 'peptides' && H0.peptide !== undefined && builders.peptides.select) {
    let pep = H0.peptide;
    if (pep === 'top' && T.peptides) pep = T.peptides.rows[[...IDX.pepRow.values()].reduce((a, b) => (T.peptides.rows[b][col(T.peptides, '-10LgP')] > T.peptides.rows[a][col(T.peptides, '-10LgP')] ? b : a), 0)][col(T.peptides, 'Peptide')];
    builders.peptides.select(pep);
  }
  if (H0.tab === 'psms' && H0.psm !== undefined && builders.psms.select) {
    let i = +H0.psm;
    if (H0.psm === 'best' && T.psms) {  // highest-scoring PSM with an embedded spectrum
      const t = T.psms, is = col(t, '-10LgP'), isc = col(t, 'Scan'), isr = col(t, 'Source File');
      i = -1; t.rows.forEach((r, k) => { if (hasSpectrum(r[isr], r[isc]) && (i < 0 || r[is] > t.rows[i][is])) i = k; });
    }
    builders.psms.select(i);
  }
  if (H0.spectrum !== undefined && T.psms) {
    const k = +H0.spectrum; if (k >= 0 && k < T.psms.rows.length) openSpectrumModal(psmObj(k));
  }
  if (H0.bench !== undefined) {  // #bench: time table operations in this browser (headless checks)
    const res = [];
    const tm = (label, fn) => { const t0 = performance.now(); fn(); res.push(`${label}: ${(performance.now() - t0).toFixed(0)} ms`); };
    tm('parse spectra JSON', () => spectra());
    for (const id of ['proteins', 'peptides', 'psms']) if (!built[id]) tm(`build ${id} tab`, () => showTab(id, false));
    if (T.psms) {
      const box = h('div'); document.body.appendChild(box);
      let dt;
      tm(`new DataTable (${T.psms.rows.length} PSMs, sorted by -10lgP)`, () => { dt = new DataTable(box, { table: T.psms, sort: ['-10LgP', true] }); });
      for (const c of ['m/z', 'Peptide', 'Accession', 'RT']) tm(`sort PSMs by ${c}`, () => { dt.sortCol = col(T.psms, c); dt.sortDesc = false; dt.sort(); dt.render(); });
      tm('text filter "lk" (cold)', () => dt.setText('lk'));
      tm('text filter "elk"', () => dt.setText('elk'));
      tm('-10lgP >= 30', () => { dt.minScore = 30; dt.refresh(); });
      tm('render 1000-row page', () => { dt.setText(''); dt.pageSize = 1000; dt.render(); });
      box.remove();
      const i0 = bestPsm(T.psms.rows[0][col(T.psms, 'Peptide')]);
      if (i0 >= 0) { const b2 = h('div'); document.body.appendChild(b2); tm('render one annotated spectrum', () => new SpectrumView(b2, psmObj(i0))); b2.remove(); }
    }
    showTab('summary', false);
    app.insertBefore(h('pre', { id: 'pp-bench', class: 'card', style: 'margin:14px 18px;font-size:15px' }, 'Benchmark (this browser)\n' + res.join('\n')), app.children[1]);
  }
  window.PP = { D, spectra, theoreticalIons, matchIons, showTab };  // for debugging from the console
})();
