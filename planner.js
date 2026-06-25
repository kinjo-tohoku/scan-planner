/* filman scan planner — client-side core (no backend).
 *
 * JavaScript port of planner/planner.py + filman_ng.crystal (lattice + kinematics).
 * Same closed-form TAS kinematics; no numpy, no server. Runs in the browser and
 * in Node (for the numerical cross-check against the Python reference).
 *
 * Exposed API mirrors the server endpoints:
 *   evaluate(cfg, scan) / grid(cfg, e, hmin,hmax,kmin,kmax, n, l)
 *   to_scn(cfg, scan, scan_no) / evaluate_map(cfg, map, trim) / map_to_scn(cfg, map, trim)
 */
(function (root) {
  "use strict";

  var CONV = 180.0 / Math.PI;
  var INV_2072 = 0.48261;       // E[meV] = (ki²-kf²)/INV_2072
  var E2K = 2.07214;            // (kept for parity; unused here)
  var TWO_PI = 2.0 * Math.PI;
  var AXES = ["C1", "A1", "C2", "A2", "C3", "A3"];
  var CHUNK = 40;              // filman scan-slot limit per GO block

  var DEFAULT_LIMITS = { C2: [-180.0, 180.0], A2: [5.0, 110.0] };
  var DEFAULT_CONFIG = {
    a: 4.0, b: 4.0, c: 4.0, alpha: 90.0, beta: 90.0, gamma: 90.0,
    d_mono: 3.355, d_ana: 3.355,
    fixed: "kf", kfix: 2.662,
    sense_m: 1, sense_s: 1, sense_a: 1,
    plane_u: [1.0, 0.0, 0.0], plane_v: [0.0, 1.0, 0.0],
    limits: DEFAULT_LIMITS
  };

  function Unreachable(msg) { this.name = "Unreachable"; this.message = msg; }
  Unreachable.prototype = Object.create(Error.prototype);
  Unreachable.prototype.constructor = Unreachable;

  // --- tiny linear algebra (3-vectors / 3x3) -----------------------------
  function dot(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
  function norm(a) { return Math.sqrt(dot(a, a)); }
  function cross(a, b) {
    return [a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]];
  }
  function matVec(M, v) {
    return [M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2],
            M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2],
            M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2]];
  }
  function inv3(m) {                // 3x3 inverse (for Q_cart → hkl)
    var a = m[0][0], b = m[0][1], c = m[0][2],
        d = m[1][0], e = m[1][1], f = m[1][2],
        g = m[2][0], h = m[2][1], i = m[2][2];
    var A = e * i - f * h, B = -(d * i - f * g), C = d * h - e * g;
    var det = a * A + b * B + c * C, id = 1.0 / det;
    return [[A * id, (c * h - b * i) * id, (b * f - c * e) * id],
            [B * id, (a * i - c * g) * id, (c * d - a * f) * id],
            [C * id, (b * g - a * h) * id, (a * e - b * d) * id]];
  }
  function clampcos(x) { return Math.max(-1.0, Math.min(1.0, x)); }
  var rad = function (d) { return d * Math.PI / 180.0; };

  // --- lattice → B matrix (Busing & Levy 1967, 2π convention) -------------
  function bMatrix(cfg) {
    var a = +cfg.a, b = +cfg.b, c = +cfg.c;
    var al = rad(cfg.alpha == null ? 90 : +cfg.alpha);
    var be = rad(cfg.beta == null ? 90 : +cfg.beta);
    var ga = rad(cfg.gamma == null ? 90 : +cfg.gamma);
    var ca = Math.cos(al), cb = Math.cos(be), cg = Math.cos(ga);
    var sa = Math.sin(al), sb = Math.sin(be), sg = Math.sin(ga);
    var V = a * b * c * Math.sqrt(Math.max(0.0,
      1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg));
    var astar = TWO_PI * b * c * sa / V;
    var bstar = TWO_PI * a * c * sb / V;
    var cstar = TWO_PI * a * b * sg / V;
    var cbs = (ca * cg - cb) / (sa * sg);  // cos β*
    var cgs = (ca * cb - cg) / (sa * sb);  // cos γ*
    var sgs = Math.sqrt(Math.max(0.0, 1.0 - cgs * cgs));
    var sbs = Math.sqrt(Math.max(0.0, 1.0 - cbs * cbs));
    return [
      [astar, bstar * cgs, cstar * cbs],
      [0.0, bstar * sgs, -cstar * sbs * ca],
      [0.0, 0.0, TWO_PI / c]
    ];
  }

  // --- spectrometer (frame + fixed config) -------------------------------
  function buildSpec(cfg) {
    var B = bMatrix(cfg);
    var pu = cfg.plane_u || [1, 0, 0], pv = cfg.plane_v || [0, 1, 0];
    var u = matVec(B, [+pu[0], +pu[1], +pu[2]]);
    var v = matVec(B, [+pv[0], +pv[1], +pv[2]]);
    var nu = norm(u);
    if (nu < 1e-9) throw new Error("plane_u must be a non-zero reflection");
    var e1 = [u[0] / nu, u[1] / nu, u[2] / nu];
    var vd = dot(v, e1);
    var vp = [v[0] - vd * e1[0], v[1] - vd * e1[1], v[2] - vd * e1[2]];
    var nv = norm(vp);
    if (nv < 1e-9) throw new Error("plane_u and plane_v are parallel (no plane)");
    var e2 = [vp[0] / nv, vp[1] / nv, vp[2] / nv];
    var n = cross(e1, e2);
    return {
      B: B, e1: e1, e2: e2, n: n,
      d_mono: +cfg.d_mono, d_ana: +cfg.d_ana,
      fixed: cfg.fixed || "kf", kfix: +cfg.kfix,
      sense_m: cfg.sense_m == null ? 1 : (+cfg.sense_m),
      sense_s: cfg.sense_s == null ? 1 : (+cfg.sense_s),
      sense_a: cfg.sense_a == null ? 1 : (+cfg.sense_a),
      tol: 1e-4
    };
  }

  function kiKf(spec, e) {
    if (spec.fixed === "ki") {
      var ki = spec.kfix, kf2 = ki * ki - INV_2072 * e;
      if (kf2 < 0) throw new Unreachable("energy transfer too large: E=" + e);
      return [ki, Math.sqrt(kf2)];
    } else if (spec.fixed === "kf") {
      var kf = spec.kfix, ki2 = kf * kf + INV_2072 * e;
      if (ki2 < 0) throw new Unreachable("energy transfer too large: E=" + e);
      return [Math.sqrt(ki2), kf];
    }
    throw new Error("fixed must be 'ki' or 'kf'");
  }

  // 6-axis angles for (h,k,l), E[meV]. Throws Unreachable when impossible.
  function angles(spec, hkl, e) {
    var kk = kiKf(spec, e), ki = kk[0], kf = kk[1];

    var tau_m = TWO_PI / spec.d_mono, tau_a = TWO_PI / spec.d_ana;
    var sm = tau_m / (2.0 * ki), sa = tau_a / (2.0 * kf);
    if (Math.abs(sm) > 1.0 || Math.abs(sa) > 1.0)
      throw new Unreachable("unobtainable takeoff angle for mono/ana");
    var c1 = CONV * Math.asin(sm) * spec.sense_m;
    var a1 = 2.0 * c1;
    var c3 = CONV * Math.asin(sa) * spec.sense_a;
    var a3 = 2.0 * c3;

    var Q = matVec(spec.B, [+hkl[0], +hkl[1], +hkl[2]]);
    var Qmag = norm(Q);
    if (Qmag <= 1e-6) return { C1: c1, A1: a1, C2: 0.0, A2: 0.0, C3: c3, A3: a3 };

    var qn = dot(Q, spec.n);
    if (Math.abs(qn) > spec.tol * Math.max(1.0, Qmag))
      throw new Unreachable("Q=(" + hkl.join(",") + ") out of scattering plane "
        + "(out-of-plane component " + qn.toPrecision(4) + " Å⁻¹); "
        + "unreachable on an in-plane instrument");

    var q1 = dot(Q, spec.e1);
    var q2 = dot(Q, spec.e2) * spec.sense_m;

    if ((ki + Qmag) < kf || (kf + Qmag) < ki || (ki + kf) < Qmag)
      throw new Unreachable("scattering triangle will not close");

    var cos2t = (kf * kf + ki * ki - Qmag * Qmag) / (2.0 * ki * kf);
    var a2 = CONV * Math.acos(clampcos(cos2t)) * spec.sense_s;

    var cos_alf = (Qmag * Qmag + ki * ki - kf * kf) / (2.0 * ki * Qmag);
    var alf = CONV * Math.acos(clampcos(cos_alf));
    var bet = CONV * Math.atan2(q2, q1);
    if (bet < -135.0) bet += 360.0;
    var c2 = 90.0 - alf + bet;
    if (spec.sense_s < 0) c2 = 180.0 - c2;
    if (c2 < -155.0) c2 += 360.0;
    if (c2 > 205.0) c2 -= 360.0;

    return { C1: c1, A1: a1, C2: c2, A2: a2, C3: c3, A3: a3 };
  }

  var r3 = function (x) { return Math.round(x * 1e3) / 1e3; };
  var r4 = function (x) { return Math.round(x * 1e4) / 1e4; };
  // |Q|, d-spacing, and in-plane Cartesian components Qx=Q·e1, Qy=Q·e2 [Å⁻¹].
  function qd(spec, hkl) {
    var Q = matVec(spec.B, [+hkl[0], +hkl[1], +hkl[2]]);
    var q = norm(Q);
    return { q: r4(q), d: q > 1e-9 ? r4(TWO_PI / q) : null,
             qx: r4(dot(Q, spec.e1)), qy: r4(dot(Q, spec.e2)) };
  }

  function build(cfg) {
    var spec = buildSpec(cfg);
    var src = cfg.limits || DEFAULT_LIMITS, lim = {};
    for (var k in src) if (Object.prototype.hasOwnProperty.call(src, k))
      lim[k] = [+src[k][0], +src[k][1]];
    return { spec: spec, limits: lim };
  }

  // Limit test. C-axes are full rotations: an angle and its ±360 equivalents are
  // the same orientation, so accept if any equivalent falls within [lo,hi] (a 360°
  // limit then reaches every azimuth). A-axes (scattering 2θ) use the plain test.
  function within(ax, v, lo, hi) {
    if (ax.charAt(0) === "C")
      return (lo <= v && v <= hi) || (lo <= v - 360 && v - 360 <= hi)
          || (lo <= v + 360 && v + 360 <= hi);
    return lo <= v && v <= hi;
  }

  function checkPoint(spec, limits, hkl, e) {
    var ang;
    try { ang = angles(spec, hkl, e); }
    catch (ex) {
      if (ex instanceof Unreachable) return { reachable: false, reason: ex.message };
      throw ex;
    }
    var vals = {};
    for (var i = 0; i < AXES.length; i++) vals[AXES[i]] = r3(ang[AXES[i]]);
    for (var ax in limits) if (Object.prototype.hasOwnProperty.call(limits, ax)) {
      if (vals[ax] !== undefined) {
        var lo = limits[ax][0], hi = limits[ax][1], v = vals[ax];
        if (!within(ax, v, lo, hi))
          return {
            reachable: false,
            reason: ax + "=" + v.toFixed(2) + "° outside [" + lo + "," + hi + "]",
            angles: vals
          };
      }
    }
    return { reachable: true, angles: vals };
  }

  function scanPoints(scan) {
    var s = scan.start.map(Number), d = scan.step.map(Number), n = scan.npts | 0;
    var out = [];
    for (var i = 0; i < n; i++)
      out.push([s[0] + i * d[0], s[1] + i * d[1], s[2] + i * d[2], s[3] + i * d[3]]);
    return out;
  }

  function evaluate(cfg, scan) {
    var bl = build(cfg), pts = [], n_ok = 0;
    var sp = scanPoints(scan);
    for (var i = 0; i < sp.length; i++) {
      var p = sp[i], r = checkPoint(bl.spec, bl.limits, [p[0], p[1], p[2]], p[3]);
      var qq = qd(bl.spec, [p[0], p[1], p[2]]);
      r.h = r4(p[0]); r.k = r4(p[1]); r.l = r4(p[2]); r.E = r4(p[3]);
      r.q = qq.q; r.d = qq.d; r.qx = qq.qx; r.qy = qq.qy;
      if (r.reachable) n_ok++;
      pts.push(r);
    }
    return { points: pts, n_reachable: n_ok, n_total: pts.length };
  }

  function grid(cfg, e, hmin, hmax, kmin, kmax, n, l) {
    var bl = build(cfg);
    n = Math.max(2, n | 0);
    var cells = [], hs = [], ks = [];
    for (var i = 0; i < n; i++) {
      var h = hmin + (hmax - hmin) * i / (n - 1), row = [];
      for (var j = 0; j < n; j++) {
        var k = kmin + (kmax - kmin) * j / (n - 1);
        row.push(checkPoint(bl.spec, bl.limits, [h, k, l], e).reachable ? 1 : 0);
      }
      cells.push(row);
    }
    for (var ii = 0; ii < n; ii++) hs.push(r4(hmin + (hmax - hmin) * ii / (n - 1)));
    for (var jj = 0; jj < n; jj++) ks.push(r4(kmin + (kmax - kmin) * jj / (n - 1)));
    return { h: hs, k: ks, z: cells };
  }

  // Reachable-region map meshed directly in the in-plane Cartesian (Qx,Qy) [Å⁻¹].
  // Each grid Q is mapped back to (h,k,l) via B⁻¹ and checked. The mesh spans the
  // full reachable disk (|Q| up to ki+kf) unless qmax is given.
  function gridQ(cfg, e, n, qmax) {
    var bl = build(cfg), Binv = inv3(bl.spec.B);
    var e1 = bl.spec.e1, e2 = bl.spec.e2;
    if (!(qmax > 0)) {
      try { var kk = kiKf(bl.spec, e); qmax = (kk[0] + kk[1]) * 1.06; }
      catch (ex) { qmax = 2.0 * bl.spec.kfix * 1.06; }
    }
    n = Math.max(2, n | 0);
    var cells = [], qxs = [], qys = [];
    for (var i = 0; i < n; i++) {
      var qx = -qmax + 2 * qmax * i / (n - 1), row = [];
      for (var j = 0; j < n; j++) {
        var qy = -qmax + 2 * qmax * j / (n - 1);
        var Qc = [qx * e1[0] + qy * e2[0], qx * e1[1] + qy * e2[1], qx * e1[2] + qy * e2[2]];
        var hkl = matVec(Binv, Qc);
        row.push(checkPoint(bl.spec, bl.limits, hkl, e).reachable ? 1 : 0);
      }
      cells.push(row);
    }
    for (var ii = 0; ii < n; ii++) qxs.push(r4(-qmax + 2 * qmax * ii / (n - 1)));
    for (var jj = 0; jj < n; jj++) qys.push(r4(-qmax + 2 * qmax * jj / (n - 1)));
    return { qx: qxs, qy: qys, z: cells };
  }

  // In-plane integer (h,k,l) Bragg reflections with |Q| ≤ qmax, as (Qx,Qy) markers.
  // No space-group absences (P1: every integer reflection).
  function reflections(cfg, qmax) {
    var bl = build(cfg), spec = bl.spec;
    if (!(qmax > 0)) qmax = 6.0;
    var col = function (j) { return [spec.B[0][j], spec.B[1][j], spec.B[2][j]]; };
    var minr = Math.min(norm(col(0)), norm(col(1)), norm(col(2)));
    var hmax = Math.min(25, Math.ceil(qmax / Math.max(minr, 1e-6)) + 1);
    var out = [], tol = 1e-4;
    for (var h = -hmax; h <= hmax; h++)
      for (var k = -hmax; k <= hmax; k++)
        for (var l = -hmax; l <= hmax; l++) {
          if (!h && !k && !l) continue;
          var Q = matVec(spec.B, [h, k, l]), Qmag = norm(Q);
          if (Qmag > qmax) continue;
          if (Math.abs(dot(Q, spec.n)) > tol * Math.max(1, Qmag)) continue;
          out.push({ h: h, k: k, l: l, qx: r4(dot(Q, spec.e1)), qy: r4(dot(Q, spec.e2)) });
        }
    return out;
  }

  // --- .scn text generation ----------------------------------------------
  function gfmt(x) {                 // mimic Python's %g for scan coordinates
    x = +x;
    if (x === 0) return "0";
    if (!isFinite(x)) return String(x);
    var s = x.toPrecision(6);
    if (s.indexOf("e") < 0 && s.indexOf(".") >= 0) s = s.replace(/\.?0+$/, "");
    return s;
  }
  function anaName(d) {
    return Math.abs(d - 3.355) < 0.01 ? "PG002"
      : Math.abs(d - 1.6775) < 0.01 ? "PG004" : "d=" + d;
  }
  function headerComment(cfg, prefix) {
    return '" ' + (prefix || "") + "a=" + cfg.a + " b=" + cfg.b + " c=" + cfg.c
      + " fixed-" + (cfg.fixed || "kf") + "=" + cfg.kfix
      + " ana=" + anaName(+cfg.d_ana);
  }
  function scanLine(no, start, step, npts, monitor, nt) {
    var h0 = +start[0], k0 = +start[1], e0 = +start[3];
    var dh = +step[0], dk = +step[1], de = +step[3];
    var f = ["NS=" + (no | 0), "HS=" + gfmt(h0), "KS=" + gfmt(k0), "ES=" + gfmt(e0),
             "DE=" + gfmt(de), "DH=" + gfmt(dh), "DK=" + gfmt(dk),
             "NP=" + (npts | 0), "MN=" + (monitor | 0)];
    if ((nt | 0) > 0) f.push("NT=" + (nt | 0));
    return f.join(",");
  }
  function to_scn(cfg, scan, scan_no) {
    scan_no = scan_no || 1;
    var l0 = +scan.start[2], dl = +scan.step[2];
    var L = [headerComment(cfg),
      scanLine(scan_no, scan.start, scan.step, scan.npts,
        scan.monitor == null ? 10000 : scan.monitor, scan.nt || 0),
      "GO " + scan_no];
    if (Math.abs(l0) > 1e-6 || Math.abs(dl) > 1e-6)
      L.push('" NOTE: l!=0 leaves the in-plane assumption (out-of-plane is unreachable)');
    return L.join("\n");
  }

  // --- maps ---------------------------------------------------------------
  function mapLines(md) {
    var start = md.start.map(Number), step = md.step.map(Number), outer = md.outer.map(Number);
    var nlines = Math.max(1, md.nlines | 0), npts = md.npts | 0;
    var mon = md.monitor == null ? 10000 : (md.monitor | 0), nt = md.nt | 0, out = [];
    for (var i = 0; i < nlines; i++)
      out.push({
        start: [start[0] + i * outer[0], start[1] + i * outer[1],
                start[2] + i * outer[2], start[3] + i * outer[3]],
        step: step.slice(), npts: npts, monitor: mon, nt: nt
      });
    return out;
  }
  function mapScans(cfg, md, trim) {
    var lines = mapLines(md);
    if (!trim) return lines;
    var bl = build(cfg), out = [];
    for (var li = 0; li < lines.length; li++) {
      var ln = lines[li], pts = scanPoints(ln);
      var ok = pts.map(function (p) {
        return checkPoint(bl.spec, bl.limits, [p[0], p[1], p[2]], p[3]).reachable;
      });
      var i = 0, n = ok.length;
      while (i < n) {
        if (!ok[i]) { i++; continue; }
        var j = i; while (j < n && ok[j]) j++;
        out.push({ start: pts[i].slice(), step: ln.step.slice(),
                   npts: j - i, monitor: ln.monitor, nt: ln.nt });
        i = j;
      }
    }
    return out;
  }
  function evaluate_map(cfg, md, trim) {
    var bl = build(cfg), scans = mapScans(cfg, md, trim);
    var out = [], n_ok = 0, n_tot = 0;
    for (var li = 0; li < scans.length; li++) {
      var ln = scans[li], pts = [], ok = 0, sp = scanPoints(ln);
      for (var s = 0; s < sp.length; s++) {
        var p = sp[s], r = checkPoint(bl.spec, bl.limits, [p[0], p[1], p[2]], p[3]);
        var qq = qd(bl.spec, [p[0], p[1], p[2]]);
        r.h = r4(p[0]); r.k = r4(p[1]); r.l = r4(p[2]); r.E = r4(p[3]); r.line = li;
        r.q = qq.q; r.d = qq.d; r.qx = qq.qx; r.qy = qq.qy;
        if (r.reachable) ok++;
        pts.push(r);
      }
      out.push({ line: li, start: ln.start.map(r4), points: pts,
                 n_reachable: ok, n_total: pts.length });
      n_ok += ok; n_tot += pts.length;
    }
    var res = { lines: out, n_reachable: n_ok, n_total: n_tot, n_lines: out.length };
    if (trim) {
      var dropped = [], all = mapLines(md);
      for (var a = 0; a < all.length; a++) {
        var ap = scanPoints(all[a]);
        for (var q = 0; q < ap.length; q++) {
          var pp = ap[q];
          if (!checkPoint(bl.spec, bl.limits, [pp[0], pp[1], pp[2]], pp[3]).reachable) {
            var dq = qd(bl.spec, [pp[0], pp[1], pp[2]]);
            dropped.push({ h: r4(pp[0]), k: r4(pp[1]), l: r4(pp[2]), E: r4(pp[3]),
                           qx: dq.qx, qy: dq.qy });
          }
        }
      }
      res.dropped = dropped;
    }
    return res;
  }
  function map_to_scn(cfg, md, trim) {
    var scans = mapScans(cfg, md, trim), n = scans.length;
    var L = [headerComment(cfg, "map " + n + " lines: ")];
    if (n === 0) { L.push('" no reachable points; revise the range'); return L.join("\n"); }
    var nblocks = Math.floor((n + CHUNK - 1) / CHUNK);
    if (nblocks > 1)
      L.push('" ' + n + " lines > " + CHUNK + "/GO: emitted as " + nblocks
        + " GO blocks (NS restarts at 1 in each block)");
    for (var c0 = 0; c0 < n; c0 += CHUNK) {
      var chunk = scans.slice(c0, c0 + CHUNK);
      for (var i = 0; i < chunk.length; i++)
        L.push(scanLine(i + 1, chunk[i].start, chunk[i].step, chunk[i].npts,
          chunk[i].monitor, chunk[i].nt));
      var k = chunk.length;
      L.push(k === 1 ? "GO 1" : "GO 1-" + k);
    }
    var oop = scans.some(function (ln) { return Math.abs(ln.start[2]) > 1e-6; })
      || Math.abs(+md.step[2]) > 1e-6 || Math.abs(+md.outer[2]) > 1e-6;
    if (oop)
      L.push('" NOTE: l!=0 leaves the in-plane assumption (out-of-plane is unreachable)');
    return L.join("\n");
  }

  var API = {
    DEFAULT_CONFIG: DEFAULT_CONFIG, DEFAULT_LIMITS: DEFAULT_LIMITS, CHUNK: CHUNK,
    Unreachable: Unreachable,
    build: build, anglesRaw: function (cfg, hkl, e) { return angles(buildSpec(cfg), hkl, e); },
    check_point: function (cfg, hkl, e) { var b = build(cfg); return checkPoint(b.spec, b.limits, hkl, e); },
    evaluate: evaluate, grid: grid, gridQ: gridQ, reflections: reflections, to_scn: to_scn,
    evaluate_map: evaluate_map, map_to_scn: map_to_scn
  };

  if (typeof module !== "undefined" && module.exports) module.exports = API;
  else root.Planner = API;
})(typeof self !== "undefined" ? self : this);
