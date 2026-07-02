/* filman scan planner — client-side core (no backend).
 *
 * Self-contained triple-axis-spectrometer engine: lattice geometry, closed-form
 * kinematics, reachability, Cooper–Nathans resolution and .scn export. No numpy,
 * no server — runs directly in the browser (and in Node, for the self-tests).
 *
 * Public API:
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
    centering: "P", km: [0.0, 0.0, 0.0],
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
    var disc = 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg;
    if (!(disc > 1e-12)) throw new Error("lattice angles are geometrically impossible (cell volume ≤ 0)");
    var V = a * b * c * Math.sqrt(disc);
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
    if (!(nu >= 1e-9)) throw new Error("scattering plane: u must be a non-zero reflection");
    var e1 = [u[0] / nu, u[1] / nu, u[2] / nu];
    var vd = dot(v, e1);
    var vp = [v[0] - vd * e1[0], v[1] - vd * e1[1], v[2] - vd * e1[2]];
    var nv = norm(vp);
    if (!(nv >= 1e-9)) throw new Error("scattering plane: v must be non-zero and not parallel to u");
    var e2 = [vp[0] / nv, vp[1] / nv, vp[2] / nv];
    var n = cross(e1, e2);
    return {
      B: B, e1: e1, e2: e2, n: n,
      d_mono: +cfg.d_mono, d_ana: +cfg.d_ana,
      fixed: cfg.fixed || "kf", kfix: +cfg.kfix,
      sense_m: cfg.sense_m == null ? 1 : (+cfg.sense_m),
      sense_s: cfg.sense_s == null ? 1 : (+cfg.sense_s),
      sense_a: cfg.sense_a == null ? 1 : (+cfg.sense_a),
      magnet: (cfg.magnet && cfg.magnet.on && cfg.magnet.ref) ? {
        ref: [+cfg.magnet.ref[0], +cfg.magnet.ref[1], +cfg.magnet.ref[2]],
        half: (cfg.magnet.half_deg != null ? +cfg.magnet.half_deg : 60) * Math.PI / 180,
        a2max: 110 * Math.PI / 180
      } : null,
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

  // window test: beam angle x [rad] passes if within ±d of a window. Two opposite windows
  // (front at 0, back at ±180°), so a beam clears if |x| < d (front) or |x| > π−d (back).
  function winNear(x, d) {
    x = Math.abs(x - 2 * Math.PI * Math.round(x / (2 * Math.PI)));   // |wrap to (−π,π]| ∈ [0,π]
    return x < d || x > Math.PI - d;
  }
  // Sample-environment magnet that rotates with the sample. Two opposite windows (the mounting
  // reference direction ± 180°), each of half-width `half`. The incident AND the scattered beam
  // must EACH pass through one of the two windows → the accessible region is two symmetric lobes
  // about the ref direction. Elastic geometry (extends addon/calculate_range.py, which only let
  // the incident use the front window). Returns true if the magnet BLOCKS Q.
  // mag = { ref:[h,k,l], half:rad, a2max:rad }. ki = incident wavevector at this point.
  function magnetBlocks(spec, hkl, ki) {
    var mag = spec.magnet;
    if (!mag || !(ki > 0)) return false;
    var Q = matVec(spec.B, [+hkl[0], +hkl[1], +hkl[2]]), qm = norm(Q);
    if (qm < 1e-9) return false;                       // |Q|=0 — beams undefined, don't block
    var rd = matVec(spec.B, mag.ref);
    var r1 = dot(rd, spec.e1), r2 = dot(rd, spec.e2);
    var rt = r1; r1 = -r2; r2 = rt;                    // rotate +90° in-plane: `ref` is the
    var rn = Math.sqrt(r1 * r1 + r2 * r2);             //   ACCESS direction (reflection you align
    if (rn < 1e-9) return false;                       //   to / centre the window on), not the beam
    var q1 = dot(Q, spec.e1), q2 = dot(Q, spec.e2);
    var sgn = (q1 * r2 - q2 * r1) > 0 ? 1 : -1;        // side of Q relative to ref (kp = ref rot −90°)
    var omega = sgn * Math.acos(clampcos((q1 * r1 + q2 * r2) / (qm * rn))) - Math.PI / 2;
    var theta = Math.asin(Math.min(1.0, qm / (2.0 * ki)));   // elastic half-scattering angle
    var C2 = omega + theta, A2 = 2.0 * theta, d = mag.half;
    if (Math.abs(A2) >= mag.a2max) return true;        // beyond the elastic 2θ ceiling
    return !(winNear(C2, d) && winNear(A2 - C2, d));   // incident AND scattered each clear a window
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
    if (spec.magnet && magnetBlocks(spec, hkl, kiKf(spec, e)[0]))
      return { reachable: false, reason: "blocked by magnet window", angles: vals };
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

  // Lattice-centering reflection condition (which reflections are present).
  function allowedCentering(cen, h, k, l) {
    switch (cen) {
      case "I": return (h + k + l) % 2 === 0;
      case "F": return (h + k) % 2 === 0 && (k + l) % 2 === 0;  // all same parity
      case "C": return (h + k) % 2 === 0;
      case "A": return (k + l) % 2 === 0;
      case "B": return (h + l) % 2 === 0;
      case "R": return ((((-h + k + l) % 3) + 3) % 3) === 0;    // obverse on hex axes
      default:  return true;                                    // P (all)
    }
  }
  // In-plane integer (h,k,l) Bragg reflections with |Q| ≤ qmax, as (Qx,Qy) markers.
  // Lattice-centering absences applied (cfg.centering, default P); no glide/screw rules.
  function reflections(cfg, qmax) {
    var bl = build(cfg), spec = bl.spec, cen = (cfg.centering || "P").toUpperCase();
    if (!(qmax > 0)) qmax = 6.0;
    var col = function (j) { return [spec.B[0][j], spec.B[1][j], spec.B[2][j]]; };
    var minr = Math.min(norm(col(0)), norm(col(1)), norm(col(2)));
    var hmax = Math.min(25, Math.ceil(qmax / Math.max(minr, 1e-6)) + 1);
    var out = [], tol = 1e-4;
    for (var h = -hmax; h <= hmax; h++)
      for (var k = -hmax; k <= hmax; k++)
        for (var l = -hmax; l <= hmax; l++) {
          if (!h && !k && !l) continue;
          if (!allowedCentering(cen, h, k, l)) continue;
          var Q = matVec(spec.B, [h, k, l]), Qmag = norm(Q);
          if (Qmag > qmax) continue;
          if (Math.abs(dot(Q, spec.n)) > tol * Math.max(1, Qmag)) continue;
          out.push({ h: h, k: k, l: l, qx: r4(dot(Q, spec.e1)), qy: r4(dot(Q, spec.e2)) });
        }
    return out;
  }

  // Magnetic satellites at G ± k_m around centering-allowed nuclear positions G.
  // Empty when k_m = 0. s = +1/-1 marks the satellite sign.
  function satellites(cfg, qmax) {
    var km = cfg.km || [0, 0, 0];
    if (!(Math.abs(km[0]) > 1e-9 || Math.abs(km[1]) > 1e-9 || Math.abs(km[2]) > 1e-9)) return [];
    var bl = build(cfg), spec = bl.spec, cen = (cfg.centering || "P").toUpperCase();
    if (!(qmax > 0)) qmax = 6.0;
    var col = function (j) { return [spec.B[0][j], spec.B[1][j], spec.B[2][j]]; };
    var minr = Math.min(norm(col(0)), norm(col(1)), norm(col(2)));
    var pad = norm(matVec(spec.B, [+km[0], +km[1], +km[2]]));
    var hmax = Math.min(25, Math.ceil((qmax + pad) / Math.max(minr, 1e-6)) + 1);
    var out = [], tol = 1e-4;
    for (var h = -hmax; h <= hmax; h++)
      for (var k = -hmax; k <= hmax; k++)
        for (var l = -hmax; l <= hmax; l++) {
          if (!allowedCentering(cen, h, k, l)) continue;
          for (var s = -1; s <= 1; s += 2) {
            var Q = matVec(spec.B, [h + s * km[0], k + s * km[1], l + s * km[2]]);
            var Qmag = norm(Q);
            if (Qmag < 1e-6 || Qmag > qmax) continue;
            if (Math.abs(dot(Q, spec.n)) > tol * Math.max(1, Qmag)) continue;
            out.push({ h: h, k: k, l: l, s: s, qx: r4(dot(Q, spec.e1)), qy: r4(dot(Q, spec.e2)) });
          }
        }
    return out;
  }

  function clipHalfplane(poly, g) {           // keep x with x·g ≤ |g|²/2
    var c = (g[0] * g[0] + g[1] * g[1]) / 2, out = [];
    for (var i = 0; i < poly.length; i++) {
      var A = poly[i], Bp = poly[(i + 1) % poly.length];
      var da = A[0] * g[0] + A[1] * g[1] - c, db = Bp[0] * g[0] + Bp[1] * g[1] - c;
      if (da <= 1e-12) out.push(A);
      if ((da < 0 && db > 0) || (da > 0 && db < 0)) {
        var t = da / (da - db);
        out.push([A[0] + t * (Bp[0] - A[0]), A[1] + t * (Bp[1] - A[1])]);
      }
    }
    return out;
  }
  // 1st Brillouin zone (Wigner–Seitz cell of the in-plane reciprocal lattice) as a
  // polygon of (qx,qy) vertices. Uses the centering-allowed (primitive) lattice.
  function brillouinZone(cfg) {
    var bl = build(cfg), spec = bl.spec, cen = (cfg.centering || "P").toUpperCase();
    var vs = [];
    for (var h = -3; h <= 3; h++)
      for (var k = -3; k <= 3; k++)
        for (var l = -3; l <= 3; l++) {
          if (!h && !k && !l) continue;
          if (!allowedCentering(cen, h, k, l)) continue;
          var Q = matVec(spec.B, [h, k, l]);
          if (Math.abs(dot(Q, spec.n)) > 1e-4 * Math.max(1, norm(Q))) continue;
          vs.push([dot(Q, spec.e1), dot(Q, spec.e2)]);
        }
    if (!vs.length) return [];
    vs.sort(function (a, b) { return (a[0] * a[0] + a[1] * a[1]) - (b[0] * b[0] + b[1] * b[1]); });
    var g1 = vs[0], g2 = null;
    for (var i = 1; i < vs.length; i++) {
      if (Math.abs(g1[0] * vs[i][1] - g1[1] * vs[i][0]) > 1e-6) { g2 = vs[i]; break; }
    }
    if (!g2) return [];
    var R = 3 * Math.sqrt(Math.max(g1[0] * g1[0] + g1[1] * g1[1], g2[0] * g2[0] + g2[1] * g2[1]));
    var poly = [[-R, -R], [R, -R], [R, R], [-R, R]];
    for (var m = -2; m <= 2; m++)
      for (var n = -2; n <= 2; n++) {
        if (!m && !n) continue;
        poly = clipHalfplane(poly, [m * g1[0] + n * g2[0], m * g1[1] + n * g2[1]]);
        if (!poly.length) return [];
      }
    poly = poly.filter(function (p, i) {        // drop consecutive duplicate vertices
      var q = poly[(i + poly.length - 1) % poly.length];
      return Math.abs(p[0] - q[0]) > 1e-6 || Math.abs(p[1] - q[1]) > 1e-6;
    });
    return poly.map(function (p) { return [r4(p[0]), r4(p[1])]; });
  }

  // === resolution (Cooper-Nathans, port of resolution_reduced_grobal.RESELP) ===
  // arcmin FWHM → radian σ  (TASIN DEGRAD): π/180 / sqrt(8 ln2) / 60
  var DEGRAD = Math.PI / 180.0 / Math.sqrt(8.0 * Math.log(2.0)) / 60.0;

  // Resolution matrix RM (3×3, frame Q∥,Q⊥,E) for |Q| and energy transfer E0.
  // p: {AKI, TM, TA, ALZ,ALM,ALA,AL3, VBETZ,VBET1,VBET2,VBET3, ETAM,ETAA,ETAS (all
  // already in radian σ), EPM, EP, ERL, IR1}. Returns null if geometry impossible.
  function resMatrix(Q, E0, p) {
    var F = 1.0 / 2.072, AOM = E0 * F;
    var AKI = p.AKI, AKF;
    if (AKI < 0.0) { AKF = -AKI; var t = AKI * AKI + AOM; if (t < 0) return null; AKI = Math.sqrt(t); }
    else { var t2 = AKI * AKI - AOM; if (t2 < 0) return null; AKF = Math.sqrt(t2); }
    var ALAM = AKI / AKF;
    var BE = -(Q * Q - 2.0 * AKI * AKI + AOM) / (2.0 * AKI * AKF);
    if (1 - BE * BE < 0) return null; var AL = Math.sqrt(1 - BE * BE);
    var B = -(Q * Q - AOM) / (2.0 * Q * AKF);
    if (1 - B * B < 0) return null; var AA = Math.sqrt(1 - B * B);
    var ALP = [B / AL, -AA / AL, 1.0 / (AL * 2.0 * AKF)];
    var SB = (Q * Q + AOM) / (2.0 * Q * AKI);
    if (1 - SB * SB < 0) return null; var SA = Math.sqrt(1 - SB * SB);
    var BET = [SB / AL, -SA / AL, BE / (AL * 2.0 * AKF)];
    var GAM = [0.0, 0.0, -1.0 / (2.0 * AKF)];
    var tm2 = AKI * AKI - (p.TM / 2) * (p.TM / 2), ta2 = AKF * AKF - (p.TA / 2) * (p.TA / 2);
    if (tm2 <= 0 || ta2 <= 0) return null;
    var TOM = p.TM / (2.0 * Math.sqrt(tm2)) * p.EPM;
    var TOA = (p.EP * p.TA / 2.0) / Math.sqrt(ta2);
    var A1 = TOM / (AKI * p.ETAM), A2 = 1.0 / (AKI * p.ETAM), A3 = 1.0 / (AKI * p.ALM),
        A4 = 1.0 / (AKF * p.ALA), A5 = TOA / (AKF * p.ETAA), A6 = -1.0 / (AKF * p.ETAA),
        A7 = 2.0 * TOM / (AKI * p.ALZ), A8 = 1.0 / (AKI * p.ALZ),
        A9 = 2.0 * TOA / (p.AL3 * AKF), A10 = -1.0 / (p.AL3 * AKF);
    if (p.IR1 === 1) { A6 = 0; A9 = 0; A10 = 0; }
    var B0 = A1 * A2 + A7 * A8, B1 = A2 * A2 + A3 * A3 + A8 * A8,
        B2 = A4 * A4 + A6 * A6 + A10 * A10, B3 = A5 * A5 + A9 * A9, B4 = A5 * A6 + A9 * A10;
    var C = -(ALAM - BE) / AL, E = -(BE * ALAM - 1.0) / AL;
    var AP = A1 * A1 + 2.0 * B0 * C + B1 * C * C + B2 * E * E + B3 * ALAM * ALAM + 2.0 * B4 * ALAM * E + A7 * A7;
    var D0 = B1 - (B0 + B1 * C) * (B0 + B1 * C) / AP;
    var D1 = B2 - (B2 * E + B4 * ALAM) * (B2 * E + B4 * ALAM) / AP;
    var D2 = B3 - (B3 * ALAM + B4 * E) * (B3 * ALAM + B4 * E) / AP;
    var D3 = 2.0 * B4 - 2.0 / AP * (B2 * E + B4 * ALAM) * (B3 * ALAM + B4 * E);
    var D4 = -2.0 / AP * (B0 + B1 * C) * (B2 * E + B4 * ALAM);
    var D5 = -2.0 / AP * (B0 + B1 * C) * (B3 * ALAM + B4 * E);
    var A = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (var i = 0; i < 3; i++) for (var j = 0; j < 3; j++) {
      var v = D0 * ALP[i] * ALP[j] + D1 * BET[i] * BET[j] + D2 * GAM[i] * GAM[j]
        + 0.5 * D3 * (BET[i] * GAM[j] + BET[j] * GAM[i])
        + 0.5 * D4 * (ALP[i] * BET[j] + ALP[j] * BET[i])
        + 0.5 * D5 * (ALP[i] * GAM[j] + ALP[j] * GAM[i]);
      if (i === 2) v *= F; if (j === 2) v *= F;
      A[i][j] = v;
    }
    var SNM = p.TM / (2.0 * AKI), SNA = p.TA / (2.0 * AKF);
    var A11 = 1.0 / (Math.pow(2 * SNM * p.ETAM, 2) + p.VBETZ * p.VBETZ) / (AKI * AKI) + Math.pow(1 / p.VBET1 / AKI, 2);
    var A12 = 1.0 / (Math.pow(2 * SNA * p.ETAA, 2) + p.VBET3 * p.VBET3) / (AKF * AKF) + Math.pow(1 / p.VBET2 / AKF, 2);
    var GAMS = (Q * p.ETAS) * (Q * p.ETAS), DENOM = A[1][1] * GAMS + 1.0;
    var RM = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (var ii = 0; ii < 3; ii++) for (var jj = 0; jj < 3; jj++)
      RM[ii][jj] = A[ii][jj] - A[ii][1] * A[jj][1] * GAMS / DENOM;
    if (p.ERL !== 1) { RM[0][1] = -RM[0][1]; RM[1][0] = -RM[1][0]; RM[2][1] = -RM[2][1]; RM[1][2] = -RM[1][2]; }
    return RM;
  }

  // 2D ellipses from RM (projection = integrated over the other axis). Each value:
  // [FWHM_1, FWHM_2, tilt_deg]. Qplane = (Q∥,Q⊥); QEplane = (Q∥,E).
  function ellipseParams(RM) {
    var F = 1.0 / Math.sqrt(8.0 * Math.log(2.0)), out = {};
    function proj(a, b, ab, ac, bc, cc) {   // project out the 3rd axis (cc) → integrated
      var AP = a - ac * ac / cc, B = b - bc * bc / cc, C = ac * bc / cc - ab;
      var V = 0.5 * Math.atan2(-2 * C, AP - B), c2 = Math.cos(V), s2 = Math.sin(V), sd = Math.sin(2 * V);
      return [1 / Math.sqrt(AP * c2 * c2 + B * s2 * s2 - C * sd) / F,
              1 / Math.sqrt(AP * s2 * s2 + B * c2 * c2 + C * sd) / F, V * 180 / Math.PI];
    }
    function slc(a, b, ab) {                  // slice (3rd axis = 0) → cross-section
      var V = 0.5 * Math.atan2(2 * ab, a - b), c2 = Math.cos(V), s2 = Math.sin(V), sd = Math.sin(2 * V);
      return [1 / Math.sqrt(a * c2 * c2 + b * s2 * s2 + ab * sd) / F,
              1 / Math.sqrt(a * s2 * s2 + b * c2 * c2 - ab * sd) / F, V * 180 / Math.PI];
    }
    // Qplane = (Q∥,Q⊥) integrating/cutting E ; QEplane = (Q∥,E) integrating/cutting Q⊥
    out.Qplane = proj(RM[0][0], RM[1][1], RM[0][1], RM[0][2], RM[1][2], RM[2][2]);
    out.Qplane_slice = slc(RM[0][0], RM[1][1], RM[0][1]);
    out.QEplane = proj(RM[0][0], RM[2][2], RM[0][2], RM[0][1], RM[1][2], RM[1][1]);
    out.QEplane_slice = slc(RM[0][0], RM[2][2], RM[0][2]);
    return out;
  }

  // Build instrument params from cfg + collimations (arcmin) and return the ellipses
  // for |Q| and energy transfer E. collim = {a1,a2,a3,a4} arcmin (ALZ,ALM,ALA,AL3).
  // Fixed GPTAS constants from 4G_O_foc_3foc_slit_*.d. Returns null if impossible.
  function resolution(cfg, Q, E, collim) {
    collim = collim || {};
    var c = function (v, d) { return (v > 0 ? v : d) * DEGRAD; };
    var kfix = +cfg.kfix;
    var p = {
      AKI: (cfg.fixed === "ki" ? 1 : -1) * kfix,
      TM: TWO_PI / (+cfg.d_mono), TA: TWO_PI / (+cfg.d_ana),
      ALZ: c(+collim.a1, 30), ALM: c(+collim.a2, 90), ALA: c(+collim.a3, 100), AL3: c(+collim.a4, 90),
      VBETZ: 10000 * DEGRAD, VBET1: 175 * DEGRAD, VBET2: 275 * DEGRAD, VBET3: 10000 * DEGRAD,
      ETAM: 30 * DEGRAD, ETAA: 30 * DEGRAD, ETAS: 0 * DEGRAD,
      EPM: 1, EP: 1, ERL: 1, IR1: 0
    };
    if (!(Q > 1e-6)) return null;
    var RM = resMatrix(Q, E, p);
    if (!RM) return null;
    var ep = ellipseParams(RM);
    ep.RM = RM;                       // reduced 3×3 (Q∥,Q⊥,E) — for dispersion-axis rotation
    return ep;
  }

  // QE resolution ellipse measured along an in-plane direction that is psi degrees from Q∥
  // (not along Q itself). Needed for the dispersion map: the swept axis is a*/b*, so when Q
  // is NOT along that axis (off-symmetry hk, or non-orthogonal lattice) the relevant width
  // and focusing tilt are along the dispersion axis, obtained by rotating RM in-plane by psi.
  // psi = azimuth(Q) − azimuth(dispersion axis). Returns [FWHM_along, FWHM_E, tilt_deg].
  // mode "slice" gives the Q⊥ = 0 cross-section instead of the (default) energy-integrated
  // projection — same rotation, then the raw 2×2 submatrix rather than the Schur complement.
  function qeAlong(RM, psiDeg, mode) {
    if (!RM) return null;
    var p = psiDeg * Math.PI / 180, c = Math.cos(p), s = Math.sin(p);
    var R = [[c, -s, 0], [s, c, 0], [0, 0, 1]], Rt = [[c, s, 0], [-s, c, 0], [0, 0, 1]];
    function mm(A, B) {
      var O = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
      for (var i = 0; i < 3; i++) for (var j = 0; j < 3; j++) {
        var t = 0; for (var k = 0; k < 3; k++) t += A[i][k] * B[k][j]; O[i][j] = t;
      }
      return O;
    }
    var M = mm(mm(R, RM), Rt);                 // M' = R · RM · Rᵀ (E axis fixed)
    var F = 1.0 / Math.sqrt(8.0 * Math.log(2.0));
    if (mode === "slice") {                    // (axis0, E) cross-section at transverse Q⊥ = 0
      var sa = M[0][0], sb = M[2][2], sab = M[0][2];
      var SV = 0.5 * Math.atan2(2 * sab, sa - sb), sc = Math.cos(SV), ss = Math.sin(SV), ssd = Math.sin(2 * SV);
      return [1 / Math.sqrt(sa * sc * sc + sb * ss * ss + sab * ssd) / F,
              1 / Math.sqrt(sa * ss * ss + sb * sc * sc - sab * ssd) / F, SV * 180 / Math.PI];
    }
    // default: project out the transverse axis (1) → (axis0, E) energy-integrated ellipse
    var a = M[0][0], b = M[2][2], ab = M[0][2], ac = M[0][1], bc = M[1][2], cc = M[1][1];
    var AP = a - ac * ac / cc, B = b - bc * bc / cc, C = ac * bc / cc - ab;
    var V = 0.5 * Math.atan2(-2 * C, AP - B), c2 = Math.cos(V), s2 = Math.sin(V), sd = Math.sin(2 * V);
    return [1 / Math.sqrt(AP * c2 * c2 + B * s2 * s2 - C * sd) / F,
            1 / Math.sqrt(AP * s2 * s2 + B * c2 * c2 + C * sd) / F, V * 180 / Math.PI];
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
  function headerComment(cfg, prefix, spec) {
    var s = '" ' + (prefix || "") + "a=" + cfg.a + " b=" + cfg.b + " c=" + cfg.c
      + " fixed-" + (cfg.fixed || "kf") + "=" + cfg.kfix
      + " ana=" + anaName(+cfg.d_ana);
    // filman's lattice constants are MOUNT-dependent: AS=|Q(u)|, BS=|Q(v)|, CG=cos∠(Q(u),Q(v))
    // (= a*, b*, cosγ* only for the default (100)/(010) plane). Print them so the operator can
    // copy them straight into the filman file for this scattering plane.
    if (spec) {
      var u = cfg.plane_u || [1, 0, 0], v = cfg.plane_v || [0, 1, 0];
      var Bu = matVec(spec.B, [+u[0], +u[1], +u[2]]), Bv = matVec(spec.B, [+v[0], +v[1], +v[2]]);
      var AS = norm(Bu), BS = norm(Bv);
      if (AS > 1e-9 && BS > 1e-9) {
        var CGv = dot(Bu, Bv) / (AS * BS);
        if (Math.abs(CGv) < 5e-5) CGv = 0;               // hide float dust (-0.0000)
        var f = function (w) { return w.map(function (x) { return gfmt(+x); }).join(" "); };
        s += '\n" set filman for plane u=(' + f(u) + ") v=(" + f(v) + "): AS=" + AS.toFixed(4)
          + " BS=" + BS.toFixed(4) + " CG=" + CGv.toFixed(4);
      }
    }
    return s;
  }
  // Express a Miller vector in the scattering-plane basis (u,v): P ≈ ξu·u + ξv·v.
  // filman specifies scans in these in-plane coordinates (HS along u, KS along v), NOT raw
  // Miller h,k. 2×2 Gram solve; exact for in-plane P (out-of-plane parts are projected out).
  // Default u=(100),v=(010) → ξu=h, ξv=k, so the usual h×k output is unchanged.
  function planeCoords(u, v, P) {
    var uu = dot(u, u), vv = dot(v, v), uv = dot(u, v);
    var det = uu * vv - uv * uv;
    if (!(Math.abs(det) > 1e-12)) return [+P[0], +P[1]];     // degenerate basis → raw h,k
    var up = u[0]*P[0] + u[1]*P[1] + u[2]*P[2];
    var vp = v[0]*P[0] + v[1]*P[1] + v[2]*P[2];
    return [(up*vv - vp*uv)/det, (vp*uu - up*uv)/det];
  }
  function scanLine(no, start, step, npts, monitor, nt, u, v) {
    u = u || [1, 0, 0]; v = v || [0, 1, 0];
    var s = planeCoords(u, v, start), d = planeCoords(u, v, step);
    var f = ["NS=" + (no | 0), "HS=" + gfmt(s[0]), "KS=" + gfmt(s[1]), "ES=" + gfmt(+start[3]),
             "DE=" + gfmt(+step[3]), "DH=" + gfmt(d[0]), "DK=" + gfmt(d[1]),
             "NP=" + (npts | 0), "MN=" + (monitor | 0)];
    if ((nt | 0) > 0) f.push("NT=" + (nt | 0));
    return f.join(",");
  }
  function to_scn(cfg, scan, scan_no) {
    scan_no = scan_no || 1;
    var u = cfg.plane_u || [1, 0, 0], v = cfg.plane_v || [0, 1, 0];
    var spec = build(cfg).spec;
    function oopc(P) { var Q = matVec(spec.B, [+P[0], +P[1], +P[2]]); return Math.abs(dot(Q, spec.n)); }
    var L = [headerComment(cfg, "", spec),
      scanLine(scan_no, scan.start, scan.step, scan.npts,
        scan.monitor == null ? 10000 : scan.monitor, scan.nt || 0, u, v),
      "GO " + scan_no];
    if (oopc(scan.start) > 1e-4 || oopc(scan.step) > 1e-4)
      L.push('" NOTE: scan leaves the scattering plane; its out-of-plane part is not in the .scn and is unreachable');
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
    var u = cfg.plane_u || [1, 0, 0], v = cfg.plane_v || [0, 1, 0];
    var scans = mapScans(cfg, md, trim), n = scans.length;
    var spec = build(cfg).spec;          // reused for the header AS/BS/CG + the out-of-plane note
    var L = [headerComment(cfg, "map " + n + " lines: ", spec)];
    if (n === 0) { L.push('" no reachable points; revise the range'); return L.join("\n"); }
    var nblocks = Math.floor((n + CHUNK - 1) / CHUNK);
    if (nblocks > 1)
      L.push('" ' + n + " lines > " + CHUNK + "/GO: emitted as " + nblocks
        + " GO blocks (NS restarts at 1 in each block)");
    for (var c0 = 0; c0 < n; c0 += CHUNK) {
      var chunk = scans.slice(c0, c0 + CHUNK);
      for (var i = 0; i < chunk.length; i++)
        L.push(scanLine(i + 1, chunk[i].start, chunk[i].step, chunk[i].npts,
          chunk[i].monitor, chunk[i].nt, u, v));
      var k = chunk.length;
      L.push(k === 1 ? "GO 1" : "GO 1-" + k);
    }
    function oopc(hkl) { var Q = matVec(spec.B, [+hkl[0], +hkl[1], +hkl[2]]); return Math.abs(dot(Q, spec.n)); }   // vs the actual plane normal
    if (oopc(md.start) > 1e-4 || oopc(md.step) > 1e-4 || oopc(md.outer) > 1e-4)
      L.push('" NOTE: scan leaves the scattering plane (out-of-plane points are unreachable)');
    return L.join("\n");
  }

  var API = {
    DEFAULT_CONFIG: DEFAULT_CONFIG, DEFAULT_LIMITS: DEFAULT_LIMITS, CHUNK: CHUNK,
    Unreachable: Unreachable,
    build: build, anglesRaw: function (cfg, hkl, e) { return angles(buildSpec(cfg), hkl, e); },
    check_point: function (cfg, hkl, e) { var b = build(cfg); return checkPoint(b.spec, b.limits, hkl, e); },
    evaluate: evaluate, grid: grid, gridQ: gridQ, reflections: reflections,
    satellites: satellites, brillouinZone: brillouinZone, resolution: resolution, qeAlong: qeAlong,
    qxy: function (cfg, hkl) { var b = build(cfg); return qd(b.spec, hkl); },
    to_scn: to_scn, evaluate_map: evaluate_map, map_to_scn: map_to_scn
  };

  if (typeof module !== "undefined" && module.exports) module.exports = API;
  else root.Planner = API;
})(typeof self !== "undefined" ? self : this);
