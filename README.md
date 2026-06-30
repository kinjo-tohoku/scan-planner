# filman scan planner

A browser-based **scan planner for the GPTAS triple-axis spectrometer** (beamline 4G, JRR-3).
Lay out your scans in (Q, E), see at a glance which points the spectrometer can actually
reach, check the instrumental resolution, and export ready-to-run **filman `.scn`** scripts.

**▶ Live app: https://kinjo-tohoku.github.io/scan-planner/**

Everything runs in your browser — there is **nothing to install** and no server to set up.
Open the page (or a local `index.html`) and start planning.

## What it does

- **Reachability maps** — for your crystal, scattering plane and fixed energy, shade the
  region the spectrometer can reach within its motor limits. Constant-E **Q–Q** maps and
  **Q–E** dispersion maps.
- **Resolution** — Cooper–Nathans resolution ellipses at the scan centre, including the
  focusing tilt on the dispersion view.
- **Sample-environment magnet** — optionally restrict the reachable region to a magnet's
  windows (vertical ±60° or horizontal ±25°), with a small schematic of the geometry.
- **`.scn` export** — turn a planned line or map into a filman scan script, ready to load
  on the instrument.
- **Estimated time** — a rough measurement-time estimate for the reachable points.

GPTAS defaults (PG002 monochromator, kf = 2.662 Å⁻¹ / 14.7 meV, standard collimations, …)
are filled in; change the crystal, plane, energy and limits to match your experiment.

## Run it locally

You don't have to — just use the live link above — but if you want a local copy:

- **Easiest:** download the three files below and double-click `index.html`.
- **Or serve the folder:** run `python -m http.server` in this directory, then open
  `http://127.0.0.1:8000/`.

```
index.html      ← the app
planner.js      ← all the crystallography / kinematics (runs in the browser)
plotly.min.js   ← plotting library (bundled, no CDN)
```

## Host it yourself

It's just three static files, so any static host works:

- **GitHub Pages** — Settings → Pages → deploy from branch, folder = root. You get a
  public URL.
- **Netlify / Cloudflare Pages** — drag-and-drop the folder for an instant URL.
- **Any campus / web server** — drop the three files into a served directory. No Python,
  proxy or daemon required.

## Notes

- The app is **read-only** — it computes locally and stores nothing, so it's safe to share
  publicly.
- `.scn` scans are written in-plane (along the two scattering-plane vectors); points out of
  the plane are treated as unreachable.
- filman accepts at most **40** scan slots per `GO`, so larger maps are split automatically.

---

Made for GPTAS (4G, JRR-3) users.
