# filman scan planner — standalone (no server)

This folder is a **self-contained, backend-free** build of the planner. All the
crystallography runs in the browser via `planner.js` (a JavaScript port of
`planner.py` + `filman_ng.crystal`, verified numerically identical). It needs
**no Python, no numpy, no server, no install** — just static file hosting.

```
static/
  index.html      ← the app (open this)
  planner.js      ← the computation (ported from Python)
  plotly.min.js   ← plotting library (bundled, no CDN)
```

## Try it locally
- Easiest: double-click `index.html` (opens in your browser).
- Or serve the folder: `python -m http.server` inside `static/`, then open
  `http://127.0.0.1:8000/`.

## Publish it for others (pick one)

**A. GitHub Pages (free, no IT request, gives a public URL)**
1. Create a GitHub repo and upload the **contents of this `static/` folder**
   (index.html, planner.js, plotly.min.js) to the repo root (or a `/docs` folder).
2. Repo → Settings → Pages → Source = your branch, folder = root (or `/docs`).
3. After a minute you get `https://<user>.github.io/<repo>/`. Link it from the
   GPTAS HP. Done — nothing runs server-side.

**B. Netlify / Cloudflare Pages (drag-and-drop)**
Drag the `static/` folder onto app.netlify.com/drop → instant public URL.

**C. On octa (or any campus web server)**
Ask whoever administers the site to drop these three files into a directory the
web server already serves (e.g. a `scan-planner/` folder under the GPTAS site).
No Python/proxy/daemon needed — they are plain static files. (This sidesteps the
numpy-crash and no-root issues entirely.)

## Notes
- Computation is identical to the Python/server version (cross-checked: angles to
  1e-9, reachability, NP/line counts, `.scn` text). If results ever look off,
  compare against the server build.
- The app is read-only and stores nothing; safe to expose publicly.
- To update after changing `planner/planner.py`, re-port the change into
  `planner.js` and re-run the cross-check before publishing.
