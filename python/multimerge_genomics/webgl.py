"""Optional 3-D views, drawn with WebGL only when asked for.

Each function returns a self-contained HTML page (no library, nothing
fetched from the web) that the interface shows in an iframe: a point cloud
(a PCA or a PCoA with three components) or a graph (the assembly graph laid
out in three dimensions). Drag to turn, wheel to fly towards the cursor (into
the cloud, to see its inner boundaries), right-drag to pan, hover a point for
its name; the legend toggles a group; a "Links" button shows or hides the
lines between points (the graph's edges, or the links a point cloud is given).

    html = scatter3d([(x, y, z, "EUR", "HG001"), ...], ("PC1", "PC2", "PC3"), "People on the panel's groups")
    html = graph3d(nodes, edges, "Assembly graph")      # nodes: (x, y, z, size, kind, label); edges: (i, j)
    xyz = layout3d(n_nodes, edges)                       # a force-directed layout, for graph3d
"""
from __future__ import annotations

import html as _html
import json
import math
from typing import Sequence

import numpy as np

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]

_PAGE = r"""<!doctype html>
<html__THEME__><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
:root { color-scheme: light dark; --bg: #f7f7f9; --fg: #1d2230; --muted: #6b7384; --line: #d8dbe2; --panel: rgba(255,255,255,.85); }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --bg: #14171f; --fg: #e7eaf0; --muted: #9aa3b5; --line: #343a46; --panel: rgba(20,23,31,.85); } }
:root[data-theme="dark"] { color-scheme: dark; --bg: #14171f; --fg: #e7eaf0; --muted: #9aa3b5; --line: #343a46; --panel: rgba(20,23,31,.85); }
:root[data-theme="light"] { color-scheme: light; }
html, body { margin: 0; height: 100%; background: var(--bg); color: var(--fg); font: 13px/1.4 "Segoe UI", system-ui, sans-serif; overflow: hidden; }
#wrap { position: relative; width: 100%; height: 100%; }
canvas { display: block; width: 100%; height: 100%; cursor: grab; touch-action: none; }
canvas:active { cursor: grabbing; }
#top { position: absolute; left: 10px; top: 8px; right: 10px; display: flex; gap: 8px; align-items: center; flex-wrap: wrap; pointer-events: none; }
#top b { font-size: 14px; }
#top .hint { color: var(--muted); font-size: 12px; }
#top button { pointer-events: auto; font: 12px "Segoe UI", system-ui, sans-serif; padding: 3px 9px; border: 1px solid var(--line); border-radius: 6px; background: var(--panel); color: var(--fg); cursor: pointer; }
#legend { position: absolute; right: 10px; bottom: 10px; background: var(--panel); border: 1px solid var(--line); border-radius: 8px; padding: 6px 10px; max-height: 45%; overflow: auto; }
#legend div { display: flex; align-items: center; gap: 7px; cursor: pointer; padding: 2px 0; user-select: none; }
#legend div.off { opacity: .4; }
#legend .sw { width: 11px; height: 11px; border-radius: 3px; flex: none; }
#axes { position: absolute; left: 10px; bottom: 10px; max-width: 50%; color: var(--muted); font-size: 12px; background: var(--panel); border: 1px solid var(--line); border-radius: 8px; padding: 4px 8px; }
#tip { position: absolute; display: none; pointer-events: none; background: var(--panel); border: 1px solid var(--line); border-radius: 6px; padding: 4px 8px; font-size: 12px; max-width: 420px; word-break: break-all; }
#nogl { display: none; padding: 20px; }
</style></head>
<body><div id="wrap">
<canvas id="c"></canvas>
<div id="top"><b>__TITLE__</b><span class="hint">drag: turn · wheel: fly towards the cursor, into the cloud · right-drag: pan · hover: name</span>
  <button id="reset">⟲ Reset</button><button id="spin">⟳ Spin</button><button id="links" style="display:none">⟋ Links</button></div>
<div id="legend"></div><div id="axes">__AXES__</div><div id="tip"></div>
<div id="nogl">WebGL is not available in this browser, so the 3-D view cannot be drawn. The 2-D figure above has the same data.</div>
</div>
<script>
const DATA = __DATA__;
(() => {
  const canvas = document.getElementById('c');
  const gl = canvas.getContext('webgl', {antialias: true, alpha: false});
  if (!gl) { document.getElementById('nogl').style.display = 'block'; canvas.style.display = 'none'; return; }
  const forced = document.documentElement.dataset.theme;
  const dark = forced ? forced === 'dark' : matchMedia('(prefers-color-scheme: dark)').matches;
  const P = DATA.points, G = DATA.groups, n = P.length;
  // links between points: the graph's edges, or (a point cloud) whatever joins consecutive points;
  // more than 65,535 points need 32-bit indices, which is an extension of WebGL 1
  const big = n > 65535, ext = big ? gl.getExtension('OES_element_index_uint') : true;
  const E = ext ? (DATA.edges || []) : [];
  const vis = G.map(() => true);
  let showLinks = DATA.links_on !== false;
  // buffers
  const pos = new Float32Array(n * 3), col = new Float32Array(n * 3), siz = new Float32Array(n);
  const hex = (h) => [parseInt(h.slice(1, 3), 16) / 255, parseInt(h.slice(3, 5), 16) / 255, parseInt(h.slice(5, 7), 16) / 255];
  const gcol = G.map((g) => hex(g[1]));
  for (let i = 0; i < n; i++) { pos[3 * i] = P[i][0]; pos[3 * i + 1] = P[i][1]; pos[3 * i + 2] = P[i][2]; }
  const fillColors = () => { for (let i = 0; i < n; i++) { const c = gcol[P[i][3]] || [.5, .5, .5]; const on = vis[P[i][3]];
    col[3 * i] = c[0]; col[3 * i + 1] = c[1]; col[3 * i + 2] = c[2]; siz[i] = on ? (P[i][5] || DATA.size || 6) : 0; } };
  fillColors();
  // the links drawn: only those whose two ends are in groups still shown
  let eidx = new Uint16Array(0), nIdx = 0;
  const fillEdges = () => { const keep = E.filter((e) => vis[P[e[0]][3]] && vis[P[e[1]][3]]);
    eidx = new (big ? Uint32Array : Uint16Array)(keep.length * 2); keep.forEach((e, k) => { eidx[2 * k] = e[0]; eidx[2 * k + 1] = e[1]; }); nIdx = eidx.length; };
  fillEdges();
  const mk = (type, src) => { const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s); return s; };
  const prog = gl.createProgram();
  // a point grows as the camera comes closer, but gently (a 0.7 power) and at most 4x,
  // so that a zoomed-in cloud stays a cloud of points, not a wall of discs
  gl.attachShader(prog, mk(gl.VERTEX_SHADER, `attribute vec3 p; attribute vec3 c; attribute float s; uniform mat4 m; uniform float k; varying vec3 vc; varying float vs;
    void main() { gl_Position = m * vec4(p, 1.0); float w = max(gl_Position.w, 0.05); gl_PointSize = s * k * clamp(pow(3.2 / w, 0.7), 0.35, 4.0); vc = c; vs = s; }`));
  gl.attachShader(prog, mk(gl.FRAGMENT_SHADER, `precision mediump float; varying vec3 vc; varying float vs; uniform float line; uniform float la;
    void main() { if (line < 0.5) { vec2 d = gl_PointCoord - vec2(0.5); if (dot(d, d) > 0.25 || vs <= 0.0) discard; }
    gl_FragColor = vec4(vc, line > 0.5 ? la : 1.0); }`));
  gl.linkProgram(prog); gl.useProgram(prog);
  const A = (name, data, size) => { const b = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, b); gl.bufferData(gl.ARRAY_BUFFER, data, gl.DYNAMIC_DRAW);
    const loc = gl.getAttribLocation(prog, name); gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, size, gl.FLOAT, false, 0, 0); return b; };
  const bp = A('p', pos, 3), bc = A('c', col, 3), bs = A('s', siz, 1);
  const be = gl.createBuffer(); gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, be); gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, eidx, gl.STATIC_DRAW);
  const um = gl.getUniformLocation(prog, 'm'), uk = gl.getUniformLocation(prog, 'k'), ul = gl.getUniformLocation(prog, 'line'), ula = gl.getUniformLocation(prog, 'la');
  // links: thin white lines on the dark theme (dark grey on the light one), translucent,
  // fainter in a point cloud (where they are a help to the eye, not the data) than in a graph
  const edgeCol = new Float32Array(n * 3).fill(dark ? 0.92 : 0.3);
  const edgeAlpha = DATA.kind === 'graph' ? 0.5 : 0.3;
  // camera: the cloud fits a cube of side 2. The camera looks at a point `cen` (world
  // coordinates) from a distance `dist`; drag turns around that point. The wheel flies
  // towards whatever is under the cursor: first by coming closer (dist shrinks), then,
  // once dist is at its floor, by carrying the point looked at forward along the ray —
  // so the camera enters the cloud, and a drag inside it looks around from where it is.
  const F = 1 / Math.tan(0.4), NEAR = 0.01, DMIN = 0.12;
  let yaw = 0.6, pitch = 0.35, dist = 3.2, cen = [0, 0, 0], spin = false, R = null;
  const reset = () => { yaw = 0.6; pitch = 0.35; dist = 3.2; cen = [0, 0, 0]; };
  window.cameraState = () => ({yaw, pitch, dist, centre: cen.slice()});   // for the tests
  const rot = () => { const cy = Math.cos(yaw), sy = Math.sin(yaw), cp = Math.cos(pitch), sp = Math.sin(pitch);
    return [cy, 0, -sy, 0,  sy * sp, cp, cy * sp, 0,  sy * cp, -sp, cy * cp, 0,  0, 0, 0, 1]; };   // yaw about y, pitch about x; column-major
  const toWorld = (v) => [R[0] * v[0] + R[1] * v[1] + R[2] * v[2], R[4] * v[0] + R[5] * v[1] + R[6] * v[2], R[8] * v[0] + R[9] * v[1] + R[10] * v[2]];   // R is orthonormal: its transpose undoes it
  const mat = () => {
    const w = canvas.width, h = canvas.height, asp = w / h, f = F, near = NEAR, far = 50;
    R = rot();
    const Tc = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, -cen[0], -cen[1], -cen[2], 1];
    const T = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, -dist, 1];
    const Pm = [f / asp, 0, 0, 0, 0, f, 0, 0, 0, 0, (far + near) / (near - far), -1, 0, 0, 2 * far * near / (near - far), 0];
    const mul = (a, b) => { const o = new Array(16).fill(0); for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) for (let k = 0; k < 4; k++) o[j * 4 + i] += a[k * 4 + i] * b[j * 4 + k]; return o; };
    return mul(Pm, mul(T, mul(R, Tc)));
  };
  let M = null;
  const draw = () => {
    const dpr = window.devicePixelRatio || 1;
    const W = Math.floor(canvas.clientWidth * dpr), H = Math.floor(canvas.clientHeight * dpr);
    if (canvas.width !== W || canvas.height !== H) { canvas.width = W; canvas.height = H; gl.viewport(0, 0, W, H); }
    M = mat();
    gl.clearColor(...(dark ? [0.08, 0.09, 0.12] : [0.97, 0.97, 0.98]), 1); gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.enable(gl.BLEND); gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA); gl.enable(gl.DEPTH_TEST);
    gl.uniformMatrix4fv(um, false, new Float32Array(M)); gl.uniform1f(uk, dpr); gl.uniform1f(ula, edgeAlpha);
    if (nIdx && showLinks) { gl.uniform1f(ul, 1); gl.bindBuffer(gl.ARRAY_BUFFER, bc); gl.bufferData(gl.ARRAY_BUFFER, edgeCol, gl.DYNAMIC_DRAW);
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, be); gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, eidx, gl.DYNAMIC_DRAW);
      gl.drawElements(gl.LINES, nIdx, big ? gl.UNSIGNED_INT : gl.UNSIGNED_SHORT, 0); }
    gl.uniform1f(ul, 0); gl.bindBuffer(gl.ARRAY_BUFFER, bc); gl.bufferData(gl.ARRAY_BUFFER, col, gl.DYNAMIC_DRAW);
    gl.bindBuffer(gl.ARRAY_BUFFER, bs); gl.bufferData(gl.ARRAY_BUFFER, siz, gl.DYNAMIC_DRAW);
    gl.drawArrays(gl.POINTS, 0, n);
  };
  const loop = () => { if (spin) yaw += 0.006; draw(); requestAnimationFrame(loop); };
  // controls
  let drag = null;
  canvas.addEventListener('pointerdown', (e) => { drag = {x: e.clientX, y: e.clientY, b: e.button, shift: e.shiftKey}; canvas.setPointerCapture(e.pointerId); });
  canvas.addEventListener('pointerup', () => { drag = null; });
  canvas.addEventListener('contextmenu', (e) => e.preventDefault());
  canvas.addEventListener('pointermove', (e) => {
    if (drag) { const dx = e.clientX - drag.x, dy = e.clientY - drag.y; drag.x = e.clientX; drag.y = e.clientY;
      if (drag.b === 2 || drag.shift) { if (!R) return; const s = toWorld([-dx * dist / canvas.clientWidth * 1.6, dy * dist / canvas.clientHeight * 1.6, 0]);
        cen = [cen[0] + s[0], cen[1] + s[1], cen[2] + s[2]]; }
      else { yaw += dx * 0.008; pitch = Math.max(-1.5, Math.min(1.5, pitch + dy * 0.008)); }
      return; }
    hover(e);
  });
  // the wheel flies towards the cursor: the point under it stays under it while the camera
  // comes closer; at the floor of dist the point looked at itself moves along the ray, so the
  // camera goes on into the cloud (wheel back: the camera pulls away from the point looked at; Reset: the start)
  canvas.addEventListener('wheel', (e) => { e.preventDefault(); if (!R) return;
    const r = canvas.getBoundingClientRect(), asp = r.width / r.height;
    const nx = ((e.clientX - r.left) / r.width) * 2 - 1, ny = 1 - ((e.clientY - r.top) / r.height) * 2;
    const X = nx * asp / F * dist, Y = ny / F * dist;           // the point under the cursor, at the depth of the point looked at (camera space)
    const ratio = e.deltaY > 0 ? 1.12 : 1 / 1.12, want = dist * ratio;
    let s;
    if (want >= DMIN && want <= 20) { s = toWorld([(1 - ratio) * X, (1 - ratio) * Y, 0]); dist = want; }       // come closer: the point stays under the cursor
    else if (want < DMIN) { const t = (1 - ratio) * Math.max(dist, 0.3) / dist;   // fly: carry the point looked at forward along the ray to the cursor (a step of ~3% of the cloud)
      s = toWorld([t * X, t * Y, -t * dist]); dist = Math.max(dist, DMIN); }
    else return;
    cen = [cen[0] + s[0], cen[1] + s[1], cen[2] + s[2]]; }, {passive: false});
  document.getElementById('reset').onclick = reset;
  document.getElementById('spin').onclick = () => { spin = !spin; };
  const lbtn = document.getElementById('links');
  if (E.length) { lbtn.style.display = ''; const paint = () => { lbtn.style.opacity = showLinks ? '1' : '.5'; lbtn.title = (DATA.links_label || 'the links between the points') + (showLinks ? ' (shown)' : ' (hidden)'); };
    paint(); lbtn.onclick = () => { showLinks = !showLinks; paint(); }; }
  const tip = document.getElementById('tip');
  const hover = (e) => {
    if (!M) return;
    const r = canvas.getBoundingClientRect(), mx = e.clientX - r.left, my = e.clientY - r.top;
    let best = -1, bd = 10 * 10;
    for (let i = 0; i < n; i++) {
      if (siz[i] <= 0) continue;
      const x = pos[3 * i], y = pos[3 * i + 1], z = pos[3 * i + 2];
      const cx = M[0] * x + M[4] * y + M[8] * z + M[12], cy = M[1] * x + M[5] * y + M[9] * z + M[13], cw = M[3] * x + M[7] * y + M[11] * z + M[15];
      if (cw <= 0) continue;
      const sx = (cx / cw * 0.5 + 0.5) * r.width, sy = (0.5 - cy / cw * 0.5) * r.height;
      const d = (sx - mx) * (sx - mx) + (sy - my) * (sy - my);
      if (d < bd) { bd = d; best = i; }
    }
    if (best < 0) { tip.style.display = 'none'; return; }
    const p = P[best];
    tip.textContent = (p[4] ? p[4] + ' · ' : '') + (G[p[3]] ? G[p[3]][0] : '');
    tip.style.display = 'block'; tip.style.left = Math.min(mx + 12, r.width - 200) + 'px'; tip.style.top = (my + 12) + 'px';
  };
  // legend
  const leg = document.getElementById('legend');
  G.forEach((g, i) => { const d = document.createElement('div'); d.innerHTML = '<span class="sw" style="background:' + g[1] + '"></span>' + g[0] + (g[2] ? ' <span style="opacity:.6">(' + g[2] + ')</span>' : '');
    d.onclick = () => { vis[i] = !vis[i]; d.classList.toggle('off', !vis[i]); fillColors(); fillEdges(); }; leg.appendChild(d); });
  if (!G.length) leg.style.display = 'none';
  loop();
})();
</script></body></html>
"""


def _normalize(xyz: np.ndarray) -> np.ndarray:
    """Centred, scaled to fit a cube of side 2 (the camera is set for that)."""
    if len(xyz) == 0:
        return xyz
    c = xyz - xyz.mean(axis=0)
    span = float(np.abs(c).max()) or 1.0
    return c / span


def _page(title: str, axes: str, data: dict, theme: str = "auto") -> str:
    payload = json.dumps(data, separators=(",", ":"), allow_nan=False)
    attr = f' data-theme="{theme}"' if theme in ("dark", "light") else ""
    return (_PAGE.replace("__THEME__", attr).replace("__TITLE__", _html.escape(title)).replace("__AXES__", _html.escape(axes))
            .replace("__DATA__", payload.replace("</", "<\\/")))


def scatter3d(points: Sequence[tuple], axes: Sequence[str] = ("x", "y", "z"), title: str = "3-D view",
              size: float = 7.0, group_notes: dict | None = None, theme: str = "auto",
              links: Sequence[tuple[int, int]] | None = None, links_label: str | None = None, links_on: bool = True) -> str:
    """A point cloud. `points`: (x, y, z, group, label) tuples; a group of
    None or "" is "other". `group_notes`: a short text per group for the legend.
    `theme`: "dark", "light" or "auto" (the viewer's own preference).
    `links`: pairs of point indices joined by a thin line (consecutive windows
    of one contig, say); a "Links" button shows and hides them (`links_on`
    says how they start) and `links_label` names them in its tooltip."""
    pts = [tuple(p) for p in points]
    edges = [[int(a), int(b)] for a, b in (links or ()) if 0 <= int(a) < len(pts) and 0 <= int(b) < len(pts) and int(a) != int(b)]
    groups: list[str] = []
    for p in pts:
        g = str(p[3]) if len(p) > 3 and p[3] not in (None, "") else "other"
        if g not in groups:
            groups.append(g)
    xyz = _normalize(np.array([[float(p[0]), float(p[1]), float(p[2]) if len(p) > 2 and p[2] is not None else 0.0] for p in pts], dtype=float))
    data = {
        "points": [[round(float(x), 4), round(float(y), 4), round(float(z), 4),
                    groups.index(str(p[3]) if len(p) > 3 and p[3] not in (None, "") else "other"),
                    str(p[4]) if len(p) > 4 and p[4] is not None else ""] for (x, y, z), p in zip(xyz, pts)],
        "groups": [[g, SERIES[i % len(SERIES)], (group_notes or {}).get(g, "")] for i, g in enumerate(groups)],
        "edges": edges,
        "size": size,
        "kind": "scatter",
        "links_label": links_label or "",
        "links_on": bool(links_on),
    }
    return _page(title, " · ".join(f"{n}: {a}" for n, a in zip(("x", "y", "z"), axes)), data, theme)


def _unfold(init: np.ndarray, n: int, rng) -> tuple[np.ndarray, np.ndarray]:
    """From a 2-D layout (centred, scaled) to a start in three dimensions.
    A graph laid out as one long line (a genome in one piece: a chain of
    unitigs with bubbles along it) would stay a stick in 3-D, so a layout
    much longer than wide is coiled into a helix: the long axis runs around
    and up the coil, the short one (the bubbles) points out of it. Anything
    else keeps its 2-D shape and gets a third coordinate as wide as its
    shorter side. The coil is presentation only: it has no meaning."""
    cov = init.T @ init / max(n - 1, 1)
    w, v = np.linalg.eigh(cov)                                   # ascending; v[:, 1] is the long axis
    axes = init @ v[:, ::-1]                                     # column 0 along the long axis
    sx = float(axes[:, 0].max() - axes[:, 0].min()) or 1.0
    sy = float(axes[:, 1].max() - axes[:, 1].min())
    if sy <= 1e-6 or sx / sy >= 3.0:
        turns = min(6.0, max(1.5, math.sqrt(n) / 12))            # ~1.5 turns for 300 nodes, ~3.7 for 2,000
        radius = sx / (2 * math.pi * turns)
        u = (axes[:, 0] - axes[:, 0].min()) / sx * 2 * math.pi * turns
        side = float(np.abs(axes[:, 1]).max()) or 1.0
        r = radius + axes[:, 1] * min(1.0, 0.3 * radius / side)   # the bubbles stick out of the coil, never through its axis
        xy = np.stack([r * np.cos(u), r * np.sin(u)], axis=1)
        z = (u / (2 * math.pi * turns) - 0.5) * 2 * radius + rng.normal(size=n) * 0.02 * radius
        return xy.astype(np.float32), z.astype(np.float32)
    z = rng.normal(size=n) * 0.12 * sy
    return init.astype(np.float32), z.astype(np.float32)


def layout3d(n: int, edges: Sequence[tuple[int, int]], iterations: int = 20, seed: int = 1, max_nodes: int = 2500,
             init: np.ndarray | None = None) -> np.ndarray:
    """A force-directed layout in three dimensions (Fruchterman-Reingold, all
    pairs, float32: a few seconds for a couple of thousand nodes). `init`, an
    (n, 2) or (n, 3) array, starts from a layout already made (the 2-D one of
    the engine, say): the third coordinate then only has to unfold. Returns an
    (n, 3) array scaled to the unit cube."""
    rng = np.random.default_rng(seed)
    if n == 0:
        return np.zeros((0, 3))
    pos = rng.normal(size=(n, 3)).astype(np.float32)
    if init is not None and len(init) == n:
        init = np.asarray(init, dtype=np.float32)
        init = init - init.mean(axis=0)
        init = init / (float(np.abs(init).max()) or 1.0) * math.sqrt(n) * 0.5
        if init.shape[1] == 2 and n >= 3:
            pos[:, :2], pos[:, 2] = _unfold(init, n, rng)
        else:
            pos[:, :init.shape[1]] = init[:, :3]
            pos[:, init.shape[1]:] *= 0.3
    if n > max_nodes or n < 2:
        return _normalize(pos)
    k = np.float32(1.0 / math.sqrt(n))
    ei = np.array([(a, b) for a, b in edges if 0 <= a < n and 0 <= b < n], dtype=int).reshape(-1, 2)
    t = 0.3
    for it in range(iterations):
        delta = pos[:, None, :] - pos[None, :, :]                      # (n, n, 3)
        d2 = (delta * delta).sum(axis=2)
        np.fill_diagonal(d2, np.inf)
        disp = (delta * (k * k / d2)[:, :, None]).sum(axis=1)          # repulsion
        if len(ei):
            dv = pos[ei[:, 0]] - pos[ei[:, 1]]
            dn = np.sqrt((dv * dv).sum(axis=1)) + 1e-9
            att = dv * (dn / k)[:, None]                                # attraction along edges
            np.add.at(disp, ei[:, 0], -att)
            np.add.at(disp, ei[:, 1], att)
        ln = np.sqrt((disp * disp).sum(axis=1)) + 1e-9
        pos += (disp / ln[:, None] * np.minimum(ln, t)[:, None]).astype(np.float32)
        t *= 0.96
    return _normalize(pos.astype(float))


def graph3d(nodes: Sequence[tuple], edges: Sequence[tuple[int, int]], title: str = "Graph", kinds: dict | None = None,
            theme: str = "auto") -> str:
    """A graph: `nodes` are (x, y, z, size, kind, label); `edges` pairs of node
    indices. `kinds`: {kind: (legend name, colour)}; a kind absent from it
    gets a series colour."""
    kinds = dict(kinds or {})
    groups: list[str] = []
    for nd in nodes:
        k = str(nd[4]) if len(nd) > 4 and nd[4] is not None else ""
        if k not in groups:
            groups.append(k)
    xyz = _normalize(np.array([[float(nd[0]), float(nd[1]), float(nd[2])] for nd in nodes], dtype=float)) if nodes else np.zeros((0, 3))
    data = {
        "points": [[round(float(x), 4), round(float(y), 4), round(float(z), 4),
                    groups.index(str(nd[4]) if len(nd) > 4 and nd[4] is not None else ""),
                    str(nd[5]) if len(nd) > 5 and nd[5] is not None else "", round(float(nd[3]) if len(nd) > 3 else 6.0, 2)]
                   for (x, y, z), nd in zip(xyz, nodes)],
        "groups": [[kinds.get(g, (g or "node", None))[0], kinds.get(g, (None, SERIES[i % len(SERIES)]))[1] or SERIES[i % len(SERIES)], ""]
                   for i, g in enumerate(groups)],
        "edges": [[int(a), int(b)] for a, b in edges if 0 <= a < len(nodes) and 0 <= b < len(nodes)],
        "size": 6,
        "kind": "graph",
        "links_label": "the edges of the graph",
        "links_on": True,
    }
    return _page(title, "force-directed layout: near = joined in the graph; a long chain is coiled to fit the view (the coil itself means nothing)", data, theme)
