// Build the ME 415 PIV-pressure deck in two formats: "edm" (ME 323 EDM deck style) and "lec" (ME 415 lecture style).
const pptxgen = require("pptxgenjs");
const path = require("path");
const fs = require("fs");
const FMT = process.argv[2];
const OUTFILE = process.argv[3];
const FIG = path.join(__dirname, "figs");
const SIZES = JSON.parse(fs.readFileSync(path.join(FIG, "sizes.json")));
const img = (n) => path.join(FIG, n + ".png");

const pres = new pptxgen();
const ALLS = [];
{ const _add = pres.addSlide.bind(pres);
  pres.addSlide = (...a) => { const sl = _add(...a); sl.__orig = sl.addNotes.bind(sl); sl.__notes = null; sl.addNotes = (t) => { sl.__notes = t; return sl; }; ALLS.push(sl); return sl; }; }
const wide = FMT === "edm";
const SHORT_FMT = FMT === "lec10";
if (wide) pres.layout = "LAYOUT_WIDE";
else { pres.defineLayout({ name: "K43", width: 10, height: 7.5 }); pres.layout = "K43"; }
pres.author = "Anany Bajpai, Arnav Oza";
pres.title = "Pressure from PIV: Algorithms and Comparison with Simulations";

// ------------------------------------------------------------------ format definitions
const BODY = "Georgia";
const TITLEF = wide ? "Georgia" : "Big Caslon";
const C = wide
  ? { body: "000000", head: "000000", note: "404040", accent: "1F4E79", accent2: "4F81BD", cap: "000000" }
  : { body: "1A5288", head: "5A554C", note: "5B5750", accent: "0076BA", accent2: "0076BA", cap: "333333" };
const S = wide ? { body: 13, head: 14, note: 12, cap: 12 } : { body: 13, head: 15, note: 11, cap: 11 };
const G = wide ? { X0: 0.6, X1: 12.73, Y0: 1.08, Y1: 7.15, W: 13.333, H: 7.5 }
               : { X0: 0.4, X1: 9.6, Y0: 1.72, Y1: 7.12, W: 10, H: 7.5 };
G.CW = G.X1 - G.X0;

let figNo = 0, eqNo = 0, pageNo = 0;
const lecSlides = []; // for page numbers (filled at end)

function newSlide(title, subtitle, opts = {}) {
  const s = pres.addSlide();
  if (wide) {
    s.background = { color: "FFFFFF" };
    s.addShape(pres.shapes.RECTANGLE, { x: 0, y: 0.8055, w: G.W, h: 0.0278, fill: { color: "000000" }, line: { type: "none" } });
    s.addText(title, { x: 0.944, y: 0.153, w: 11.444, h: 0.528, fontFace: TITLEF, fontSize: 30, color: "000000", align: "center", valign: "top", margin: 0, isTextBox: true });
  } else {
    s.background = { path: path.join(FIG, "gridpaper.jpg") };
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 0.74, y: -0.4, w: 8.52, h: 1.6, rectRadius: 0.28, fill: { color: "9EC2C9" }, line: { type: "none" } });
    s.addText(title, { x: 0.8, y: 0.05, w: 8.4, h: 1.1, fontFace: TITLEF, fontSize: opts.titleSize || 46, color: "000000", align: "center", valign: "middle", margin: 0, isTextBox: true, fit: "shrink" });
    if (subtitle) s.addText(subtitle, { x: 0.6, y: 1.22, w: 8.8, h: 0.42, fontFace: BODY, fontSize: 19, color: "5B5750", align: "center", valign: "middle", margin: 0, isTextBox: true });
    // footer
    s.addShape(pres.shapes.RECTANGLE, { x: 0, y: 7.295, w: 3.15, h: 0.205, fill: { color: "9EC2C9" }, line: { type: "none" } });
    s.addShape(pres.shapes.RECTANGLE, { x: 3.15, y: 7.295, w: 3.7, h: 0.205, fill: { color: "D5D5D5" }, line: { type: "none" } });
    s.addShape(pres.shapes.RECTANGLE, { x: 6.85, y: 7.295, w: 3.15, h: 0.205, fill: { color: "9EC2C9" }, line: { type: "none" } });
    const ft = { y: 7.295, h: 0.205, fontFace: BODY, fontSize: 9, margin: 0, valign: "middle", isTextBox: true };
    s.addText("Anany Bajpai & Arnav Oza", { ...ft, x: 0.3, w: 2.8, color: "6E6E6E" });
    s.addText("ME 415: CFD & HT", { ...ft, x: 3.15, w: 3.7, color: "000000", align: "center" });
    s.addText("Autumn ’26", { ...ft, x: 6.88, w: 1.5, color: "000000" });
    s.addText("Term Project", { ...ft, x: 8.4, w: 1.55, color: "000000", align: "right" });
    lecSlides.push(s);
  }
  return s;
}

// paragraphs: items = [{k:'head'|'b'|'note'|'text'|'num', t: string | [[text, opts],...]}]
function para(s, items, box, o = {}) {
  const sz = o.size || S.body;
  const runs = [];
  items.forEach((it, i) => {
    const segs = typeof it.t === "string" ? [[it.t, {}]] : it.t;
    let base = { fontFace: BODY, color: C.body, fontSize: sz, lineSpacingMultiple: o.ls || 1.08 };
    let pOpts = { align: o.align || (wide ? "justify" : "left") };
    if (it.k === "head") { base = { ...base, fontSize: o.headSize || S.head, color: C.head, underline: wide ? undefined : { style: "sng" } }; pOpts = { align: "left", paraSpaceBefore: i ? (o.gap ?? 7) : 0, paraSpaceAfter: 2 }; }
    else if (it.k === "b") { pOpts.bullet = wide ? { indent: 12 } : { code: "27A1", indent: 16 }; pOpts.paraSpaceBefore = o.bgap ?? 2; }
    else if (it.k === "num") { pOpts.bullet = o.startAt ? { type: "number", indent: 16, numberStartAt: o.startAt } : { type: "number", indent: 16 }; pOpts.paraSpaceBefore = o.bgap ?? 3; }
    else if (it.k === "note") { base = { ...base, fontSize: o.noteSize || S.note, color: C.note, italic: false }; pOpts.paraSpaceBefore = 2; pOpts.align = "left"; }
    else { pOpts.paraSpaceBefore = i ? (o.pgap ?? 4) : 0; }
    if (it.color) base.color = it.color;
    segs.forEach(([txt, so], j) => {
      const r = { text: txt, options: { ...base, ...so } };
      Object.assign(r.options, pOpts);
      if (j === segs.length - 1 && i < items.length - 1) r.options.breakLine = true;
      runs.push(r);
    });
  });
  s.addText(runs, { x: box.x, y: box.y, w: box.w, h: box.h, margin: 0, valign: o.valign || "top", isTextBox: true });
}

function fit(name, w, h) { const [pw, ph] = SIZES[name]; const a = pw / ph; return (w / h > a) ? { w: h * a, h } : { w, h: w / a }; }
function fig(s, name, box, caption, o = {}) {
  const capH = caption ? (o.capH || 0.3) : 0;
  const d = fit(name, box.w, box.h - capH - (caption ? 0.06 : 0));
  const x = box.x + (o.align === "left" ? 0 : (box.w - d.w) / 2);
  const y = box.y + (o.valign === "bottom" ? (box.h - capH - d.h - 0.06) : 0);
  s.addImage({ path: o.file || img(name), x, y, w: d.w, h: d.h, altText: o.alt || caption || name });
  if (wide && o.border !== false) s.addShape(pres.shapes.RECTANGLE, { x, y, w: d.w, h: d.h, fill: { type: "none" }, line: { color: "000000", width: 0.5 } });
  if (caption) {
    figNo += 1;
    s.addText(`Fig. ${figNo} ${caption}`, { x: box.x, y: y + d.h + 0.06, w: box.w, h: capH, fontFace: BODY, fontSize: o.capSize || S.cap, color: C.cap, align: "center", valign: "top", margin: 0, isTextBox: true });
  }
  return { x, y, w: d.w, h: d.h };
}
function eq(s, name, x, y, scale = 0.62, label) {
  const [pw, ph] = SIZES[name];
  const w = pw / 300 * scale, h = ph / 300 * scale;
  let lx = x;
  if (label) {
    eqNo += 1;
    s.addText(`Eq. (${eqNo})${label === true ? "" : ": " + label}`, { x, y, w: 3.5, h: 0.24, fontFace: BODY, fontSize: S.note, color: C.note, margin: 0, isTextBox: true });
    y += 0.25;
  }
  s.addImage({ path: img(name), x: lx, y, w, h, altText: name });
  return { w, h: h + (label ? 0.25 : 0) };
}
function placeholder(s, box, text) {
  s.addShape(pres.shapes.RECTANGLE, { x: box.x, y: box.y, w: box.w, h: box.h, fill: wide ? { color: "F7F7F7" } : { color: "FFFFFF", transparency: 35 }, line: { color: "8C8C8C", width: 1, dashType: "dash" } });
  s.addText(text, { x: box.x + 0.1, y: box.y + 0.1, w: box.w - 0.2, h: box.h - 0.2, fontFace: BODY, fontSize: 11, italic: true, color: "7F7F7F", align: "center", valign: "middle", margin: 0, isTextBox: true });
}
function table(s, rows, o) {
  const hdr = wide ? { fill: { color: "4F81BD" }, color: "FFFFFF", bold: true } : { fill: { color: "C94779" }, color: "FFFFFF", bold: false, fontFace: TITLEF };
  const body = wide ? { color: "000000" } : { color: "3D3D3D" };
  const border = wide ? { type: "solid", pt: 0.75, color: "CBCBCB" } : { type: "solid", pt: 0.75, color: "3B6EA5" };
  const data = rows.map((r, i) => r.map((c) => {
    const cell = typeof c === "object" && c !== null && c.text !== undefined ? c : { text: String(c) };
    const op = { fontFace: BODY, fontSize: o.size || 11, border, valign: "middle", margin: [2, 4, 2, 4], ...(i === 0 ? hdr : body), ...(cell.options || {}) };
    if (i > 0 && wide && i % 2 === 0) op.fill = op.fill || { color: "EEF3FA" };
    return { text: cell.text, options: op };
  }));
  s.addTable(data, { x: o.x, y: o.y, w: o.w, colW: o.colW, rowH: o.rowH || 0.26, autoPage: false });
}
// flow box
function fbox(s, text, x, y, w, h, o = {}) {
  if (wide) {
    s.addShape(pres.shapes.RECTANGLE, { x, y, w, h, fill: { color: o.fill || "4F81BD" }, line: { color: o.fill || "4F81BD" } });
    s.addText(text, { x, y, w, h, fontFace: BODY, fontSize: o.size || 11, color: o.color || "FFFFFF", align: "center", valign: "middle", margin: 3, isTextBox: true });
  } else {
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.08, fill: { color: o.fillL || "FFFFFF", transparency: o.fillL ? 0 : 20 }, line: { color: "000000", width: 1 } });
    s.addText(text, { x, y, w, h, fontFace: BODY, fontSize: o.size || 11, color: "000000", align: "center", valign: "middle", margin: 3, isTextBox: true });
  }
}
function arrowDown(s, x, y, h) {
  if (wide) s.addShape(pres.shapes.LINE, { x, y, w: 0, h, line: { color: "404040", width: 1.25, endArrowType: "triangle" } });
  else s.addShape(pres.shapes.DOWN_ARROW, { x: x - 0.11, y, w: 0.22, h, fill: { color: "A39E96" }, line: { color: "403C38", width: 0.75 } });
}
function arrowRight(s, x, y, w) {
  if (wide) s.addShape(pres.shapes.LINE, { x, y, w, h: 0, line: { color: "404040", width: 1.25, endArrowType: "triangle" } });
  else s.addShape(pres.shapes.RIGHT_ARROW, { x, y: y - 0.12, w, h: 0.24, fill: { color: "A39E96" }, line: { color: "403C38", width: 0.75 } });
}
const R = (t, o = {}) => [t, o];
const B = (t) => [t, { bold: true }];
const I = (t) => [t, { italic: true }];

// ================================================================== SLIDES
const TEAM = "Anany Bajpai 24B2135 · Arnav Oza 24B2217";

const BL = {};
const SHORT = FMT === "lec10";
BL.cover = () => {
if (!wide) {
  const s = pres.addSlide();
  s.background = { path: path.join(FIG, "gridpaper.jpg") };
  s.addShape(pres.shapes.RECTANGLE, { x: 0.96, y: 0.36, w: 1.48, h: 0.58, fill: { type: "none" }, line: { color: "000000", width: 2.5 } });
  s.addText("ME 415", { x: 0.96, y: 0.36, w: 1.48, h: 0.58, fontFace: TITLEF, fontSize: 28, color: "000000", align: "center", valign: "middle", margin: 0, isTextBox: true });
  const blk = "000000", blu = "45A7DE";
  const L = (pairs, size) => pairs.map(([t, c]) => ({ text: t, options: { color: c, fontSize: size, fontFace: TITLEF } }));
  s.addText(L([["P", blk], ["ressure", blu]], 62), { x: 0.96, y: 0.95, w: 6, h: 1.0, margin: 0, valign: "middle", isTextBox: true });
  s.addText(L([["f", blk], ["rom ", blu], ["PIV", blk]], 44), { x: 0.96, y: 1.85, w: 6, h: 0.75, margin: 0, valign: "middle", isTextBox: true });
  s.addText(L([["A", blk], ["lgorithms ", blu], ["&", blk]], 32), { x: 0.96, y: 2.6, w: 6, h: 0.55, margin: 0, valign: "middle", isTextBox: true });
  s.addText(L([["C", blk], ["omparison with ", blu], ["S", blk], ["imulations", blu]], 32), { x: 0.96, y: 3.1, w: 6.5, h: 0.6, margin: 0, valign: "middle", isTextBox: true });
  s.addText([{ text: "Project: ", options: { color: "8A8A8A" } }, { text: "6", options: { color: "000000", underline: { style: "sng" } } }], { x: 6.2, y: 2.55, w: 3.05, h: 0.6, fontFace: BODY, fontSize: 30, align: "right", margin: 0, isTextBox: true });
  s.addText("Term project · Autumn ’26", { x: 6.2, y: 3.12, w: 3.05, h: 0.3, fontFace: BODY, fontSize: 13, color: "8A8A8A", align: "right", margin: 0, isTextBox: true });
  s.addText("By Anany Bajpai (24B2135) & Arnav Oza (24B2217)", { x: 0.5, y: 3.95, w: 9, h: 0.4, fontFace: BODY, fontSize: 19, color: "8A8A8A", align: "center", margin: 0, isTextBox: true });
  s.addText("Instructor: Prof. Anjishnu Choudhury, MED-IITB", { x: 0.5, y: 4.35, w: 9, h: 0.3, fontFace: BODY, fontSize: 13, color: "8A8A8A", align: "center", margin: 0, isTextBox: true });
  const brown = "8E5A0C";
  const outline = SHORT_FMT ? ["Problem & objectives", "Governing equations & algorithms", "Solution strategy & test cases", "Results: grid, method, boundaries, time step", "CFD cylinder wake (OpenFOAM)", "Literature & conclusions"] : ["Introduction, problem & objectives", "Governing equations & algorithms", "Solution strategy & test cases", "Results & discussion", "CFD comparison & literature", "Conclusions"];
  s.addText(outline.map((t, i) => ({ text: t, options: { bullet: { type: "number", indent: 22 }, color: brown, fontSize: 16, fontFace: BODY, paraSpaceBefore: 6, breakLine: i < outline.length - 1 } })), { x: 0.98, y: 4.9, w: 6, h: 2.3, margin: 0, valign: "top", isTextBox: true });
  s.addNotes("Title. Introduce the team and the project: we rebuild the pressure field from PIV-type velocity data, test three algorithms, and compare against exact solutions, a CFD simulation and published studies.");
}
};

BL.intro = () => {
{
  const s = wide ? newSlide("ME 415 Term Project · Project 6") : newSlide("Introduction", "What is pressure from PIV?");
  const intro = [
    { k: "b", t: [R("Particle Image Velocimetry ("), B("PIV"), R(") photographs tracer particles in a laser sheet twice, Δt apart, and turns their displacement into a velocity field u(x, y). It does not measure pressure.")] },
    { k: "b", t: "Pressure can be recovered from the momentum equation, because its gradient depends only on velocity and how velocity changes in time." },
    { k: "b", t: "We wrote a Python code that rebuilds pressure from PIV-like velocity with three algorithms (direct, omni-directional, Poisson)." },
    { k: "b", t: "We tested it on exact Navier–Stokes solutions with a synthetic PIV error model, ran an OpenFOAM cylinder wake (Re = 100) as a realistic test, and checked our trends against published studies." },
  ];
  if (wide) {
    s.addText("Pressure from PIV: Algorithms and Comparison with Simulations", { x: 0, y: 1.32, w: G.W, h: 0.5, fontFace: BODY, fontSize: 24, color: "000000", align: "center", margin: 0, isTextBox: true });
    fig(s, "f_piv_schematic", { x: 0.6, y: 2.2, w: 4.2, h: 3.05 }, "Planar PIV (schematic)", { capSize: 10 });
    fig(s, "f_hero", { x: 4.95, y: 2.2, w: 2.95, h: 3.05 }, "Our reconstruction", { capSize: 10 });
    para(s, intro, { x: 8.15, y: 2.2, w: 4.6, h: 3.1 }, { size: 12.5, bgap: 4 });
    s.addText("Team: " + TEAM, { x: 0, y: 5.75, w: G.W, h: 0.4, fontFace: BODY, fontSize: 18, color: "000000", align: "center", margin: 0, isTextBox: true });
    s.addText("ME 415: Computational Fluid Dynamics & Heat Transfer · Instructor: Prof. Anjishnu Choudhury", { x: 0, y: 6.25, w: G.W, h: 0.3, fontFace: BODY, fontSize: 11, color: "404040", align: "center", margin: 0, isTextBox: true });
  } else {
    fig(s, "f_piv_schematic", { x: 0.4, y: 1.75, w: 5.2, h: 3.0 }, "Planar PIV (schematic)", { capSize: 10 });
    fig(s, "f_hero", { x: 5.85, y: 1.75, w: 3.75, h: 3.0 }, "Our pressure reconstruction", { capSize: 10 });
    para(s, intro, { x: 0.45, y: 4.85, w: 9.1, h: 2.35 }, { size: 12.5, bgap: 3 });
  }
  s.addNotes("PIV gives velocity only. Pressure follows from Navier–Stokes. Our work: a code with three integration algorithms, verified on exact solutions with a realistic PIV error model, plus a CFD wake and a comparison with the literature.");
}
};

BL.problem = () => {
{
  const s = wide ? newSlide("Problem Statement and Objectives") : newSlide("Problem & Objectives", "Why is pressure from PIV hard, and what did we set out to do?");
  // diagram: why the problem is hard
  const dx = wide ? 0.6 : 0.4, dw = wide ? 3.4 : 2.75, bh = wide ? 0.62 : 0.6;
  let y = G.Y0 + (wide ? 0.15 : 0.05);
  const steps = [
    ["PIV velocity u, v\n(with noise σ)", "4F81BD"],
    ["Differentiate: ∂u/∂t, (u·∇)u, ∇²u\nnoise grows ~σ/Δt, σ/h, σ/h²", "C0504D"],
    ["Pressure gradient\n∇p = G", "4F81BD"],
    ["Integrate G\n(path, rays or Poisson)", "4F81BD"],
    ["Pressure p\n(up to a constant)", "1F4E79"],
  ];
  const gap = wide ? 0.32 : 0.3;
  steps.forEach(([t, col], i) => {
    fbox(s, t, dx, y, dw, bh, { fill: col, size: wide ? 11 : 10, fillL: i === 1 ? "F8D7D3" : undefined });
    if (i < steps.length - 1) arrowDown(s, dx + dw / 2, y + bh, gap);
    y += bh + gap;
  });
  figNo += 1;
  s.addText(`Fig. ${figNo} Where the noise is amplified`, { x: dx - 0.1, y: y - gap + 0.04, w: dw + 0.2, h: 0.28, fontFace: BODY, fontSize: S.cap, color: C.cap, align: "center", margin: 0, isTextBox: true });
  const tx = dx + dw + (wide ? 0.45 : 0.3);
  const items = [
    { k: "head", t: "Problem" },
    { k: "b", t: "Pressure sets forces, noise and cavitation, but taps give only wall values and probes disturb the flow. PIV measures the whole velocity field without touching the flow, but not pressure." },
    { k: "b", t: "Pressure must be inferred: differentiate the measured velocity, then integrate. Differentiation amplifies measurement noise (Fig. " + figNo + ")." },
    { k: "head", t: "Why it matters" },
    { k: "b", t: "Full-field pressure from one non-intrusive measurement: loads on wings, blades and bluff bodies, aeroacoustic sources, and flows where probes cannot go (e.g. heart valves)." },
    { k: "head", t: "Objectives" },
    { k: "num", t: "Implement three algorithms: direct, omni-directional and Poisson." },
    { k: "num", t: "Verify them on exact Navier–Stokes solutions (order of accuracy)." },
    { k: "num", t: "Quantify the effect of spacing, noise, Δt, boundary conditions, window size and denoising." },
    { k: "num", t: "Compare with a CFD cylinder wake (Re = 100) and published studies." },
    { k: "head", t: "Expected output" },
    { k: "text", t: "Error trends and guidelines (method, boundary data, spacing, Δt) for planning PIV-pressure experiments." },
  ];
  para(s, items, { x: tx, y: G.Y0, w: G.X1 - tx, h: G.Y1 - G.Y0 }, { size: wide ? 13 : 12, headSize: wide ? 14 : 14, gap: wide ? 8 : 5 });
  s.addNotes("The core difficulty: velocity has to be differentiated, which amplifies noise, before it can be integrated into pressure. All four objectives are met; the CFD comparison is in the results section.");
}
};

BL.equations = () => {
{
  const s = wide ? newSlide("Governing Equations and Assumptions") : newSlide("Governing Equations", "From velocity to pressure");
  const ex = G.X0, sc = wide ? 0.7 : 0.5;
  let y = G.Y0 + 0.02;
  const eqs = [["eq_ns", "Incompressible Navier–Stokes"], ["eq_G", "Pressure gradient from velocity"], ["eq_euler", "Eulerian material derivative"], ["eq_lagr", "Lagrangian material derivative"], ["eq_poisson", "Pressure-Poisson problem"], ["eq_compat", "Solvability (compatibility)"]];
  if (wide) {
    eqs.forEach(([n, l]) => { const d = eq(s, n, ex, y, sc, l); y += d.h + 0.16; });
    const items = [
      { k: "head", t: "From velocity to pressure" },
      { k: "b", t: "Rearranging momentum (Eq. 1) gives the pressure gradient G from velocity alone (Eq. 2). Every term in G is a derivative of measured data." },
      { k: "b", t: "Du/Dt needs three snapshots. Eulerian: at a fixed point (Eq. 3). Lagrangian: following the fluid particle (Eq. 4), which removes the fast change seen at a point as a pattern is carried past." },
      { k: "b", t: "Taking ∇· of Eq. 2 gives a Poisson equation with Neumann edges (Eq. 5). It is solvable only if Eq. 6 holds, and then only up to a constant, so we always compare p − p̄." },
      { k: "head", t: "Assumptions" },
      { k: "b", t: "Incompressible, Newtonian fluid with constant ρ and μ." },
      { k: "b", t: "Two-dimensional planar flow: out-of-plane terms neglected (true for our test cases and the 2-D wake at Re = 100)." },
      { k: "b", t: "Velocity known on a uniform grid at three instants, n − 1, n, n + 1." },
    ];
    para(s, items, { x: 7.2, y: G.Y0, w: G.X1 - 7.2, h: G.Y1 - G.Y0 }, { gap: 9, bgap: 4 });
  } else {
    const col2 = 5.05;
    eqs.forEach(([n, l], i) => {
      const x = i % 2 === 0 ? ex : col2; const yy = G.Y0 + Math.floor(i / 2) * 0.98;
      eq(s, n, x, yy, i === 3 ? 0.42 : sc, l);
    });
    const items = [
      { k: "b", t: "Every term in G (Eq. 2) is a derivative of measured velocity. Du/Dt comes from three snapshots, at a fixed point (Eq. 3) or along the particle path (Eq. 4)." },
      { k: "b", t: "∇· of Eq. 2 gives a Poisson problem with Neumann edges (Eq. 5): solvable only if Eq. 6 holds, and only up to a constant, so we compare p − p̄." },
      { k: "head", t: "Assumptions" },
      { k: "text", t: "Incompressible Newtonian fluid, constant ρ and μ; 2-D planar flow (out-of-plane terms neglected); velocity on a uniform grid at instants n − 1, n, n + 1." },
    ];
    para(s, items, { x: 0.45, y: 4.72, w: 9.1, h: 2.45 }, { size: 12.5, gap: 5, bgap: 3 });
  }
  s.addNotes("Eq. 2 is the key: pressure gradient from velocity. Eulerian vs Lagrangian matters when structures are convected. The Poisson problem with pure Neumann data is singular; we fix the constant with a zero-mean condition.");
}
};

BL.algorithms = () => {
{
  const s = wide ? newSlide("Reconstruction Algorithms") : newSlide("Algorithms", "Three ways to integrate ∇p = G");
  const fh = wide ? 2.9 : 2.55;
  fig(s, "f_methods", { x: G.X0, y: G.Y0, w: G.CW, h: fh + (wide ? 0.25 : 0) }, "The three integration algorithms on the PIV grid", { border: false });
  const cw = (G.CW - 0.5) / 3, ty = G.Y0 + fh + (wide ? 0.4 : 0.12), th = G.Y1 - ty;
  const cols = [
    [{ k: "head", t: "(a) Direct integration" }, { k: "text", t: "March G·dl from a reference corner with the trapezoidal rule, x-then-y and y-then-x; average the two paths." }, { k: "note", t: "Fastest (0.3 ms at 64²), but every error on the path is carried to every later point." }],
    [{ k: "head", t: "(b) Omni-directional" }, { k: "text", t: "Integrate along 32 ray directions into each node (Liu & Katz 2006). Edge values come from a consistent linear system; interior values are ray averages." }, { k: "note", t: "Averaging cancels path errors; slowest (≈0.4 s)." }],
    [{ k: "head", t: "(c) Pressure-Poisson (main)" }, { k: "text", t: "5-point Laplacian as a sparse matrix; Neumann edges by ghost-node elimination; singular system closed by a zero-mean constraint, solved directly." }, { k: "note", t: "15 ms; also Jacobi, GS, SOR and CG (backup)." }],
  ];
  cols.forEach((c, i) => para(s, c, { x: G.X0 + i * (cw + 0.25), y: ty, w: cw, h: th }, { size: wide ? 13 : 11.5, headSize: wide ? 14 : 13, align: "left" }));
  s.addNotes("Direct: simplest but path dependent. Omni: average over many paths, robust but expensive. Poisson: the standard method; elliptic, so errors are smoothed. All three use the same G.");
}
};

BL.strategy = () => {
{
  const s = wide ? newSlide("Solution Strategy: Logic Chart and Pseudo-code") : newSlide("Solution Strategy", "Logic chart and pseudo-code");
  const fx = G.X0, fw = wide ? 4.6 : 3.7, bh = wide ? 0.5 : 0.5, gap = wide ? 0.2 : 0.17;
  const steps = [
    "Velocity source: exact solutions (Taylor–Green, Kovasznay) or CFD wake",
    "Synthetic PIV: box filter W → subsample to h → noise σU",
    "Optional denoising: Gaussian / Savitzky–Golay",
    "Material derivative Du/Dt: Eulerian or Lagrangian",
    "Pressure gradient G = −ρ Du/Dt + μ∇²u",
    "Integrate: Direct | Omni | Poisson",
    "Compare with exact p: ε_L2, error map, core pressure",
  ];
  let y = G.Y0;
  steps.forEach((t, i) => {
    fbox(s, t, fx, y, fw, bh, { size: wide ? 11 : 10, fill: i === 5 ? "1F4E79" : "4F81BD", fillL: i === 5 ? "DCEBF2" : undefined });
    if (i < steps.length - 1) arrowDown(s, fx + fw / 2, y + bh, gap);
    y += bh + gap;
  });
  // loop back arrow (sweep)
  const lx = fx + fw + 0.12, yTop = G.Y0 + bh + gap + bh / 2, yBot = y - gap - bh / 2;
  s.addShape(pres.shapes.LINE, { x: fx + fw, y: yBot, w: 0.2, h: 0, line: { color: "C0504D", width: 1.5 } });
  s.addShape(pres.shapes.LINE, { x: fx + fw + 0.2, y: yTop, w: 0, h: yBot - yTop, line: { color: "C0504D", width: 1.5 } });
  s.addShape(pres.shapes.LINE, { x: fx + fw, y: yTop, w: 0.2, h: 0, line: { color: "C0504D", width: 1.5, beginArrowType: "triangle" } });
  s.addText("sweep h, σ, Δt, BC, W, filter;\n10 noise seeds each", { x: fx + fw + 0.26, y: (yTop + yBot) / 2 - 0.4, w: wide ? 1.25 : 1.0, h: 0.8, fontFace: BODY, fontSize: 9.5, color: "C0504D", margin: 0, isTextBox: true });
  figNo += 1;
  s.addText(`Fig. ${figNo} Logic chart of the code`, { x: fx, y: y - gap + 0.05, w: fw, h: 0.26, fontFace: BODY, fontSize: S.cap, color: C.cap, align: "center", margin: 0, isTextBox: true });
  // pseudo-code
  const px = fx + fw + (wide ? 1.6 : 1.3), pw = G.X1 - px;
  const code = [
    "for case in [TaylorGreen, Kovasznay, CFD]:",
    "  truth = exact_or_cfd(case)   # fine grid",
    "  for param in sweep(h, σ, Δt, BC, W, filter):",
    "    for seed in range(10):",
    "      F = degrade(truth, W, h, σ, seed)",
    "      F = denoise(F)            # optional",
    "      a = DuDt(F, eulerian | lagrangian)",
    "      G = -ρ*a + μ*lap(u)",
    "      p = poisson(G) | direct(G) | omni(G)",
    "      ε = L2(p - p_ref) / L2(p_ref)  # zero mean",
    "    save mean, std of ε -> CSV -> figures",
  ];
  const ch = wide ? 2.45 : 2.0;
  s.addShape(pres.shapes.RECTANGLE, { x: px, y: G.Y0, w: pw, h: ch, fill: { color: wide ? "F2F2F2" : "FFFFFF", transparency: wide ? 0 : 15 }, line: { color: "7F7F7F", width: 0.75 } });
  s.addText(code.map((t, i) => ({ text: t, options: { breakLine: i < code.length - 1 } })), { x: px + 0.1, y: G.Y0 + 0.08, w: pw - 0.2, h: ch - 0.16, fontFace: "Courier New", fontSize: wide ? 11 : 9, color: "1F1F1F", margin: 0, valign: "top", isTextBox: true });
  const mods = [
    { k: "head", t: "Main functions (Python, NumPy/SciPy)" },
    { k: "b", t: [B("fields/"), R(" exact solutions, CFD loader · "), B("piv/"), R(" degrade, denoise")] },
    { k: "b", t: [B("core/"), R(" derivatives, material derivative, Poisson (sparse direct + iterative), direct and omni integration, reconstruct()")] },
    { k: "b", t: [B("analysis/"), R(" metrics, plots · "), B("studies/"), R(" one script per sweep")] },
    { k: "b", t: [B("tests/"), R(" 87 unit tests, all passing (manufactured solutions, order checks)")] },
  ];
  para(s, mods, { x: px, y: G.Y0 + ch + 0.15, w: pw, h: G.Y1 - G.Y0 - ch - 0.15 }, { size: wide ? 13 : 11, headSize: wide ? 14 : 12.5, align: "left", bgap: wide ? 5 : 3 });
  s.addNotes("One function, reconstruct(), runs the whole chain. Every study script loops over a parameter and 10 noise seeds and writes a CSV; the figures are generated from those CSVs.");
}
};

BL.testcases = () => {
{
  const s = wide ? newSlide("Test Cases, Synthetic PIV and Error Metric") : newSlide("Test Cases", "Exact flows, synthetic PIV and the error metric");
  const items = [
    { k: "head", t: "Exact Navier–Stokes solutions (Fig. " + (figNo + 1) + ")" },
    { k: "b", t: [B("Taylor–Green vortex: "), R("unsteady and decaying; the viscous term cancels, so ∇p balances convection. A copy carried at U = 1 tests the time derivative.")] },
    { k: "b", t: [B("Kovasznay flow (Re = 20, 40): "), R("steady flow behind a grid; here the viscous term μ∇²u matters.")] },
    { k: "b", t: [B("CFD cylinder wake (Re = 100): "), R("OpenFOAM solution sampled on an 81 × 81 window behind the cylinder (Fig. " + (figNo + 2) + ").")] },
    { k: "head", t: "Synthetic PIV model" },
    { k: "text", t: "Exact field on a 4× finer grid → box filter of width W (interrogation window) → keep every 4th point (spacing h) → add Gaussian noise σ·U to u and v." },
  ];
  if (wide) {
    fig(s, "f_testcases", { x: G.X0, y: G.Y0 + 0.05, w: 6.4, h: 3.6 }, "Test flows: pressure (colour) and streamlines");
    para(s, items, { x: 7.35, y: G.Y0, w: G.X1 - 7.35, h: 3.75 }, { gap: 8 });
    para(s, [{ k: "head", t: "Error metric" }], { x: 7.35, y: 4.95, w: 4, h: 0.3 });
    eq(s, "eq_eps", 7.45, 5.32, 0.62);
    para(s, [{ k: "note", t: "Baseline: 64 × 64 vectors, W = h, σ = 1 % of U, Δt = 0.05; mean ± std over 10 noise seeds." }], { x: 7.35, y: 6.35, w: G.X1 - 7.35, h: 0.5 });
    fig(s, "c_cfdp", { x: G.X0, y: 4.85, w: 6.4, h: 2.3 }, "CFD wake window (x/D 1–5), pressure at t = 187.5", { capH: 0.25, capSize: 11, border: false });
  } else {
    fig(s, "f_testcases", { x: 0.4, y: G.Y0, w: 4.6, h: 2.55 }, "Test flows: pressure and streamlines", { capSize: 10 });
    fig(s, "c_cfdp", { x: 0.4, y: 4.45, w: 4.6, h: 1.95 }, "CFD wake window (x/D 1–5)", { border: false, capSize: 10, capH: 0.25 });
    para(s, items, { x: 5.2, y: G.Y0, w: 4.4, h: 3.25 }, { size: 11, headSize: 12.5, gap: 4 });
    para(s, [{ k: "head", t: "Error metric" }], { x: 5.2, y: 5.0, w: 4, h: 0.3 }, { headSize: 13 });
    eq(s, "eq_eps", 5.3, 5.38, 0.45);
    para(s, [{ k: "note", t: "Baseline: 64 × 64 vectors, W = h, σ = 1 % of U, Δt = 0.05; mean ± std over 10 noise seeds." }], { x: 5.2, y: 6.35, w: 4.4, h: 0.7 });
  }
  s.addNotes("Exact solutions give the true pressure everywhere, so the error is known exactly. The degradation model mimics PIV: spatial averaging over the window, finite spacing, random error.");
}
};

BL.verif = () => {
{
  const s = wide ? newSlide("Verification: Grid Convergence on Clean Data") : newSlide("Verification", "Grid convergence on clean (exact) data");
  const fh = wide ? 3.65 : 3.0;
  fig(s, "f_conv_clean", { x: G.X0, y: G.Y0, w: G.CW, h: fh }, "Error vs PIV spacing h, exact point samples (W = 0, σ = 0)", { border: false });
  const items = [
    { k: "head", t: "1. All three methods are second-order accurate" },
    { k: "b", t: "Fitted slopes 1.82–2.13: central differences and the trapezoidal rule are O(h²), and the Poisson solve adds no lower-order error." },
    { k: "head", t: "2. Boundaries do not spoil the order" },
    { k: "b", t: "Second-order one-sided stencils and ghost-node Neumann conditions at the edges; a first-order version gives slope ≈ 1 (backup)." },
    { k: "head", t: "3. The viscous term matters at low Re" },
    { k: "b", t: "Without μ∇²u, Kovasznay stops converging (unit test). The code is verified before noise is added." },
  ];
  if (wide) {
    para(s, items.slice(0, 4), { x: G.X0, y: G.Y0 + fh + 0.1, w: G.CW / 2 - 0.2, h: 2.2 });
    para(s, items.slice(4), { x: G.X0 + G.CW / 2 + 0.2, y: G.Y0 + fh + 0.1, w: G.CW / 2 - 0.2, h: 2.2 });
  } else {
    para(s, items, { x: 0.45, y: G.Y0 + fh + 0.08, w: 9.1, h: G.Y1 - G.Y0 - fh - 0.08 }, { size: 12, headSize: 13, gap: 4, bgap: 1 });
  }
  s.addNotes("Before studying noise we check the code is correct. All methods converge at second order on all three flows.");
}
};

BL.noiseres = () => {
{
  const s = wide ? newSlide("Noise vs. Resolution: An Optimum Spacing") : newSlide("Noise vs. Resolution", "Finer is not always better");
  const items = [
    { k: "head", t: "1. With 1 % noise, refining the grid stops helping" },
    { k: "b", t: "Poisson error is lowest at h ≈ 0.032–0.035 (n = 64) for all three flows and rises on finer grids: Taylor–Green 0.082 → 0.20 at h = 0.009." },
    { k: "head", t: "2. Why: two competing errors" },
    { k: "b", t: "Truncation falls as h², but differentiation amplifies noise: σ/h in (u·∇)u and μσ/h² in the viscous term. The sum has a minimum." },
    { k: "head", t: "3. Direct integration only gets worse" },
    { k: "b", t: "0.26 → 2.2 as h shrinks: a path crosses ~1/h noisy steps, so its error grows like a random walk." },
  ];
  if (wide) {
    fig(s, "f_conv_noisy", { x: G.X0, y: G.Y0, w: 8.7, h: 3.75 }, "Error vs spacing with σ = 1 %, W = h (mean ± std, 10 seeds)", { border: false });
    fig(s, "f_tradeoff", { x: 9.55, y: G.Y0, w: 3.15, h: 3.0 }, "Error budget (schematic)", { border: false });
    para(s, items.slice(0, 4), { x: G.X0, y: 5.05, w: 6.0, h: 2.1 });
    para(s, items.slice(4), { x: 6.85, y: 5.05, w: 5.88, h: 2.1 });
  } else {
    fig(s, "f_conv_noisy", { x: 0.4, y: G.Y0, w: 9.2, h: 3.0 }, "Error vs spacing with σ = 1 %, W = h (mean ± std, 10 seeds)", { border: false });
    fig(s, "f_tradeoff", { x: 7.15, y: 4.82, w: 2.45, h: 2.35 }, "Error budget", { border: false, capSize: 10 });
    para(s, items, { x: 0.45, y: 4.78, w: 6.55, h: 2.42 }, { size: 11.5, headSize: 12.5, gap: 3, bgap: 1 });
  }
  s.addNotes("Main practical result: there is an optimum vector spacing. Below it, noise amplification by the derivatives dominates truncation error.");
}
};

BL.noisemethod = () => {
{
  const s = wide ? newSlide("Noise Sensitivity and Choice of Method") : newSlide("Noise & Method", "How each algorithm handles measurement noise");
  const items = [
    { k: "b", t: [B("Error grows linearly with noise "), R("(log–log slopes 0.97–1.01): the reconstruction is a linear operation on the velocity error.")] },
    { k: "b", t: [B("Direct integration is 5–7× worse: "), R("one path carries every error along it, while Poisson (elliptic) and omni average over many paths.")] },
    { k: "b", t: [B("Omni ≈ Poisson in accuracy "), R("(0.063 vs 0.082 Taylor–Green; 0.067 vs 0.059 Kovasznay) but ≈27× slower (0.41 s vs 15 ms).")] },
    { k: "b", t: [B("Velocity-source Poisson "), R("∇²p = −ρ(u_x² + 2u_y v_x + v_y²) is 5–8× worse with noise: it multiplies noisy derivatives.")] },
  ];
  if (wide) {
    fig(s, "f_noise", { x: G.X0, y: G.Y0, w: 6.9, h: 3.6 }, "Error vs noise level (n = 64, W = h)", { border: false });
    fig(s, "f_method_bar", { x: 7.45, y: G.Y0, w: 5.28, h: 3.6 }, "Method comparison at σ = 1 %", { border: false });
    para(s, items, { x: G.X0, y: 4.95, w: G.CW, h: 2.1 }, { bgap: 6, size: 14 });
  } else {
    fig(s, "f_noise", { x: 0.4, y: G.Y0, w: 4.85, h: 2.6 }, "Error vs noise (n = 64, W = h)", { border: false, capSize: 10 });
    fig(s, "f_method_bar", { x: 5.35, y: G.Y0, w: 4.25, h: 2.6 }, "Methods at σ = 1 %", { border: false, capSize: 10 });
    para(s, items, { x: 0.45, y: 4.45, w: 9.1, h: 2.75 }, { size: 12, bgap: 4 });
  }
  s.addNotes("Linear growth with sigma is expected. The important comparison is between methods: averaging methods are far more robust than a single integration path.");
}
};

BL.bc = () => {
{
  const s = wide ? newSlide("Boundary Conditions for the Poisson Solver") : newSlide("Boundary Conditions", "What we impose on the edges matters");
  const items = [
    { k: "head", t: "1. Known pressure on edges helps" },
    { k: "b", t: "Exact p on 2 opposite edges halves the error (Taylor–Green 0.082 → 0.044; Kovasznay 0.059 → 0.033); 1 edge gives little gain." },
    { k: "head", t: "2. Pure Neumann: fixing the constant matters" },
    { k: "b", t: "Pinning one node puts that node’s error into the whole field (0.112 vs 0.082 with a zero-mean constraint)." },
    { k: "head", t: "3. Bernoulli edges fail in rotational flow" },
    { k: "b", t: "p + ½ρ|u|² = const holds only in irrotational flow; both test flows are rotational, giving 0.30–0.81. Use it only at free-stream edges." },
    { k: "note", t: "Why: Neumann data n·G are themselves noisy derivatives, while a Dirichlet edge anchors the solution." },
  ];
  if (wide) {
    fig(s, "f_bc_bar", { x: G.X0, y: G.Y0 + 0.2, w: 6.6, h: 4.3 }, "Poisson error for six boundary treatments, σ = 1 %", { border: false });
    para(s, items, { x: 7.45, y: G.Y0 + 0.2, w: G.X1 - 7.45, h: 5.6 }, { gap: 10, bgap: 3 });
  } else {
    fig(s, "f_bc_bar", { x: 0.4, y: G.Y0, w: 5.1, h: 3.3 }, "Poisson error, six boundary treatments", { border: false, capSize: 10 });
    para(s, items.slice(0, 4), { x: 5.65, y: G.Y0, w: 3.95, h: 3.35 }, { size: 11.5, headSize: 12.5, gap: 5, align: "left" });
    para(s, items.slice(4), { x: 0.45, y: 4.65, w: 9.1, h: 2.0 }, { size: 12, headSize: 13, gap: 4, align: "left" });
  }
  s.addNotes("Boundary data are the single most effective improvement. But Bernoulli is only valid where the flow is irrotational; inside a vortex it is simply wrong.");
}
};

BL.dt = () => {
{
  const s = wide ? newSlide("Time Step: Eulerian vs. Lagrangian Du/Dt") : newSlide("Time Step", "Eulerian vs. Lagrangian material derivative");
  const items = [
    { k: "head", t: "1. Convected vortex, clean data" },
    { k: "b", t: "Eulerian error rises 10× by Δt = 0.2 (1.8e-3 → 1.8e-2): at a fixed point the pattern sweeps past, so ∂u/∂t changes fast. Lagrangian follows the particle and stays 3× lower (5.8e-3)." },
    { k: "head", t: "2. Stationary vortex" },
    { k: "b", t: "Eulerian error does not depend on Δt; Lagrangian grows (4.2e-3 at Δt = 0.2) because curved paths are traced with a straight step." },
    { k: "head", t: "3. Noisy data" },
    { k: "b", t: "Noise enters ∂u/∂t as σ/Δt, so both improve with larger Δt (0.27 → 0.07). The best Δt balances noise against truncation." },
  ];
  if (wide) {
    fig(s, "f_dt", { x: G.X0, y: G.Y0, w: G.CW, h: 3.95 }, "Error vs time between snapshots", { border: false });
    [0, 2, 4].forEach((k, i) => para(s, items.slice(k, k + 2), { x: G.X0 + i * 4.13, y: 5.15, w: 3.87, h: 2.0 }, { bgap: 3, align: "left" }));
  } else {
    fig(s, "f_dt", { x: 0.4, y: G.Y0, w: 6.2, h: 2.95 }, "Error vs time between snapshots", { border: false, capSize: 10 });
    para(s, items.slice(0, 2), { x: 6.75, y: G.Y0 + 0.05, w: 2.85, h: 2.9 }, { size: 11, headSize: 12, align: "left" });
    para(s, items.slice(2), { x: 0.45, y: 4.82, w: 9.1, h: 2.3 }, { size: 11.5, headSize: 12.5, gap: 5, align: "left" });
  }
  s.addNotes("Lagrangian wins when structures are convected through the window, as in a wake. With noise, a longer time separation reduces the amplification of noise in the time derivative.");
}
};

BL.window = () => {
{
  const s = wide ? newSlide("Interrogation Window and Denoising") : newSlide("Window & Denoising", "Smoothing trades noise for bias");
  const items = [
    { k: "head", t: "1. Interrogation window W" },
    { k: "b", t: "Averaging over W smooths velocity. Taylor–Green core suction drops only 0.3 % at W = 4h (vortex ≫ window); Kovasznay ε rises 0.015 → 0.061 (sharper gradients)." },
    { k: "head", t: "2. Denoising cuts the noisy error about 3×" },
    { k: "b", t: "σ = 1 %: 0.081 → 0.027 (Gaussian), 0.031 (Savitzky–Golay); but on clean data the Gaussian filter adds bias (0.0002 → 0.007)." },
    { k: "head", t: "3. Noise biases the vortex core, not just scatters it" },
    { k: "b", t: "Without denoising the core suction is over-predicted (1.39× at σ = 10 %): (u·∇)u multiplies noisy terms, and products of noise do not average to zero." },
  ];
  if (wide) {
    fig(s, "f_denoise", { x: G.X0, y: G.Y0, w: G.CW, h: 3.95 }, "Taylor–Green core window: suction ratio and field error vs noise", { border: false });
    [0, 2, 4].forEach((k, i) => para(s, items.slice(k, k + 2), { x: G.X0 + i * 4.13, y: 5.15, w: 3.87, h: 2.0 }, { bgap: 3, align: "left" }));
  } else {
    fig(s, "f_denoise", { x: 0.4, y: G.Y0, w: 6.2, h: 2.95 }, "Core suction ratio and field error vs noise", { border: false, capSize: 10 });
    para(s, items.slice(0, 2), { x: 6.75, y: G.Y0 + 0.05, w: 2.85, h: 2.9 }, { size: 11, headSize: 12, align: "left" });
    para(s, items.slice(2), { x: 0.45, y: 4.82, w: 9.1, h: 2.3 }, { size: 11.5, headSize: 12.5, gap: 5, align: "left" });
  }
  s.addNotes("Filtering is a trade-off: it removes random error but introduces systematic smoothing error. Worth it once noise is above about 0.5 %.");
}
};

BL.cfdsetup = () => {
{
  const s = wide ? newSlide("CFD Simulation: Cylinder Wake at Re = 100") : newSlide("CFD Simulation", "Cylinder wake at Re = 100 (OpenFOAM)");
  const setRows = [["Item", "Setting"], ["Equations", "Incompressible NS (Eq. 1), laminar"], ["Solver", "pimpleFoam (PISO mode)"], ["Time / convection", "backward (2nd order) / linearUpwind"], ["Time step", "max Courant 0.5, run to t = 200 D/U"], ["Re = UD/ν", "100 (U = D = 1, ν = 0.01)"], ["Output", "Cd, Cl (forces); PIV window every 0.05 D/U"]];
  const meshRows = [["Mesh", "Cells", "1st cell", "St", "C̄d", "Cl amp."], ["Coarse", "6,690", "0.022", "–", "–", "–"], ["Medium", "26,760", "0.011", "–", "–", "–"], ["Fine", "107,040", "0.0055", "–", "–", "–"], [{ text: "Park et al. (1998)", options: { italic: true } }, "–", "–", "≈0.165", "≈1.33", "–"]];
  if (wide) {
    fig(s, "f_cfd_domain", { x: G.X0, y: G.Y0, w: 5.9, h: 3.85 }, "Domain, boundary conditions and PIV window (not to scale)", { border: false });
    para(s, [{ k: "b", t: "2-D laminar is valid: the wake stays two-dimensional below Re ≈ 190." }, { k: "b", t: "A small transverse inlet kick (t < 6) starts shedding sooner; after t = 6 the inlet is exactly U = (1, 0)." }], { x: G.X0, y: 5.05, w: 5.9, h: 1.0 }, { size: 12 });
    table(s, setRows, { x: 6.85, y: G.Y0, w: 5.88, colW: [1.75, 4.13], size: 11, rowH: 0.27 });
    fig(s, "c_paraview", { x: 6.85, y: 3.2, w: 5.88, h: 3.95 }, "ParaView 6.2: fine-mesh pressure at t = 200, Kármán vortex street", { capH: 0.28 });
  } else {
    fig(s, "f_cfd_domain", { x: 0.4, y: G.Y0, w: 4.4, h: 2.85 }, "Domain, BCs and PIV window", { border: false, capSize: 10 });
    table(s, setRows, { x: 5.0, y: G.Y0, w: 4.6, colW: [1.35, 3.25], size: 9.5, rowH: 0.25 });
    fig(s, "c_paraview", { x: 4.6, y: 3.95, w: 5.0, h: 3.2 }, "ParaView 6.2: fine-mesh pressure at t = 200", { border: false, capSize: 10, capH: 0.25 });
    para(s, [{ k: "b", t: "2-D laminar is valid: the wake stays 2-D below Re ≈ 190." }, { k: "b", t: "A small inlet kick (t < 6) starts shedding sooner; from t = 6 the inlet is exactly U = (1, 0)." }, { k: "b", t: "Three meshes (6,690 / 26,760 / 107,040 cells) for grid independence." }], { x: 0.45, y: 4.75, w: 4.0, h: 2.3 }, { size: 11, bgap: 4, align: "left" });
  }
  s.addNotes("OpenFOAM set-up. The ParaView screenshot shows the actual fine-mesh solution at t = 200: high pressure at the stagnation point, low pressure on the sides and in the shed vortex cores, which alternate above and below the centreline.");
}
};

BL.cfdgrid = () => {
{
  const s = wide ? newSlide("CFD: Grid Independence and Periodicity") : newSlide("CFD Verification", "Grid independence and periodic shedding");
  const gRows = [["Mesh", "Cells", "St", "C̄d", "Cl amp."],
    ["Coarse", "6,690", "0.1644", "1.3754", "0.3216"],
    ["Medium", "26,760", "0.1688", "1.3880", "0.3451"],
    [{ text: "Fine", options: { bold: true } }, { text: "107,040", options: { bold: true } }, { text: "0.1695", options: { bold: true } }, { text: "1.3879", options: { bold: true } }, { text: "0.3451", options: { bold: true } }],
    ["Richardson extrap.", "–", "0.1697", "conv.", "conv."],
    [{ text: "2-D simulations, Re = 100", options: { italic: true } }, "–", "≈0.164–0.165", "≈1.33", "≈0.33"]];
  const items = [
    { k: "num", t: [B("St converges: "), R("0.1644 → 0.1688 → 0.1695 (changes 2.7 %, then 0.46 %). Observed order p = 2.48, extrapolated 0.1697, GCI fine 0.13 %; asymptotic ratio 1.005, so the meshes are in the asymptotic range.")] },
    { k: "num", t: [B("C̄d and Cl amplitude "), R("are already converged on the medium mesh (changes 0.00 % and 0.01 %), so a Richardson order cannot be defined.")] },
    { k: "num", t: [B("Periodic from t = 76.3: "), R("cycle-to-cycle changes of amplitude and period stay below 10⁻³; the PIV window (t = 175–200) is well inside. Residuals level off near 10⁻³: periodic, not steady.")] },
    { k: "num", t: [B("vs literature: "), R("St and C̄d are 3–4 % higher, as expected from the 5 % blockage of the ±10 D symmetry sides.")] },
  ];
  if (wide) {
    table(s, gRows, { x: G.X0, y: G.Y0 + 0.05, w: 6.0, colW: [2.1, 0.95, 1.05, 0.95, 0.95], size: 11.5, rowH: 0.29 });
    fig(s, "c_st", { x: G.X0, y: 3.05, w: 6.0, h: 2.6 }, "St and mean Cd vs mesh; error bars = GCI", { border: false });
    fig(s, "c_periodic", { x: 6.95, y: G.Y0 + 0.05, w: 5.78, h: 2.75 }, "Fine-mesh lift and the periodicity criterion", { border: false });
    para(s, items, { x: 6.95, y: 4.0, w: 5.78, h: 3.1 }, { size: 12, bgap: 4, align: "left" });
    para(s, [{ k: "note", t: "Literature values approximate (Park, Kwon & Choi 1998; Williamson 1996); check before quoting." }], { x: G.X0, y: 6.15, w: 6.0, h: 0.5 }, { noteSize: 10.5 });
  } else {
    table(s, gRows, { x: 0.4, y: G.Y0, w: 4.95, colW: [1.75, 0.8, 0.85, 0.8, 0.75], size: 9, rowH: 0.24 });
    fig(s, "c_st", { x: 0.4, y: 3.35, w: 4.95, h: 2.15 }, "St and Cd vs mesh (bars: GCI)", { border: false, capSize: 10, capH: 0.25 });
    fig(s, "c_periodic", { x: 5.5, y: G.Y0, w: 4.1, h: 1.95 }, "Lift and periodicity check", { border: false, capSize: 10, capH: 0.25 });
    para(s, items, { x: 5.5, y: 3.8, w: 4.1, h: 3.3 }, { size: 10.5, bgap: 3, align: "left" });
    para(s, [{ k: "note", t: "Literature values approximate; check before quoting." }], { x: 0.45, y: 5.6, w: 4.9, h: 0.4 }, { noteSize: 10 });
  }
  s.addNotes("Grid independence: the Strouhal number converges with observed order about 2.5 and a fine-mesh GCI of 0.13 %. Drag and lift amplitude stop changing at the medium mesh. The flow is periodic after t = 76, so the sampling window from t = 175 is fully developed. Our values are a few per cent above the unbounded-domain literature because the side walls block 5 % of the flow.");
}
};

BL.cfdrecon = () => {
{
  const s = wide ? newSlide("Pressure Reconstructed from the CFD Wake") : newSlide("CFD Wake: Reconstruction", "Pressure from PIV-like samples of the OpenFOAM velocity");
  const items = [
    { k: "num", t: [B("The wake pressure is recovered: "), R("from the 81 × 81 CFD velocity samples (h = 0.05 D, no noise) the Poisson error is ε = 0.5–0.75 % over a shedding period.")] },
    { k: "num", t: [B("Vortex-core suction "), R("(−0.22 ρU² at x/D ≈ 2.7) is matched by every method within about 1 %.")] },
    { k: "num", t: [B("Error sits on the CFD mesh's block edges "), R("(x/D = 2 and the diagonals to (2, ±2), the edges of the O-grid). The interpolated velocity has small slope jumps there, and differentiation amplifies them; direct integration shows a kink at x/D = 2. This is an error of the \"measured\" data, not of the method.")] },
  ];
  if (wide) {
    fig(s, "c_fields", { x: G.X0, y: G.Y0, w: G.CW, h: 3.45 }, "CFD pressure, reconstruction (Poisson, Eulerian) and error map at t = 187.5", { border: false });
    fig(s, "c_line", { x: G.X0, y: 4.65, w: 4.4, h: 2.5 }, "Line through a vortex core (y/D = 0.35)", { border: false, capH: 0.25, capSize: 11 });
    para(s, items, { x: 5.3, y: 4.7, w: G.X1 - 5.3, h: 2.45 }, { size: 12.5, bgap: 4, align: "left" });
  } else {
    fig(s, "c_fields", { x: 0.4, y: G.Y0, w: 9.2, h: 2.95 }, "CFD pressure, reconstruction and error map (t = 187.5)", { border: false, capSize: 10 });
    fig(s, "c_line", { x: 0.4, y: 4.75, w: 4.1, h: 2.4 }, "Line through a vortex core", { border: false, capSize: 10, capH: 0.25 });
    para(s, items, { x: 4.65, y: 4.75, w: 4.95, h: 2.4 }, { size: 10.5, bgap: 3, align: "left" });
  }
  s.addNotes("Reconstruction on the CFD wake. Without noise the error is below 1 %. Note where the error lives: on the block boundaries of the OpenFOAM mesh, where the sampled velocity has small kinks; differentiation turns kinks into steps. This is the same mechanism as noise amplification, applied to a systematic error.");
}
};

BL.cfddt = () => {
{
  const s = wide ? newSlide("CFD Wake: Time Step, Scheme and Noise") : newSlide("CFD Wake: Δt and Noise", "Eulerian vs Lagrangian on a real wake");
  const items = [
    { k: "num", t: [B("Error pulses twice per period "), R("(T = 5.90, i.e. St = 0.1695): it peaks each time a vortex forms on one side.")] },
    { k: "num", t: [B("At Δt = 0.05 Lagrangian is 15–18 % better "), R("throughout the period: the vortices are carried through the window.")] },
    { k: "num", t: [B("At Δt ≥ 0.1 Eulerian is slightly better, "), R("unlike the convected Taylor–Green. Likely: in the near wake (x/D 1–5) vortices are still forming, move slower than U and spin fast (w ≈ c), so straight-line tracking of curved paths costs more; more tracked points also leave the window.")] },
    { k: "num", t: [B("With 1 % noise "), R("ε ≈ 0.13–0.15 for both schemes, nearly flat in Δt. Higher than Taylor–Green (0.08) because wake pressure varies only about ±0.2 ρU², so the same velocity noise is larger relative to the signal.")] },
  ];
  if (wide) {
    fig(s, "c_period", { x: G.X0, y: G.Y0, w: 5.7, h: 3.65 }, "Error over one shedding period (clean, Δt = 0.05)", { border: false });
    fig(s, "c_dt_cfd", { x: 6.45, y: G.Y0, w: 2.95, h: 5.95 }, "Error vs Δt on the wake", { border: false, capSize: 11 });
    para(s, items.slice(0, 2), { x: G.X0, y: 4.95, w: 5.7, h: 2.15 }, { size: 12.5, bgap: 5, align: "left" });
    para(s, items.slice(2), { x: 9.6, y: G.Y0, w: G.X1 - 9.6, h: 6.0 }, { size: 12, bgap: 8, align: "left", startAt: 3 });
  } else {
    fig(s, "c_period", { x: 0.4, y: G.Y0, w: 5.6, h: 2.85 }, "Error over one shedding period (clean, Δt = 0.05)", { border: false, capSize: 10, capH: 0.25 });
    fig(s, "c_dt_cfd", { x: 6.2, y: G.Y0, w: 3.4, h: 5.35 }, "Error vs Δt on the wake", { border: false, capSize: 10 });
    para(s, items, { x: 0.45, y: 4.65, w: 5.6, h: 2.45 }, { size: 10, bgap: 2, align: "left" });
  }
  s.addNotes("Time-resolved behaviour on the wake. The Lagrangian advantage seen on the convected Taylor–Green holds only at the smallest time step here; in the near wake the vortices are not yet simply convected, so the theory's assumption (w much smaller than c) does not hold. With noise, noise amplification dominates and the scheme matters little.");
}
};

BL.lit = () => {
{
  const s = wide ? newSlide("Comparison with Published Results") : newSlide("Literature", "Do our trends match published studies?");
  const ok = { text: "Agrees", options: { color: "1E7B34", bold: true } };
  const rows = [
    ["Our result", "Published finding", "Match"],
    ["Exact p on 2 edges: ε 0.044 vs 0.082 all-Neumann (≈2×)", "Pan et al. [3]: Dirichlet edges give less error than Neumann (≈2× for large domains); avoid pure Neumann if possible", ok],
    ["Bernoulli edges: ε 0.30–0.81 in rotational flow", "Pan et al. [3]: Bernoulli-based Dirichlet data are valid only in irrotational regions", ok],
    ["Error from 0.02 % (clean) to >100 % (direct, σ ≥ 2 %); depends on method, spacing and filtering", "Charonko et al. [2]: errors from <1 % to >100 %; sensitive to all these factors; no single best method", ok],
    ["Lagrangian 3× better at large Δt in convected flow; noise ∝ 1/Δt", "Violato et al. [6]: Eulerian limited to short Δt; short Δt amplifies precision error, long Δt adds truncation error", ok],
    ["Omni and Poisson recover the main features; direct integration is very noise-sensitive", "van Gent et al. [4]: all methods tested recover the main features; noise lowers accuracy", { text: "Partly", options: { color: "B26B00", bold: true } }],
    ["CFD wake (fine mesh): St 0.1695, C̄d 1.388, Cl amplitude 0.345", "2-D simulations at Re = 100 (e.g. Park et al. [7]): St ≈ 0.164–0.165, C̄d ≈ 1.33; experiments St ≈ 0.164 [8]", { text: "Agrees (+3–4 %, blockage)", options: { color: "1E7B34", bold: true } }],
  ];
  if (wide) table(s, rows, { x: G.X0, y: G.Y0 + 0.1, w: G.CW, colW: [4.75, 5.98, 1.4], size: 13, rowH: 0.72 });
  else table(s, rows, { x: 0.4, y: G.Y0, w: 9.2, colW: [3.55, 4.55, 1.1], size: 10, rowH: 0.62 });
  para(s, [{ k: "note", t: "Reference numbers refer to the References slide." }], { x: G.X0, y: G.Y1 - 0.35, w: G.CW, h: 0.3 });
  s.addNotes("This is the research part of the project: our controlled synthetic experiments reproduce the main published conclusions, and the CFD wake matches published Strouhal and drag values within the few per cent expected from blockage. Partial: van Gent et al. found no method critically sensitive to noise, whereas our direct integration is; their noise levels and methods differ.");
}
};

BL.concl = () => {
{
  const s = wide ? newSlide("Conclusions") : newSlide("Conclusions", "Key findings and takeaway");
  const left = [
    { k: "head", t: "Key findings" },
    { k: "b", t: "All three algorithms are second-order accurate on exact data (slopes 1.8–2.1): the code is verified." },
    { k: "b", t: "Noise, not discretisation, limits PIV pressure: error grows linearly with σ, and there is an optimum spacing (h ≈ 0.03 here)." },
    { k: "b", t: "Poisson and omni are 5–7× more robust than direct integration; Poisson is ≈27× faster than omni." },
    { k: "b", t: "Known pressure on two edges halves the error; Bernoulli edges hold only in irrotational flow." },
    { k: "b", t: "Lagrangian Du/Dt tolerates larger Δt in convected flow; with noise, larger Δt helps both." },
    { k: "b", t: "Filtering trades noise for bias: it helps at σ ≥ 1 % and hurts on clean data." },
    { k: "b", t: "On the OpenFOAM wake (St 0.1695, GCI 0.13 %) the clean reconstruction error is 0.5–0.75 %; mesh-block kinks in the data are amplified like noise." },
  ];
  const right = [
    { k: "head", t: "Objectives" },
    { k: "text", t: [B("All four met: "), R("algorithms implemented, verified on exact flows, studied parametrically, and compared with a CFD wake and published studies.")] },
    { k: "head", t: "Recommended practice" },
    { k: "text", t: "Poisson solver; Dirichlet edges only where p is known (free stream); Lagrangian Du/Dt with the largest Δt that keeps truncation small; spacing near the noise–truncation optimum; light denoising above σ ≈ 0.5 %." },
    { k: "head", t: "Main takeaway" },
    { k: "text", t: "PIV pressure is limited by how velocity noise is amplified by differentiation. Method, boundary data, spacing and Δt decide whether that noise is averaged out or accumulated." },
    { k: "note", t: "Way ahead: make noise depend on W (larger windows, less noise in real PIV); correlated noise; 3-D time-resolved data." },
  ];
  if (wide) {
    para(s, left, { x: G.X0, y: G.Y0 + 0.1, w: 5.9, h: 5.9 }, { bgap: 6, size: 13.5, headSize: 15 });
    para(s, right, { x: 6.85, y: G.Y0 + 0.1, w: 5.88, h: 5.9 }, { gap: 12, size: 14, headSize: 15 });
  } else {
    para(s, left, { x: 0.45, y: G.Y0, w: 4.45, h: 5.4 }, { size: 11.5, headSize: 14, bgap: 4, align: "left" });
    para(s, right, { x: 5.15, y: G.Y0, w: 4.45, h: 5.4 }, { size: 12.5, headSize: 14, gap: 7, align: "left" });
  }
  s.addNotes("Summarise in one sentence: the accuracy of PIV pressure is set by how measurement noise propagates, and our study shows which choices control that.");
}
};

BL.refs = () => {
{
  const s = wide ? newSlide("References") : newSlide("References", "Papers cited in this talk");
  const refs = [
    "B. W. van Oudheusden, PIV-based pressure measurement, Meas. Sci. Technol. 24 (2013) 032001.",
    "J. J. Charonko, C. V. King, B. L. Smith, P. P. Vlachos, Assessment of pressure field calculations from particle image velocimetry measurements, Meas. Sci. Technol. 21 (2010) 105401.",
    "Z. Pan, J. Whitehead, S. Thomson, T. Truscott, Error propagation dynamics of PIV-based pressure field calculations: how well does the pressure Poisson solver perform inherently?, Meas. Sci. Technol. 27 (2016) 084012.",
    "P. L. van Gent et al., Comparative assessment of pressure field reconstructions from particle image velocimetry measurements and Lagrangian particle tracking, Exp. Fluids 58 (2017).",
    "X. Liu, J. Katz, Instantaneous pressure and material acceleration measurements using a four-exposure PIV system, Exp. Fluids 41 (2006) 227–240.",
    "V. Violato, P. Moore, F. Scarano, Lagrangian and Eulerian pressure field evaluation of rod-airfoil flow from time-resolved tomographic PIV, Exp. Fluids 50 (2011) 1057–1070.",
    "J. Park, K. Kwon, H. Choi, Numerical solutions of flow past a circular cylinder at Reynolds numbers up to 160, KSME Int. J. 12 (1998) 1200–1205.",
    "C. H. K. Williamson, Vortex dynamics in the cylinder wake, Annu. Rev. Fluid Mech. 28 (1996) 477–539.",
    "L. I. G. Kovasznay, Laminar flow behind a two-dimensional grid, Proc. Camb. Philos. Soc. 44 (1948) 58–62.",
  ];
  const runs = refs.map((r, i) => ({ text: `[${i + 1}]  ${r}`, options: { fontFace: BODY, fontSize: wide ? 13 : 11.5, color: i < 6 ? C.body : (wide ? "000000" : C.accent), paraSpaceBefore: 7, breakLine: i < refs.length - 1, indentLevel: 0 } }));
  s.addText(runs, { x: G.X0, y: G.Y0 + 0.1, w: G.CW, h: G.Y1 - G.Y0 - 0.1, margin: 0, valign: "top", isTextBox: true });
  s.addNotes("References. The cylinder values quoted from Park et al. are approximate; check them in the paper before quoting exact numbers.");
}
};

BL.backupindex = () => {
{
  const s = wide ? newSlide("Backup Slides") : newSlide("Backup Slides", "Supporting detail for questions");
  if (wide) s.addText("Supporting detail for questions", { x: 0, y: 1.6, w: G.W, h: 0.4, fontFace: BODY, fontSize: 16, color: "404040", align: "center", margin: 0, isTextBox: true });
  const list = SHORT_FMT ? ["CFD set-up, grid study, wake Δt and noise", "Noise sensitivity, window and denoising", "Iterative solvers and scheme order", "Kovasznay reconstruction", "CFD lift, spectrum and residuals"] : ["Iterative Poisson solvers: Jacobi, Gauss–Seidel, SOR, CG", "First- vs second-order schemes", "Kovasznay reconstruction and error map", "CFD lift history, spectrum and solver residuals"];
  para(s, list.map((t) => ({ k: "b", t })), { x: wide ? 4.0 : 2.0, y: wide ? 2.4 : 2.6, w: 7.5, h: 2.0 }, { size: wide ? 16 : 17, bgap: 8, align: "left" });
}
};

BL.b_solvers = () => {
{
  const s = wide ? newSlide("Backup: Iterative Poisson Solvers") : newSlide("Iterative Solvers", "Backup: same Poisson system, five solvers");
  const rows = [["Solver (64², tol 1e-8)", "Iterations", "Time"], ["Jacobi (ω = 0.95)", "25,366", "1.3 s"], ["Gauss–Seidel", "12,328", "0.93 s"], ["SOR (ω = 1.932)", "314", "0.024 s"], ["Conjugate gradient", "180", "0.013 s"], ["Jacobi, classic (ω = 1)", "stalls at 4.2e-8", "–"]];
  const items = [
    { k: "b", t: "All solve the same matrix as the direct solve; SOR at ω_opt and CG need 80–140× fewer iterations than Jacobi." },
    { k: "b", t: "Iterations grow as N^1.0 (Jacobi, GS), N^0.52 (SOR) and N^0.43 (CG)." },
    { k: "b", t: "Classic Jacobi stalls on the all-Neumann problem (checkerboard mode); ω = 0.95 fixes it. SOR does not converge at ω = 2 and diverges at 2.1." },
  ];
  if (wide) {
    fig(s, "b_residuals", { x: G.X0, y: G.Y0, w: G.CW, h: 3.6 }, "Scaled residual history, Taylor–Green 64²", { border: false });
    table(s, rows, { x: G.X0, y: 4.95, w: 5.6, colW: [2.6, 1.8, 1.2], size: 11, rowH: 0.27 });
    para(s, items, { x: 6.55, y: 4.95, w: 6.18, h: 2.1 }, { size: 12.5, bgap: 4 });
  } else {
    fig(s, "b_residuals", { x: 0.4, y: G.Y0, w: 9.2, h: 3.15 }, "Scaled residual history, Taylor–Green 64²", { border: false, capSize: 10 });
    table(s, rows, { x: 0.4, y: 5.0, w: 4.4, colW: [2.0, 1.45, 0.95], size: 9.5, rowH: 0.25 });
    para(s, items, { x: 5.0, y: 4.95, w: 4.6, h: 2.2 }, { size: 11, bgap: 3, align: "left" });
  }
}
};

BL.b_order = () => {
{
  const s = wide ? newSlide("Backup: First- vs. Second-Order Schemes") : newSlide("Scheme Order", "Backup: first- vs second-order differences");
  const items = [
    { k: "b", t: "Fitted slopes: 0.98 vs 1.96 (Taylor–Green), 1.00 vs 1.96 (Kovasznay Re = 20)." },
    { k: "b", t: "1 % error is reached at Δx ≈ 0.153 with second order but needs Δx ≈ 0.012 with first order: a 12× finer grid." },
    { k: "b", t: "With 1 % noise both schemes level off and rise again on fine grids; the second-order scheme stays lower." },
  ];
  if (wide) {
    fig(s, "b_order", { x: G.X0, y: G.Y0, w: 7.6, h: 5.6 }, "Grid convergence of the reconstructed pressure, Taylor–Green", { border: false });
    para(s, items, { x: 8.5, y: G.Y0 + 0.4, w: G.X1 - 8.5, h: 4.0 }, { bgap: 8 });
  } else {
    fig(s, "b_order", { x: 0.4, y: G.Y0, w: 9.2, h: 4.0 }, "Grid convergence, Taylor–Green", { border: false, capSize: 10 });
    para(s, items, { x: 0.45, y: 5.85, w: 9.1, h: 1.3 }, { size: 11.5, bgap: 3 });
  }
}
};

BL.b_kov = () => {
{
  const s = wide ? newSlide("Backup: Kovasznay Reconstruction") : newSlide("Kovasznay Flow", "Backup: noisy reconstruction and error map");
  if (wide) {
    fig(s, "f_fields_kov", { x: G.X0, y: G.Y0 + 0.2, w: G.CW, h: 4.0 }, "Kovasznay Re = 20: exact, reconstructed (σ = 1 %, 64², Poisson) and error", { border: false });
    para(s, [{ k: "b", t: "ε_L2 = 0.056 for this seed. The error is mostly random scatter from the noise; the largest values sit at the edges and corners, where one-sided stencils and the noisy Neumann data n·G act." }], { x: G.X0, y: 5.6, w: G.CW, h: 0.8 });
  } else {
    fig(s, "f_fields_kov", { x: 0.4, y: G.Y0 + 0.2, w: 9.2, h: 3.2 }, "Kovasznay Re = 20: exact, reconstructed (σ = 1 %) and error", { border: false, capSize: 10 });
    para(s, [{ k: "b", t: "ε_L2 = 0.056 for this seed. The error is mostly random scatter from the noise; the largest values sit at the edges and corners, where one-sided stencils and the noisy Neumann data n·G act." }], { x: 0.45, y: 5.4, w: 9.1, h: 0.8 }, { size: 12.5 });
  }
}
};

BL.b_cfdlift = () => {
{
  const s = wide ? newSlide("Backup: CFD Lift, Spectrum and Residuals") : newSlide("CFD Lift & Residuals", "Backup: forces, spectrum and solver residuals");
  const items = [
    { k: "b", t: "Cl from all three meshes over the last 30 time units; the FFT peak gives St (Hann window, uniform resampling)." },
    { k: "b", t: "Coarse mesh: lower Cl amplitude (0.322 vs 0.345) and lower St (0.1644): too little resolution at the cylinder." },
    { k: "b", t: "Initial residuals fall from 1, then level off near 10⁻³ (p, Uy) and 3 × 10⁻⁴ (Ux): each time step starts from a changed, periodic flow, so they never reach zero." },
  ];
  if (wide) {
    fig(s, "c_clfft", { x: G.X0, y: G.Y0, w: G.CW, h: 3.6 }, "Lift coefficient and its spectrum, three meshes", { border: false });
    fig(s, "c_resid_fine", { x: G.X0, y: 4.9, w: 3.4, h: 2.25 }, "Fine-mesh initial residuals", { border: false, capH: 0.25, capSize: 11 });
    para(s, items, { x: 4.4, y: 4.95, w: G.X1 - 4.4, h: 2.2 }, { size: 12.5, bgap: 5 });
  } else {
    fig(s, "c_clfft", { x: 0.4, y: G.Y0, w: 9.2, h: 3.0 }, "Lift coefficient and its spectrum, three meshes", { border: false, capSize: 10 });
    fig(s, "c_resid_fine", { x: 0.4, y: 4.85, w: 3.2, h: 2.25 }, "Fine-mesh residuals", { border: false, capSize: 10, capH: 0.25 });
    para(s, items, { x: 3.8, y: 4.9, w: 5.8, h: 2.2 }, { size: 11, bgap: 4, align: "left" });
  }
}
};

BL.fin = () => {
if (!wide) {
  const s = pres.addSlide();
  s.background = { path: path.join(FIG, "gridpaper.jpg") };
  s.addText("Fin.", { x: 0, y: 2.75, w: 10, h: 1.3, fontFace: TITLEF, fontSize: 88, color: "45A7DE", align: "center", valign: "middle", margin: 0, isTextBox: true });
  s.addText("Thank you. Questions?", { x: 0, y: 4.05, w: 10, h: 0.55, fontFace: TITLEF, fontSize: 28, color: "8E5A0C", align: "center", margin: 0, isTextBox: true });
}
};

// SHORT-deck slides are defined in short_slides.js
// ================================================================== 10-MINUTE (SHORT) DECK — lecture format only
const SCRIPT = {
  cover: "[ANANY · 0:00–0:15 · 15 s]\nGood morning. I'm Anany, and this is Arnav. Our project is pressure from PIV: we rebuild the pressure field from velocity measurements, compare three algorithms, and test them against exact solutions and a CFD simulation.",
  s_problem: "[ANANY · 0:15–1:05 · 50 s]\nPIV photographs tracer particles in a laser sheet twice, a short time apart, and gives the velocity over a whole plane without touching the flow. But it does not measure pressure, which sets forces, noise and cavitation. So pressure has to be inferred: differentiate the measured velocity, then integrate, and differentiation amplifies measurement noise.\nOur four objectives: implement three integration algorithms, verify them on exact Navier–Stokes solutions, find what controls the error, and compare with an OpenFOAM cylinder wake and published studies.",
  equations: "[ANANY · 1:05–2:00 · 55 s]\nWe start from incompressible Navier–Stokes, equation 1, and rearrange momentum for the pressure gradient, equation 2: G equals minus rho times the material acceleration plus the viscous term. Every term is a derivative of measured velocity.\nThe material derivative needs three snapshots: Eulerian, equation 3, at a fixed point, or Lagrangian, equation 4, following the particle.\nTaking the divergence of G gives the pressure Poisson equation, 5, with Neumann boundary data. It is solvable only if the compatibility condition, 6, holds, and only up to a constant, so we fix the mean to zero. Physically it is a least-squares fit that keeps the gradient part of G and discards the curl, which is pure error. We assume incompressible, two-dimensional flow.",
  algorithms: "[ANANY · 2:00–2:50 · 50 s]\nThree ways to integrate G. Direct integration marches G along a path from a corner with the trapezoidal rule, x-then-y and y-then-x, and averages the two. It is the fastest, but every error is carried along the path.\nOmni-directional integration, after Liu and Katz, averages line integrals from 32 directions into each node. It is robust but slow.\nThe Poisson solver is our main method: a sparse 5-point Laplacian, Neumann edges by ghost-point elimination, and a zero-mean constraint for the singular system, solved directly in about 15 milliseconds.",
  strategy: "[ANANY · 2:50–3:30 · 40 s]\nOur logic chart: take a velocity source, degrade it like a real PIV system, optionally denoise, compute the material derivative, form G, integrate, and compare with the exact pressure. Each study varies one parameter around a baseline, with ten noise seeds.\nThe code is Python, with 87 unit tests that check every stencil and solver against exact answers.",
  testcases: "[ANANY · 3:30–4:05 · 35 s]\nWe need flows whose pressure is known exactly. Taylor–Green is unsteady and its viscous term cancels, so it tests convection; a convected copy tests the time derivative. Kovasznay flow is steady with a strong viscous term. And the CFD wake.\nThe synthetic PIV model box-filters the exact field over the window, keeps every fourth point to get spacing h, and adds 1 percent noise. The error is the relative L2 norm of the demeaned pressure difference.",
  s_accuracy: "[ANANY · 4:05–5:00 · 55 s]\nOn clean data, top, all three methods converge at second order on all three flows, with slopes 1.8 to 2.1, so the code is verified.\nWith 1 percent noise, bottom, the error stops falling: Poisson has a minimum near h = 0.03 and rises again on finer grids. Truncation error falls as h squared, but noise is amplified as one over h by the convective term and one over h squared by the viscous term. So noisy PIV has an optimum vector spacing, and the finest grid can make pressure worse.\nDirect integration only gets worse: its errors add up like a random walk.",
  s_methodbc: "[ANANY · 5:00–5:50 · 50 s]\nAt 1 percent noise, left, direct integration is five to seven times worse than Poisson or omni: one path accumulates every error, while Poisson and omni average over many paths. Omni and Poisson are equally accurate, but omni is about 27 times slower, so Poisson is the practical choice.\nRight, boundary conditions: true pressure on two opposite edges halves the error, because Neumann data are themselves noisy derivatives. Bernoulli edges fail, up to 0.8, because Bernoulli needs irrotational flow and both test flows carry vorticity.\nI'll hand over to Arnav.",
  s_dt: "[ARNAV · 5:50–6:40 · 50 s]\nThanks. Theory says the Lagrangian derivative wins when structures are carried along faster than they spin. On the convected vortex, left, its error is three times lower at large time step. With noise, middle, both schemes improve as the time step grows, because noise in the time derivative scales as sigma over delta t.\nOn the CFD wake, right, Lagrangian wins only at the smallest step. Our window sits where vortices are still forming and spin about as fast as they move, so the assumption breaks.",
  s_cfd: "[ARNAV · 6:40–7:30 · 50 s]\nWe simulated the laminar cylinder wake at Reynolds number 100 in OpenFOAM, on three meshes up to 107,000 cells. This is the fine-mesh pressure in ParaView, with the Kármán vortex street.\nThe Strouhal number converges to 0.1695, with observed order 2.5 and a grid convergence index of 0.13 percent. Drag and lift converge on the medium mesh, and shedding is periodic from t = 76, well before our sampling window. Our values are 3 to 4 percent above the literature, plausibly from the 5 percent blockage of the side walls.",
  cfdrecon: "[ARNAV · 7:30–8:20 · 50 s]\nFrom an 81 by 81 PIV-like sample of the CFD velocity, the pressure is recovered to within 0.5 to 0.75 percent, and every method matches the vortex-core suction within about 1 percent.\nThe remaining error sits on straight lines. These are the block edges of the CFD mesh, where the sampled velocity has small kinks that differentiation amplifies, just like noise. With 1 percent noise the wake error rises to about 14 percent, because the wake's pressure variations are small compared with the velocity noise.",
  lit: "[ARNAV · 8:20–8:55 · 35 s]\nOur trends reproduce the literature: Pan and co-workers on Dirichlet versus Neumann boundaries and on Bernoulli edges, Charonko and co-workers on errors from below 1 to over 100 percent, and Violato and co-workers on the Eulerian time-step limit. The CFD wake matches published Strouhal and drag values within the blockage effect.",
  concl: "[ARNAV · 8:55–9:50 · 55 s]\nTo conclude. The code is verified at second order. Noise, not discretisation, limits PIV pressure: the error grows linearly with noise and there is an optimum vector spacing. Poisson integration is robust and fast, known pressure on two edges halves the error, and the Lagrangian derivative helps only where structures are simply convected. All four objectives are met.\nThe main takeaway: the accuracy of pressure from PIV is set by how velocity noise is amplified by differentiation, and the choice of method, boundary data, spacing and time step decides whether that noise is averaged out or accumulated. Thank you.",
  fin: "[BOTH · 9:50–10:00]\nThank you. We're happy to take questions. Backup slides cover the CFD set-up and grid study, noise and window effects, the solvers, scheme order and CFD residuals.",
};

BL.s_problem = () => {
  const s = newSlide("Problem & Objectives", "PIV gives velocity, not pressure");
  fig(s, "f_piv_schematic", { x: 0.4, y: G.Y0, w: 4.9, h: 2.95 }, "Planar PIV: two frames Δt apart → velocity field", { border: false, capSize: 10, capH: 0.25 });
  para(s, [
    { k: "b", t: "PIV measures the whole velocity field without touching the flow, but not pressure." },
    { k: "b", t: "Momentum gives the pressure gradient from velocity alone:" },
  ], { x: 5.5, y: G.Y0 + 0.05, w: 4.1, h: 1.35 }, { size: 12, bgap: 4, align: "left" });
  eq(s, "eq_G", 5.6, 3.15, 0.5);
  para(s, [{ k: "b", t: [B("Differentiate, then integrate: "), R("differentiation amplifies measurement noise.")] }], { x: 5.5, y: 3.85, w: 4.1, h: 0.75 }, { size: 12, align: "left" });
  para(s, [
    { k: "head", t: "Objectives" },
    { k: "num", t: "Implement three integration algorithms: direct, omni-directional, Poisson." },
    { k: "num", t: "Verify them on exact Navier–Stokes solutions." },
    { k: "num", t: "Find what controls the error: spacing, noise, Δt, boundary conditions." },
    { k: "num", t: "Compare with an OpenFOAM cylinder wake and published studies." },
  ], { x: 0.45, y: 4.85, w: 9.1, h: 2.25 }, { size: 12.5, headSize: 14, bgap: 2, align: "left" });
};

BL.s_method = () => {
  const s = newSlide("Method", "Pressure-Poisson and two alternatives");
  eq(s, "eq_poisson", 0.45, G.Y0 + 0.1, 0.6);
  para(s, [{ k: "text", t: "Solvable only if source and boundary flux balance; fixed up to a constant → zero mean. A least-squares fit: keeps the gradient part of G, discards the curl (pure error)." }], { x: 5.6, y: G.Y0, w: 4.0, h: 1.0 }, { size: 11, align: "left" });
  fig(s, "f_methods", { x: 0.4, y: 2.8, w: 9.2, h: 2.75 }, "(a) direct along two paths, (b) rays from 32 directions, (c) 5-point Poisson stencil with a ghost node", { border: false, capSize: 10, capH: 0.25 });
  const cw = 2.95;
  [
    [{ k: "head", t: "Direct" }, { k: "text", t: "March G·dl from a corner. Fastest (0.3 ms); errors carried forward." }],
    [{ k: "head", t: "Omni-directional" }, { k: "text", t: "Average of 32 ray integrals. Robust; slow (0.4 s)." }],
    [{ k: "head", t: "Poisson (main)" }, { k: "text", t: "Sparse 5-point Laplacian, ghost-point Neumann edges, zero-mean. 15 ms." }],
  ].forEach((c, i) => para(s, c, { x: 0.45 + i * (cw + 0.12), y: 5.65, w: cw, h: 1.45 }, { size: 11, headSize: 12.5, align: "left" }));
};

BL.s_accuracy = () => {
  const s = newSlide("Accuracy vs Grid Spacing", "Second order when clean; an optimum spacing when noisy");
  fig(s, "f_conv_clean", { x: 0.4, y: G.Y0, w: 5.95, h: 2.65 }, "Clean data: slopes 1.8–2.1 (second order)", { border: false, capSize: 10, capH: 0.25 });
  fig(s, "f_conv_noisy", { x: 0.4, y: 4.45, w: 5.95, h: 2.65 }, "σ = 1 % noise: Poisson optimum at h ≈ 0.03 (★)", { border: false, capSize: 10, capH: 0.25 });
  para(s, [
    { k: "num", t: [B("Verified: "), R("all methods, all flows converge at second order on clean data.")] },
    { k: "num", t: [B("With noise the error rises again "), R("below h ≈ 0.03 (TG: 0.082 → 0.20).")] },
    { k: "num", t: [B("Why: "), R("truncation ∝ h², but noise is amplified ∝ σ/h (convective) and μσ/h² (viscous).")] },
    { k: "num", t: [B("Direct "), R("only gets worse: errors add like a random walk along the path.")] },
  ], { x: 6.5, y: G.Y0 + 0.1, w: 3.1, h: 5.3 }, { size: 11.5, bgap: 8, align: "left" });
};

BL.s_methodbc = () => {
  const s = newSlide("Method & Boundary Conditions", "Errors at σ = 1 % noise, 10 seeds");
  fig(s, "f_method_bar", { x: 0.4, y: G.Y0, w: 4.55, h: 3.15 }, "Integration methods", { border: false, capSize: 10, capH: 0.25 });
  fig(s, "f_bc_bar", { x: 5.05, y: G.Y0, w: 4.55, h: 3.15 }, "Poisson boundary conditions", { border: false, capSize: 10, capH: 0.25 });
  para(s, [
    { k: "b", t: [B("Direct is 5–7× worse: "), R("one path accumulates every error; Poisson and omni average over many paths.")] },
    { k: "b", t: [B("Omni ≈ Poisson in accuracy, "), R("but ≈27× slower → Poisson is the practical choice.")] },
    { k: "b", t: [B("Known p on 2 opposite edges halves the error "), R("(0.082 → 0.044): Neumann data are themselves noisy derivatives.")] },
    { k: "b", t: [B("Bernoulli edges fail (up to 0.81): "), R("valid only in irrotational flow; both test flows carry vorticity.")] },
  ], { x: 0.45, y: 5.0, w: 9.1, h: 2.1 }, { size: 12, bgap: 3, align: "left" });
};

BL.s_dt = () => {
  const s = newSlide("Time Step: Eulerian vs Lagrangian", "Following the particle helps only when structures are convected");
  fig(s, "f_dt", { x: 0.4, y: G.Y0, w: 6.0, h: 3.0 }, "Convected vortex (clean) and stationary vortex (σ = 1 %)", { border: false, capSize: 10, capH: 0.25 });
  fig(s, "c_dt_cfd_clean", { x: 6.55, y: G.Y0, w: 3.05, h: 3.0 }, "CFD wake (clean)", { border: false, capSize: 10, capH: 0.25 });
  para(s, [
    { k: "num", t: [B("Convected vortex: "), R("Lagrangian error 3× lower at Δt = 0.2: the vortex is carried faster than it spins.")] },
    { k: "num", t: [B("With noise "), R("both improve with larger Δt: noise in ∂u/∂t scales as σ/Δt.")] },
    { k: "num", t: [B("CFD wake: "), R("Lagrangian wins only at Δt = 0.05 (15–18 %). In the vortex-formation region vortices spin about as fast as they move, so the assumption breaks.")] },
  ], { x: 0.45, y: 5.0, w: 9.1, h: 2.1 }, { size: 12, bgap: 4, align: "left" });
};

BL.s_cfd = () => {
  const s = newSlide("CFD: Cylinder Wake at Re = 100", "OpenFOAM, three meshes, verified");
  fig(s, "c_paraview", { x: 0.4, y: G.Y0, w: 5.0, h: 3.65 }, "ParaView 6.2: fine-mesh pressure, t = 200", { border: false, capSize: 10, capH: 0.25 });
  const gRows = [["Mesh", "Cells", "St", "C̄d"],
    ["Coarse", "6,690", "0.1644", "1.375"],
    ["Medium", "26,760", "0.1688", "1.388"],
    [{ text: "Fine", options: { bold: true } }, { text: "107,040", options: { bold: true } }, { text: "0.1695", options: { bold: true } }, { text: "1.388", options: { bold: true } }],
    ["Richardson", "–", "0.1697", "conv."],
    [{ text: "Literature", options: { italic: true } }, "–", "≈0.165", "≈1.33"]];
  table(s, gRows, { x: 5.6, y: G.Y0, w: 4.0, colW: [1.25, 0.95, 0.9, 0.9], size: 10, rowH: 0.27 });
  para(s, [
    { k: "b", t: "pimpleFoam (PISO), 2nd-order schemes, Co ≤ 0.5, t = 200." },
    { k: "b", t: "St: order p = 2.48, GCI 0.13 %; periodic from t = 76." },
    { k: "b", t: "3–4 % above literature: 5 % side-wall blockage." },
  ], { x: 5.6, y: 3.55, w: 4.0, h: 1.9 }, { size: 11, bgap: 4, align: "left" });
  para(s, [{ k: "text", t: [B("Window for PIV: "), R("81 × 81 points, x/D 1–5, y/D −2 to 2, every 0.05 D/U for t = 175–200 (about 4 shedding periods), from the fine mesh.")] }], { x: 0.45, y: 5.65, w: 9.1, h: 0.9 }, { size: 11.5, align: "left" });
};

BL.__shortOrder = ["cover", "s_problem", "equations", "algorithms", "strategy", "testcases", "s_accuracy", "s_methodbc", "s_dt", "s_cfd", "cfdrecon", "lit", "concl", "fin",
  "refs", "backupindex", "cfdsetup", "cfdgrid", "cfddt", "noisemethod", "window", "b_solvers", "b_order", "b_kov", "b_cfdlift"];

const ORDER = SHORT ? BL.__shortOrder : ['cover','intro','problem','equations','algorithms','strategy','testcases','verif','noiseres','noisemethod','bc','dt','window','cfdsetup','cfdgrid','cfdrecon','cfddt','lit','concl','refs','backupindex','b_solvers','b_order','b_kov','b_cfdlift','fin'];
ORDER.forEach((k) => {
  if (!BL[k]) throw new Error("no slide " + k);
  const n0 = ALLS.length; BL[k]();
  ALLS.slice(n0).forEach((sl, i) => { const t = (SHORT && i === 0 && SCRIPT[k]) ? SCRIPT[k] : sl.__notes; if (t) sl.__orig(t); });
});

// page numbers (lecture style: "n /N")
if (!wide) {
  const N = ALLS.length;
  lecSlides.forEach((s) => {
    const pos = ALLS.indexOf(s) + 1;
    s.addText([{ text: String(pos) + " ", options: { color: "8A8A8A", fontSize: 11 } }, { text: "/" + N, options: { color: "000000", fontSize: 14 } }], { x: 8.9, y: 7.0, w: 1.05, h: 0.28, fontFace: TITLEF, align: "right", valign: "bottom", margin: 0, isTextBox: true });
  });
}

pres.writeFile({ fileName: OUTFILE }).then(() => console.log("wrote", OUTFILE, "figs", figNo, "eqs", eqNo));
