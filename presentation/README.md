# Presentation material

| File | What it is |
|---|---|
| `PIV_Pressure_ME415_Lecture-format_10min.pptx` | 10-minute talk (14 slides + backups), lecture format; speaker script in the notes |
| `PIV_Pressure_ME415_Lecture-format.pptx` | Full deck, lecture format |
| `PIV_Pressure_ME415_EDM-format.pptx` | Full deck, ME 323 report format |
| `PIV_Pressure_Complete_Explainer.pdf` | 59-page explainer: theory, numerics, code, every result, CFD, viva Q&A |

## Rebuilding

- **Figures** (from the project's CSVs): `python source/figs.py <path-to-repo-root> source/figs` and
  `python source/figs2.py <path-to-repo-root> source/figs`.
- **Decks** (Node.js + `pptxgenjs`), from `source/`:
  - `node build.js lec10 out.pptx` (10-minute lecture deck)
  - `node build.js lec out.pptx` (full lecture deck)
  - `node build.js edm out.pptx` (EDM-format deck)
- **Explainer PDF**: `cd source/report && pdflatex main.tex` (run twice for the table of contents).

`source/figs/c_*.png` are the OpenFOAM results and the ParaView screenshot from the CFD run
(also in `results/cfd/`).
