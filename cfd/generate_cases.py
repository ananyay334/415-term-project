# piv_pressure/cfd/generate_cases.py
"""Generate the three OpenFOAM cases (coarse / medium / fine) for the
2-D laminar cylinder wake at Re = 100.

Target: OpenFOAM (ESI / openfoam.com) v2306 or newer, solver ``pimpleFoam``.
(Foundation OpenFOAM >= 11 uses ``foamRun -solver incompressibleFluid`` and
``constant/physicalProperties`` + ``constant/momentumTransport``; see
README [verify for your version].)

Non-dimensional set-up: D = 1, U_inf = 1, nu = 0.01 (Re = 100), rho = 1
(OpenFOAM's incompressible p is the kinematic pressure p/rho).

Geometry and blocking (12 hex blocks, one cell thick in z)
---------------------------------------------------------
Domain x in [-10, 20], y in [-10, 10]; cylinder of radius 0.5 at the origin.
An O-grid of four blocks joins the cylinder (arcs) to the square
[-a, a]^2 (a = 2); eight rectangular blocks fill the rest:

      13 ---- 12 -------- 11 ------------------ 10      y = 10
      |  TL   |    TM     |          TR          |
      14 ---- 6 --------- 5 ------------------- 9       y = a
      |  LM   |  O-grid   |          RM          |
      15 ---- 7 --------- 4 ------------------- 8       y = -a
      |  BL   |    BM     |          BR          |
      16 ---- 17 -------- 18 ------------------ 19      y = -10
    x = -10   -a          a                    x = 20

Point indices are those of the z-min plane; the z-max copy is index + 20.
Running this script rewrites ``cases/{coarse,medium,fine}`` and checks the
block topology (orientation, matching cell counts and grading on shared
edges, every outer edge assigned to a patch) before writing anything.

Usage:  python cfd/generate_cases.py
"""
from __future__ import annotations

import json
import math
import stat
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CASES = HERE / "cases"

# ----------------------------------------------------------------- geometry
R = 0.5                      # cylinder radius (D = 1)
A = 2.0                      # half-width of the O-grid square
X_IN, X_OUT, H = -10.0, 20.0, 10.0
Z0, Z1 = -0.05, 0.05         # one cell thick in z; span 0.1 -> Aref = D * 0.1
NU, U_INF = 0.01, 1.0
C45 = R / math.sqrt(2.0)

# PIV window exported from the fine mesh (x/D in [1, 5], y/D in [-2, 2], 0.05 D).
WINDOW = {"x0": 1.0, "x1": 5.0, "y0": -2.0, "y1": 2.0, "nx": 81, "ny": 81}
SAMPLE_START, END_TIME, SAMPLE_DT = 175.0, 200.0, 0.05


@dataclass
class MeshLevel:
    """Cell counts of one mesh level (see module docstring for the blocks)."""

    name: str
    n_t: int        # cells per O-grid quarter (tangential) = cells across [-a, a]
    n_r: int        # radial cells in the O-grid
    n_up: int       # cells upstream (x in [-10, -a])
    n_down: int     # cells downstream (x in [a, 20])
    n_side: int     # cells to the side walls (|y| in [a, 10])
    e_r: float = 8.0      # O-grid radial expansion (outer / inner cell)
    e_up: float = 8.0     # upstream expansion (far / near cell)
    e_down: float = 3.0   # downstream expansion (far / near cell)
    e_side: float = 8.0   # lateral expansion (far / near cell)
    sample: bool = False  # add the PIV-window sampling function object
    nu: float = NU                      # kinematic viscosity (Re = U D / nu)
    end_time: float = END_TIME          # upper bound of the simulated time
    sample_start: float = SAMPLE_START  # first sampling time of the window
    sample_dt: float = SAMPLE_DT        # sampling interval of the window
    kick: bool = True                   # transverse inlet kick to start shedding
    residual_stop: float | None = None  # stop when all initial residuals < this

    def n_cells(self) -> int:
        """Total cell count."""
        o = 4 * self.n_t * self.n_r
        mid = self.n_t * (self.n_up + self.n_down + 2 * self.n_side)
        corners = 2 * self.n_side * (self.n_up + self.n_down)
        return o + mid + corners


LEVELS = [
    MeshLevel("coarse", n_t=24, n_r=20, n_up=15, n_down=60, n_side=15),
    MeshLevel("medium", n_t=48, n_r=40, n_up=30, n_down=120, n_side=30),
    MeshLevel("fine", n_t=96, n_r=80, n_up=60, n_down=240, n_side=60, sample=True),
]
#: Optional steady case (Re = 20 is below the shedding threshold, so the flow
#: becomes steady): medium mesh, no kick, window sampled every 0.5 D/U from
#: t = 0 and an automatic stop when every initial residual is below 1e-5.
OPTIONAL_LEVELS = [
    MeshLevel("steady_re20", n_t=48, n_r=40, n_up=30, n_down=120, n_side=30, sample=True,
              nu=0.05, end_time=300.0, sample_start=0.0, sample_dt=0.5, kick=False,
              residual_stop=1e-5),
]

POINTS_2D = {
    0: (C45, -C45), 1: (C45, C45), 2: (-C45, C45), 3: (-C45, -C45),
    4: (A, -A), 5: (A, A), 6: (-A, A), 7: (-A, -A),
    8: (X_OUT, -A), 9: (X_OUT, A), 10: (X_OUT, H), 11: (A, H), 12: (-A, H),
    13: (X_IN, H), 14: (X_IN, A), 15: (X_IN, -A), 16: (X_IN, -H), 17: (-A, -H),
    18: (A, -H), 19: (X_OUT, -H),
}
NP = len(POINTS_2D)

# Cylinder arcs: (start, end, mid-point on the circle)
ARCS = [(0, 1, (R, 0.0)), (1, 2, (0.0, R)), (2, 3, (-R, 0.0)), (3, 0, (0.0, -R))]

# Outer patches: 2-D edges traversed with the fluid on the LEFT (outward normal).
PATCH_EDGES = {
    "inlet": ("patch", [(13, 14), (14, 15), (15, 16)]),
    "outlet": ("patch", [(19, 8), (8, 9), (9, 10)]),
    "top": ("symmetryPlane", [(10, 11), (11, 12), (12, 13)]),
    "bottom": ("symmetryPlane", [(16, 17), (17, 18), (18, 19)]),
    "cylinder": ("wall", [(1, 0), (2, 1), (3, 2), (0, 3)]),
}


def blocks(m: MeshLevel) -> list[tuple[str, tuple[int, int, int, int], tuple[int, int], tuple[float, float]]]:
    """(name, base quad v0..v3 counter-clockwise, (n1, n2), (grading1, grading2))."""
    return [
        ("O_right", (0, 4, 5, 1), (m.n_r, m.n_t), (m.e_r, 1.0)),
        ("O_top", (1, 5, 6, 2), (m.n_r, m.n_t), (m.e_r, 1.0)),
        ("O_left", (2, 6, 7, 3), (m.n_r, m.n_t), (m.e_r, 1.0)),
        ("O_bottom", (3, 7, 4, 0), (m.n_r, m.n_t), (m.e_r, 1.0)),
        ("RM", (4, 8, 9, 5), (m.n_down, m.n_t), (m.e_down, 1.0)),
        ("TR", (5, 9, 10, 11), (m.n_down, m.n_side), (m.e_down, m.e_side)),
        ("TM", (6, 5, 11, 12), (m.n_t, m.n_side), (1.0, m.e_side)),
        ("TL", (14, 6, 12, 13), (m.n_up, m.n_side), (1.0 / m.e_up, m.e_side)),
        ("LM", (15, 7, 6, 14), (m.n_up, m.n_t), (1.0 / m.e_up, 1.0)),
        ("BL", (16, 17, 7, 15), (m.n_up, m.n_side), (1.0 / m.e_up, 1.0 / m.e_side)),
        ("BM", (17, 18, 4, 7), (m.n_t, m.n_side), (1.0, 1.0 / m.e_side)),
        ("BR", (18, 19, 8, 4), (m.n_down, m.n_side), (m.e_down, 1.0 / m.e_side)),
    ]


# ------------------------------------------------------------------ checks
def check_topology(m: MeshLevel) -> None:
    """Assert orientation, edge matching and patch coverage of the blocking."""
    edge_info: dict[tuple[int, int], list[tuple[int, float, str]]] = {}
    for name, quad, (n1, n2), (g1, g2) in blocks(m):
        pts = np.array([POINTS_2D[k] for k in quad])
        area = 0.5 * np.sum(pts[:, 0] * np.roll(pts[:, 1], -1) - np.roll(pts[:, 0], -1) * pts[:, 1])
        assert area > 0, f"block {name} is not counter-clockwise"
        v0, v1, v2, v3 = quad
        # edges in local directions: x1 along v0->v1 and v3->v2, x2 along v0->v3 and v1->v2
        for (a, b), n, g in (((v0, v1), n1, g1), ((v3, v2), n1, g1),
                             ((v0, v3), n2, g2), ((v1, v2), n2, g2)):
            key = (min(a, b), max(a, b))
            g_norm = g if a < b else 1.0 / g          # grading expressed from low to high index
            edge_info.setdefault(key, []).append((n, g_norm, name))
    patch_edges = {(min(a, b), max(a, b)) for _, (_, e) in PATCH_EDGES.items() for a, b in e}
    for key, uses in edge_info.items():
        if len(uses) == 2:
            (n_a, g_a, b_a), (n_b, g_b, b_b) = uses
            assert n_a == n_b, f"cell counts differ on edge {key}: {b_a} {n_a} vs {b_b} {n_b}"
            assert math.isclose(g_a, g_b, rel_tol=1e-12), f"grading differs on edge {key}"
            assert key not in patch_edges, f"interior edge {key} listed as a patch"
        elif len(uses) == 1:
            assert key in patch_edges, f"outer edge {key} of {uses[0][2]} not in any patch"
        else:
            raise AssertionError(f"edge {key} used by {len(uses)} blocks")
    for patch, (_, edges) in PATCH_EDGES.items():
        for a, b in edges:
            key = (min(a, b), max(a, b))
            assert key in edge_info and len(edge_info[key]) == 1, f"{patch} edge {key} invalid"
            # fluid on the left of a->b  <=>  outward normal = (dy, -dx) points out of the block
            (xa, ya), (xb, yb) = POINTS_2D[a], POINTS_2D[b]
            normal = np.array([yb - ya, -(xb - xa)])
            mid = 0.5 * (np.array(POINTS_2D[a]) + np.array(POINTS_2D[b]))
            inside = _block_centroid_of_edge(m, key)
            assert np.dot(normal, mid - inside) > 0, f"{patch} edge {a}->{b} has inward normal"


def _block_centroid_of_edge(m: MeshLevel, key: tuple[int, int]) -> np.ndarray:
    """Centroid of the (unique) block that owns boundary edge ``key``."""
    for _, quad, _, _ in blocks(m):
        pairs = {(min(quad[i], quad[(i + 1) % 4]), max(quad[i], quad[(i + 1) % 4])) for i in range(4)}
        if key in pairs:
            return np.mean([POINTS_2D[k] for k in quad], axis=0)
    raise KeyError(key)


def first_last_cell(length: float, n: int, ratio: float) -> tuple[float, float]:
    """First and last cell sizes of a geometric grading (last/first = ratio)."""
    if math.isclose(ratio, 1.0):
        return length / n, length / n
    k = ratio ** (1.0 / (n - 1))
    first = length * (k - 1.0) / (k**n - 1.0)
    return first, first * ratio


# ---------------------------------------------------------------- writers
HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
  Generated by cfd/generate_cases.py  -  2-D cylinder wake, Re = 100
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    {location}object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""


def foam_header(cls: str, obj: str, location: str | None = None) -> str:
    """Standard OpenFOAM dictionary header."""
    loc = f'location    "{location}";\n    ' if location else ""
    return HEADER.format(cls=cls, obj=obj, location=loc)


def block_mesh_dict(m: MeshLevel) -> str:
    """Text of system/blockMeshDict for mesh level ``m``."""
    lines = [foam_header("dictionary", "blockMeshDict"), "convertToMeters 1;\n", "vertices\n("]
    for z in (Z0, Z1):
        for k in range(NP):
            x, y = POINTS_2D[k]
            lines.append(f"    ({x:.10g} {y:.10g} {z:g})    // {k + (0 if z == Z0 else NP)}")
    lines.append(");\n\nblocks\n(")
    for name, (v0, v1, v2, v3), (n1, n2), (g1, g2) in blocks(m):
        hexv = f"{v0} {v1} {v2} {v3} {v0 + NP} {v1 + NP} {v2 + NP} {v3 + NP}"
        lines.append(f"    hex ({hexv}) ({n1} {n2} 1) simpleGrading ({g1:.10g} {g2:.10g} 1)    // {name}")
    lines.append(");\n\nedges\n(")
    for a, b, (xm, ym) in ARCS:
        for off, z in ((0, Z0), (NP, Z1)):
            lines.append(f"    arc {a + off} {b + off} ({xm:.10g} {ym:.10g} {z:g})")
    lines.append(");\n\nboundary\n(")
    for patch, (ptype, edges) in PATCH_EDGES.items():
        lines.append(f"    {patch}\n    {{\n        type {ptype};\n        faces\n        (")
        for a, b in edges:
            lines.append(f"            ({a} {b} {b + NP} {a + NP})")
        lines.append("        );\n    }")
    lines.append("    frontAndBack\n    {\n        type empty;\n        faces\n        (")
    for _, (v0, v1, v2, v3), _, _ in blocks(m):
        lines.append(f"            ({v0} {v3} {v2} {v1})")                        # z-min, normal -z
        lines.append(f"            ({v0 + NP} {v1 + NP} {v2 + NP} {v3 + NP})")    # z-max, normal +z
    lines.append("        );\n    }\n);\n\nmergePatchPairs\n(\n);\n")
    return "\n".join(lines)


INLET_KICK = """    inlet
    {
        // U = (1, 0) as specified, except for a small transverse "kick"
        // (v = 0.1 for t < 5, ramped to 0 at t = 6) that breaks the
        // symmetry so that vortex shedding starts sooner.
        type            uniformFixedValue;
        uniformValue    table
        (
            (0     (1 0.1 0))
            (5     (1 0.1 0))
            (6     (1 0   0))
            (10000 (1 0   0))
        );
    }"""

INLET_PLAIN = """    inlet
    {
        type            fixedValue;
        value           uniform (1 0 0);
    }"""


def u_file(kick: bool = True) -> str:
    """Text of 0/U (with or without the start-up kick at the inlet)."""
    return U_FILE_TEMPLATE.replace("@INLET@", INLET_KICK if kick else INLET_PLAIN)


U_FILE_TEMPLATE = foam_header("volVectorField", "U", "0") + """
dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (1 0 0);

boundaryField
{
@INLET@
    outlet
    {
        type            zeroGradient;
    }
    top
    {
        type            symmetryPlane;
    }
    bottom
    {
        type            symmetryPlane;
    }
    cylinder
    {
        type            noSlip;
    }
    frontAndBack
    {
        type            empty;
    }
}
"""

P_FILE = foam_header("volScalarField", "p", "0") + """
dimensions      [0 2 -2 0 0 0 0];     // kinematic pressure p/rho

internalField   uniform 0;

boundaryField
{
    inlet
    {
        type            zeroGradient;
    }
    outlet
    {
        type            fixedValue;
        value           uniform 0;
    }
    top
    {
        type            symmetryPlane;
    }
    bottom
    {
        type            symmetryPlane;
    }
    cylinder
    {
        type            zeroGradient;
    }
    frontAndBack
    {
        type            empty;
    }
}
"""

def transport(nu: float = NU) -> str:
    """Text of constant/transportProperties for viscosity ``nu``."""
    return foam_header("dictionary", "transportProperties", "constant") + f"""
transportModel  Newtonian;

nu              [0 2 -1 0 0 0 0] {nu:g};    // Re = U D / nu = {U_INF * 2 * R / nu:g}
"""


RUN_TIME_CONTROL = """// Stop the run automatically once the steady-state criterion is met: every
// initial residual (U components and p) below the threshold [verify the
// runTimeControl syntax for your OpenFOAM version].
runTimeControl1
{{
    type            runTimeControl;
    libs            (utilityFunctionObjects);
    conditions
    {{
        steadyState
        {{
            type            equationInitialResidual;
            fields          (U p);
            value           {value:g};
            mode            maximum;
        }}
    }}
    satisfiedAction end;
}}
"""

TURBULENCE = foam_header("dictionary", "turbulenceProperties", "constant") + """
simulationType  laminar;
"""

FV_SCHEMES = foam_header("dictionary", "fvSchemes", "system") + """
ddtSchemes
{
    default         backward;                  // 2nd-order implicit
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;
    div(phi,U)      Gauss linearUpwind grad(U);   // 2nd-order upwind-biased
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}
"""

FV_SOLUTION = foam_header("dictionary", "fvSolution", "system") + """
solvers
{
    p
    {
        solver          GAMG;
        smoother        DIC;
        tolerance       1e-06;
        relTol          0.01;
    }

    pFinal
    {
        $p;
        tolerance       1e-07;
        relTol          0;
    }

    "(U|UFinal)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-08;
        relTol          0;
    }
}

PIMPLE
{
    momentumPredictor        yes;
    nOuterCorrectors         1;       // = PISO
    nCorrectors              2;
    nNonOrthogonalCorrectors 1;
}
"""

DECOMPOSE = foam_header("dictionary", "decomposeParDict", "system") + """
numberOfSubdomains  4;

method              scotch;
"""

FORCE_COEFFS = """// Force coefficients on the cylinder (included from controlDict)
forceCoeffs1
{
    type            forceCoeffs;
    libs            (forces);
    writeControl    timeStep;
    writeInterval   1;
    log             false;

    patches         (cylinder);
    rho             rhoInf;          // incompressible solver: use rhoInf
    rhoInf          1;
    CofR            (0 0 0);
    liftDir         (0 1 0);
    dragDir         (1 0 0);
    pitchAxis       (0 0 1);
    magUInf         1;
    lRef            1;               // D
    Aref            0.1;             // D x span (span = cell thickness 0.1)
}
"""


def sample_window_dict(start: float = SAMPLE_START, interval: float = SAMPLE_DT) -> str:
    """Function object sampling U and p on the 81 x 81 PIV window."""
    x = np.linspace(WINDOW["x0"], WINDOW["x1"], WINDOW["nx"])
    y = np.linspace(WINDOW["y0"], WINDOW["y1"], WINDOW["ny"])
    pts = "\n".join(f"                ({xi:.6g} {yj:.6g} 0)" for yj in y for xi in x)
    return f"""// Run-time sampling of the PIV window (fine mesh only; included from controlDict).
// Set type "cloud" is the ESI (openfoam.com) name; Foundation versions call it
// "points" (with "ordered no;") [verify for your version].
sampleWindow
{{
    type            sets;
    libs            (sampling);
    timeStart       {start:g};
    writeControl    adjustableRunTime;   // also adjusts dt to hit the sample times
    writeInterval   {interval:g};
    interpolationScheme cellPoint;
    setFormat       raw;
    fields          (p U);

    sets
    {{
        window
        {{
            type    cloud;
            axis    xyz;
            points
            (
{pts}
            );
        }}
    }}
}}
"""


def control_dict(m: MeshLevel) -> str:
    """Text of system/controlDict."""
    extra = '    #include "sampleWindow"\n' if m.sample else ""
    if m.residual_stop is not None:
        extra += '    #include "runTimeControl"\n'
    return foam_header("dictionary", "controlDict", "system") + f"""
application     pimpleFoam;

startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {m.end_time:g};

deltaT          0.001;               // initial value; adjusted for maxCo
adjustTimeStep  yes;
maxCo           0.5;
maxDeltaT       0.01;

writeControl    adjustableRunTime;
writeInterval   25;
purgeWrite      2;                   // keep only the last two field dumps

writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   8;

runTimeModifiable true;

functions
{{
    #include "forceCoeffs"
{extra}}}
"""


ALLRUN = """#!/bin/sh
# Mesh and run one case.  Usage: ./Allrun [nProcs]   (default 1 = serial)
cd "${0%/*}" || exit 1
NPROCS=${1:-1}
blockMesh > log.blockMesh 2>&1 || { echo "blockMesh failed (see log.blockMesh)"; exit 1; }
checkMesh > log.checkMesh 2>&1
if [ "$NPROCS" -gt 1 ]; then
    foamDictionary -entry numberOfSubdomains -set "$NPROCS" system/decomposeParDict
    decomposePar -force > log.decomposePar 2>&1
    mpirun -np "$NPROCS" pimpleFoam -parallel > log.pimpleFoam 2>&1
    reconstructPar -latestTime > log.reconstructPar 2>&1
else
    pimpleFoam > log.pimpleFoam 2>&1
fi
echo "$(basename "$PWD"): finished (see log.pimpleFoam)"
"""

ALLCLEAN = """#!/bin/sh
# Remove results of this case (keeps the set-up files).
cd "${0%/*}" || exit 1
rm -rf processor* postProcessing constant/polyMesh log.* [1-9]* 0.[0-9]*
"""


def write_case(m: MeshLevel) -> Path:
    """Write the complete case directory for one mesh level."""
    check_topology(m)
    root = CASES / m.name
    files = {
        "0/U": u_file(m.kick), "0/p": P_FILE,
        "constant/transportProperties": transport(m.nu),
        "constant/turbulenceProperties": TURBULENCE,
        "system/blockMeshDict": block_mesh_dict(m),
        "system/controlDict": control_dict(m),
        "system/fvSchemes": FV_SCHEMES, "system/fvSolution": FV_SOLUTION,
        "system/decomposeParDict": DECOMPOSE, "system/forceCoeffs": FORCE_COEFFS,
        "Allrun": ALLRUN, "Allclean": ALLCLEAN,
    }
    if m.sample:
        files["system/sampleWindow"] = sample_window_dict(m.sample_start, m.sample_dt)
    if m.residual_stop is not None:
        files["system/runTimeControl"] = RUN_TIME_CONTROL.format(value=m.residual_stop)
    for rel, text in files.items():
        path = root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    for script in ("Allrun", "Allclean"):
        p = root / script
        p.chmod(p.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    first_r, last_r = first_last_cell(A - R, m.n_r, m.e_r)
    first_d, _ = first_last_cell(X_OUT - A, m.n_down, m.e_down)
    info = {
        **asdict(m), "cells": m.n_cells(),
        "cylinder_tangential_cell": math.pi * R / 2 / m.n_t,
        "cylinder_first_radial_cell": first_r,
        "ogrid_outer_radial_cell_on_axis": last_r,
        "wake_cell_dy": 2 * A / m.n_t, "wake_first_dx": first_d,
        "U_inf": U_INF, "D": 2 * R, "Re": U_INF * 2 * R / m.nu,
        "window": WINDOW if m.sample else None,
    }
    (root / "mesh_info.json").write_text(json.dumps(info, indent=2))
    return root


def main() -> None:
    """Write the three Re = 100 cases and the optional steady Re = 20 case."""
    print(f"{'case':12s} {'Re':>5s} {'cells':>8s} {'cyl ds':>8s} {'1st dr':>8s} {'wake dy':>8s}")
    for m in LEVELS + OPTIONAL_LEVELS:
        root = write_case(m)
        info = json.loads((root / "mesh_info.json").read_text())
        print(f"{m.name:12s} {info['Re']:5.0f} {info['cells']:8d} {info['cylinder_tangential_cell']:8.4f} "
              f"{info['cylinder_first_radial_cell']:8.4f} {info['wake_cell_dy']:8.4f}   -> {root}")


if __name__ == "__main__":
    main()
