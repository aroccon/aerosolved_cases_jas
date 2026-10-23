#!/usr/bin/env python3
"""Particles for the frozen-HIT case: one positions file, particles initially at rest.

EL: <case>/constant/kinematicCloudPosition with N random positions in the periodic box,
    used by a single manualInjection injector with U0 = (0 0 0) and a uniform size
    distribution (set in constant/kinematicCloudProperties).
EE: <case>/0/V = -U, i.e. particle velocity relative to the gas for particles at rest,
    so that both frameworks start from the same state.

Old per-particle files (kinematicCloudPosition<i>, constInjectorData) are removed.

Usage:
    python3 make_particles.py hit_frozen_150 --N 50000 --L 0.01
"""
import argparse
import glob
import os
import re

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('case', help='case folder (reads <case>/0/U)')
ap.add_argument('--N', type=int, default=50000, help='number of particles')
ap.add_argument('--L', type=float, default=0.01, help='box size [m]')
ap.add_argument('--seed', type=int, default=1)
a = ap.parse_args()

cdir = os.path.join(a.case, 'constant')

# --- positions --------------------------------------------------------------------
rng = np.random.default_rng(a.seed)
x = a.L * rng.random((a.N, 3))

for f in glob.glob(os.path.join(cdir, 'kinematicCloudPosition*')) + \
         glob.glob(os.path.join(cdir, 'constInjectorData')):
    os.remove(f)
with open(os.path.join(cdir, 'kinematicCloudPosition'), 'w') as fh:
    fh.write('FoamFile\n{\n    version     2.0;\n    format      ascii;\n'
             '    class       vectorField;\n    location    "constant";\n'
             '    object      kinematicCloudPosition;\n}\n\n'
             f'{a.N}\n(\n')
    np.savetxt(fh, x, fmt='(%.8e %.8e %.8e)')
    fh.write(');\n')
print(f'Written {a.N} positions in [0, {a.L:g}]^3 to {cdir}/kinematicCloudPosition')

# --- 0/V = -U ---------------------------------------------------------------------
txt = open(os.path.join(a.case, '0', 'U')).read()
m = re.search(r'internalField\s+nonuniform\s+List<vector>\s*(\d+)\s*\(', txt)
nc = int(m.group(1))
end = txt.index('\n)', m.end())
U = np.array(txt[m.end():end].replace('(', ' ').replace(')', ' ').split(), dtype=float)
if U.size != 3 * nc:
    raise SystemExit(f'0/U: expected {3 * nc} values, read {U.size}')
U = U.reshape(nc, 3)
with open(os.path.join(a.case, '0', 'V'), 'w') as fh:
    fh.write('FoamFile\n{\n    version     2.0;\n    format      ascii;\n'
             '    class       volVectorField;\n    object      V;\n}\n\n'
             'dimensions      [0 1 -1 0 0 0 0];\n\n'
             f'internalField   nonuniform List<vector>\n{nc}\n(\n')
    np.savetxt(fh, -U, fmt='(%.8e %.8e %.8e)')
    fh.write(');\n\nboundaryField\n{\n    ".*"\n    {\n        type cyclic;\n    }\n}\n')
print(f'Written 0/V = -U ({nc} cells, rms = {np.sqrt(np.mean(U ** 2)):.4g} m/s)')
