#!/usr/bin/env python3
"""Exact prose correction: geometric mu annihilator, unchanged binary FI."""
from pathlib import Path
import json,hashlib,sys
H=Path(__file__).resolve().parent;R=H.parents[2];C=H.parent;m=json.loads((H/('manifest.candidate.json'if '--candidate'in sys.argv else'manifest.json')).read_text());assert m['point_correction_manifest_sha256']==hashlib.sha256((C/'point_convention_correction_20261027/manifest.json').read_bytes()).hexdigest()
for a in m['files']:assert hashlib.sha256((R/a['path']).read_bytes()).hexdigest()==a['sha256']
s=(C/'point_convention_correction_20261027/case85_corrected.tex').read_text();needle='Compose this map with the standard torsion annihilator $T_3-1-3\\langle3\\rangle$ in the canonical cotangent convention.';new='On the geometric $\\mu_{85}$-level model compose this map with the torsion annihilator $T_3-\\langle3\\rangle-3$. Since $3\\equiv1\\pmod2$, its binary cotangent matrix is exactly the retained $T_3+I+D_3$ matrix.';assert s.replace(needle,new)==(H/'case85_final.tex').read_text()
r=json.loads((C/'composite85/FI_positive_independent_verified.json').read_text());assert r['rank']==122 and r['minimum_cusp_prefix_cost']==32 and r['all_degree9_prefix_configurations_independent'];print('C85_MU_FI_CONVENTION_READONLY_PASS binary133x225_rank122_prefix32_unchanged=true')
