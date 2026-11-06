#!/usr/bin/env python3
"""Read-only point-convention correction, preserving every existing matrix."""
from pathlib import Path
import json,hashlib,sys
H=Path(__file__).resolve().parent;R=H.parents[2];m=json.loads((H/('manifest.candidate.json'if '--candidate'in sys.argv else'manifest.json')).read_text());assert m['status']=='CONVENTION-SAFE POINT-HECKE PROOF CORRECTION'
for a in m['files']:assert hashlib.sha256((R/a['path']).read_bytes()).hexdigest()==a['sha256'],a['path']
for a in m['source_manifests']:
 p=R/a['path'];assert hashlib.sha256(p.read_bytes()).hexdigest()==a['sha256'];old=json.loads(p.read_text())
 for f in old['files']:assert hashlib.sha256((R/f['path']).read_bytes()).hexdigest()==f['sha256'],f['path']
assert not any('P83/fricke_convention_supplement'in a['path']for a in m['files']+m['source_manifests'])
C=H.parent
for n,folder,cases,rank,needed in[(83,'P83/canonical_units',410,247,205),(85,'composite85',320,225,160),(91,'composite91',288,265,144),(65,'composite65',336,121,168)]:
 r=json.loads((C/folder/'canonical_R1_sieve_verified.json').read_text());assert r['status']=='PASS'and r['cases']==cases and not r['unresolved'];assert all(v['rank']==rank for v in r['results']);assert sum(v['convention']==2 for v in r['results'])==needed
 p=H/(f'case{n}_corrected.tex'if n!=65 else'case65_ordinary_corrected.tex');s=p.read_text();assert('A_3=T_3-1-3\\tau'in s)and'Fricke-transported point-level operator'not in s and'F T_3^{\\rm point}'not in s
 if n==83:assert'B=\\sigma T-T-\\sigma P+P-3\\sigma^{11}P+3\\sigma^{10}P'in s
 elif n==85:assert'B=\\sigma T-T-\\sigma P+P-3\\sigma\\tau P+3\\tau P'in s
 elif n==91:assert'B=\\sigma T-T-\\sigma P-2P+3\\tau P'in s
 else:assert'B_2=\\sigma T-T-\\sigma P-2P+3\\tau P.'in s
print('NONIC_POINT_HECKE_CONVENTIONS_READONLY_PASS all_original_B2_cases677 unchanged_fullrank_matrices1354=true')
