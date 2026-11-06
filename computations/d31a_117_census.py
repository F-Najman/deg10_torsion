from pathlib import Path
import ast,json,hashlib,subprocess,sys,time
from sage.all import GF,PolynomialRing,EllipticCurve,divisors,euler_phi
H=Path(__file__).resolve().parents[1];C=H/'continuation_20261031a/C117';t0=time.time()
raw=C/'inputs/all_curves.m';meta=json.loads((C/'inputs/independent_point_counts.json').read_text());assert meta['status']=='PASS' and hashlib.sha256(raw.read_bytes()).hexdigest()==meta['input_sha256']
data=ast.literal_eval(raw.read_text().replace('[*','[').replace('*]',']'));R=PolynomialRing(GF(2),'x');rows=[]
for d,poly,curves in data:
 F=GF(2**d,'a',modulus=R(poly));fv=lambda v:R(v)(F.gen());eligible=[]
 for jj,aa,card in curves:
  if card%117:continue
  E=EllipticCurve(F,[fv(v)for v in aa]);assert E.cardinality()==card
  inv=list(map(int,E.abelian_group().invariants()));assert inv[-1]%39==0
  eligible.append({'j':jj,'cardinality':int(card),'group':inv,'has117':inv[-1]%117==0,'j_degree':int(next(e for e in divisors(d)if fv(jj)**(2**e)==fv(jj)))})
 rows.append({'degree':int(d),'eligible_curve_orders':sorted(set(e['cardinality']for e in eligible)),'eligible_curves':eligible})
 print('D31A_C117_FINITE_DEGREE',d,len(eligible),rows[-1]['eligible_curve_orders'],flush=True)
assert all(not r['eligible_curves']for r in rows[:7]);assert rows[7]['eligible_curve_orders']==[234] and rows[8]['eligible_curve_orders']==[468] and not rows[9]['eligible_curves']
assert len(rows[8]['eligible_curves'])==3 and all(r['j_degree']==3 for r in rows[8]['eligible_curves'])
script=C/'inputs/work/degree9/cyclic/continuation_20261027/C117_X139/certify.py'
r=subprocess.run([sys.executable,str(script)],text=True,capture_output=True,timeout=400);(C/'inherited_replay.log').write_text(r.stdout+r.stderr);assert r.returncode==0 and 'C117_D9_COMPLETE_READ_ONLY_CERTIFICATE_PASS' in r.stdout,(r.returncode,r.stdout,r.stderr)
result={'status':'PASS','scope':'Degree-independent X1(39) dyadic injection, intact cusps, lower A5 tests, known degree8 and9 lift bindings; new degree10 necessary117 finite census','all_original_census_curves_present':True,'finite_fields':rows,'source39_degree10_reduction_cannot_support117':True,'old_packet_actual_exit_code':r.returncode,'elapsed_seconds':time.time()-t0}
(C/'finite_and_inherited_verified.json').write_text(json.dumps(result,indent=2)+'\n');print('D31A_C117_CENSUS_AND_INHERITED_COMPLETE',time.time()-t0,flush=True)
