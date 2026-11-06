from pathlib import Path
from math import gcd,isqrt
import json,hashlib,re,ast
C=Path(__file__).resolve().parent;B=C.parents[1];m=json.loads((C/'final_certificate_summary.json').read_text());assert m['status']=='PASS' and m['degree']==10 and m['excluded_groups']==[[1,90]] and m['parent_mathematical_review']=='PASS'
for n,h in m['files_sha256'].items():assert hashlib.sha256((B/n).read_bytes()).hexdigest()==h,n
r=json.loads((C/'independent_cycle_and_boundary_controls.json').read_text());assert r['status']=='PASS' and r['independent_full_matrix_quartic_patterns']==94420 and r['all_ranks']==10
for D in[-119,-319]:
 forms=[]
 for a in range(1,isqrt((-D)//3)+1):
  for b in range(-a,a+1):
   if (b*b-D)%(4*a):continue
   c=(b*b-D)//(4*a)
   if a>c or gcd(gcd(a,b),c)!=1 or((abs(b)==a or a==c)and b<0):continue
   forms.append([a,b,c])
 assert len(forms)==10
 assert next(v for v in r['CM_class_numbers']if v['discriminant']==D)['forms']==forms
for h in[2,10]:assert 2**(h//2)%45 not in[1,44]
g=json.loads((C/'independent_genuine_and_known_lift.json').read_text());assert g['status']=='PASS' and g['independent_source_closed_places']==108 and g['independent_source_diamond_orbits']==10
raw=ast.literal_eval((C/'inputs/all_curves.m').read_text().replace('[*','[').replace('*]',']'));d10=next(row for row in raw if row[0]==10);assert sorted({row[2]for row in d10[2]if row[2]%45==0})==[990,1080]
old=C/'inputs/old45';log=(old/'N45_quotient_injectivity_control.log').read_text();assert 'EQUAL_TWO_PRIMARY_ORDERS 8192 8192'in log and 'INJECTIVE_REDUCTION_3_15_DONE'in log
st=json.loads((old/'N45_quotient_injectivity_control_status.json').read_text());assert st['native_exit_code']==0 and st['log_sha256']==hashlib.sha256((old/'N45_quotient_injectivity_control.log').read_bytes()).hexdigest()
log=(old/'N45_good9_and_pure9_retry.log').read_text();assert len(re.findall(r'NONIC45_GOOD9_CASE index=\d+ q=7 principal=false',log))==5 and 'NONIC45_GOOD9_AND_PURE9_DONE'in log
for row in m['jobs']:
 st=json.loads((B/row['status_path']).read_text());assert st['state']=='EXITED' and st['actual_exit_code']==row['actual_exit_code']
 if row['accepted_completion_marker']:assert row['accepted_completion_marker']in(B/st['log']).read_text()
assert m['new_decic_A7_A11_RR_not_premises']
print('D31A_C90_FINAL_CERTIFICATE_COMPLETE files='+str(len(m['files_sha256']))+' quarticpatterns94420=true CMclassnumbers10and10=true markedcycle=true')
