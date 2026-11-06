from pathlib import Path
from math import gcd,isqrt
from collections import Counter
import ast,re,json,time
from sage.all import QuadraticField
C=Path(__file__).resolve().parents[1]/'continuation_20261031a/C90';start=time.monotonic();classes=[]
for D in[-119,-319]:
 forms=[]
 for a in range(1,isqrt((-D)//3)+1):
  for b in range(-a,a+1):
   if (b*b-D)%(4*a):continue
   c=(b*b-D)//(4*a)
   if a>c or gcd(gcd(a,b),c)!=1 or ((abs(b)==a or a==c)and b<0):continue
   forms.append([a,b,c])
 F=QuadraticField(D);assert F.discriminant()==D and F.class_number()==len(forms)==10;classes.append({'discriminant':D,'forms':forms,'class_number':10,'independent_PARI_class_number':10})
assert QuadraticField(35**2-4096).discriminant()==-319 and QuadraticField((-55)**2-4096).discriminant()==-119
for h in[1,2,5,10]:
 if h%2:continue
 assert 2**(h//2)%45 not in[1,44]
# Independent full rank at each leaf, without incremental basis reuse.
s=(C/'inputs/old45/N45_geometry/cyclic_sieves/N45_completion_jets.m').read_text();ds=ast.literal_eval(re.findall(r'\bds:=(.*?);',s,re.S)[0]);blocks=[]
for i,d in enumerate(ds,1):
 m=8//d;nr,nc,ss=re.findall(r'localjets\[<'+str(i)+','+str(m)+r'>\]:=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S)[0];nr=int(nr);nc=int(nc);es=ast.literal_eval(ss);assert nr==39;blocks.append([sum(es[nc*r+j]<<r for r in range(nr))for j in range(nc)])
def rank(rs):
 p={}
 for x in rs:
  while x:
   i=x.bit_length()-1
   if i in p:x^=p[i]
   else:p[i]=x;break
 return len(p)
count=0
def enum(i,left,rows,has4):
 global count
 if i==len(ds):
  if left==0 and has4:assert len(rows)==10 and rank(rows)==10;count+=1
  return
 for m in range(min(left//ds[i],len(blocks[i])//ds[i])+1):enum(i+1,left-m*ds[i],rows+blocks[i][:m*ds[i]],has4 or(bool(m)and ds[i]==4))
enum(0,10,[],False);assert count==94420
# Native old degree9 success and degree8 accepted nonprincipality evidence.
old=C/'inputs/old45';log=(old/'N45_good9_and_pure9_retry.log').read_text();assert 'NONIC45_GOOD9_AND_PURE9_DONE' in log and len(re.findall(r'NONIC45_GOOD9_CASE index=\d+ q=7 principal=false',log))==5
res={'status':'PASS','CM_class_numbers':classes,'source_decic_orders':[990,1080],'Frobenius_traces':[35,-55],'possible_automorphism_orders':[1,2,5,10],'marked_cycle_even_cases':[[2,2],[10,32]],'marked_scalars_never_pm1_mod45':True,'independent_full_matrix_quartic_patterns':94420,'all_ranks':10,'old_degree9_A7_tests_replayed_from_accepted_native_log':5,'elapsed_seconds':time.monotonic()-start}
(C/'independent_cycle_and_boundary_controls.json').write_text(json.dumps(res,indent=2)+'\n');print('D31A_C90_CYCLE_BOUNDARY_CONTROLS_COMPLETE',10,10,94420,time.monotonic()-start,flush=True)
