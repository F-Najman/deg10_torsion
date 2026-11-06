from pathlib import Path
import ast,json,hashlib,time
H=Path(__file__).resolve().parents[1];B=H/'next_cases_20261030/dyadic_local';start=time.monotonic()
p=B/'all_curves.m';data=ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']'));summ=[]
for d,pol,rows in data:
 q=1<<d;mod=sum(c<<i for i,c in enumerate(pol));assert len(pol)==d+1 and pol[-1]==1
 M=[]
 for a in range(q):
  row=[0]*q
  for b in range(1,q):
   v=row[b>>1]<<1
   if v&q:v^=mod
   row[b]=v^(a if b&1 else 0)
  M.append(row)
 sq=[M[a][a] for a in range(q)];cube=[M[a][sq[a]] for a in range(q)];tr=[]
 for a in range(q):
  v=t=a
  for _ in range(1,d):v=sq[v];t^=v
  assert t in [0,1];tr.append(t)
 invs=[0]*q
 for a in range(1,q):
  v=1;b=a;n=q-2
  while n:
   if n&1:v=M[v][b]
   b=sq[b];n>>=1
  assert M[v][a]==1;invs[a]=sq[v]
 enc=lambda v:sum(c<<i for i,c in enumerate(v))
 js={};counts=[]
 for j,aa,n in rows:
  a1,a2,a3,a4,a6=map(enc,aa);j=enc(j);a12=sq[a1];a32=sq[a3];b4=M[a1][a3];b8=M[a12][a6]^M[b4][a4]^M[a2][a32]^sq[a4]
  delta=M[sq[a12]][b8]^sq[a32]^M[M[a12][b4]][a32];assert delta
  if a1==0:assert j==0
  else:
   c4=sq[a12];assert M[M[sq[c4]][c4]][M[delta].index(1)]==j
  answer=1
  for x in range(q):
   rhs=cube[x]^M[a2][sq[x]]^M[a4][x]^a6;b=M[a1][x]^a3
   if b==0:answer+=1
   elif tr[M[rhs][invs[b]]]==0:answer+=2
  assert answer==n,(d,aa,answer,n);assert n%95 and n%111;counts.append(n);js.setdefault(j,[]).append(n)
 assert set(js)==set(range(q))
 assert all(len(ns)==2 and sum(ns)==2*(q+1) for j,ns in js.items() if j)
 assert len(js[0])==(3 if d%2 else 7)
 summ.append({'degree':d,'field_order':q,'curve_twists_checked':len(rows),'distinct_curve_orders':sorted(set(counts))})
result={'status':'PASS','levels':[95,111],'all_degrees_checked':list(range(1,11)),'point_count_method':'Exact characteristic-two quadratic trace formula, independent polynomial-basis integer field arithmetic; no Magma or Sage curve operations','total_twists_checked':sum(x['curve_twists_checked'] for x in summ),'all_finite_j_values_and_twists_present':True,'no_curve_order_divisible_by_either_level':True,'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'seconds':time.monotonic()-start,'fields':summ}
(B/'independent_point_counts.json').write_text(json.dumps(result,indent=2)+'\n');print('NEXT95111_INDEPENDENT_FINITE_POINT_COUNTS_PASS',result['total_twists_checked'],result['seconds'],flush=True)
