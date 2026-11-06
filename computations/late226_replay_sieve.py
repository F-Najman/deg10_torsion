from pathlib import Path
import json,ast,time
from sage.all import ZZ,matrix
H=Path('work/deg10');P=H/'continuation_20261029b';tt=time.monotonic()
exec((H/'computations/late_replay_filters_v2.py').read_text().split('# The two low-genus filters')[0])
group,source,small,cover,images,tests,remaining=read(P/'C226/filter26.m');iv,rel,hg,ds,cc,native=group;G26=Quot(rel);H26,A26=allowed(rel,hg,ds,cc);F={d:field(d,p) for d,p,c in cover};own={(r[0],r[1]):r for r in source};sm={(r[0],r[1]):r for r in small};cache={}
for kk,e,w,cl in images:
 m,d=kk[:2];E,T,Q=src(d,kk[2:],F[d]);assert 26*Q==E(0) and 13*Q!=E(0) and 2*Q!=E(0) and T!=13*Q
 if m<=2:EP=E;U=Q if m==1 else T+Q
 else:EP,(U,)=iso2(E,T if m==3 else T+13*Q,[Q])
 assert 26*U==EP(0) and 13*U!=EP(0) and 2*U!=EP(0) and deg(tate(EP,U),d)==e and w==d//e;cache[tuple(kk)]=cl
supports=[[1,d,i] for d,i in own]+[[2,5,i,5,j] for d,i in own if d==5 for v,j in sm if v==5]
for k,desc,cd,checks,passed in tests:
 assert desc==supports[k-1];rs=[own[desc[1],desc[2]]]+([sm[desc[3],desc[4]]] if desc[0]==2 else []);assert cd==10-sum(r[0] for r in rs)
 for m,cl,yes in checks:
  actual=sums([cache[tuple([m,r[0]]+r[2])] for r in rs]);assert G26.contains([a-b for a,b in zip(actual,cl)]) and bool(yes)==(H26.residue([-x for x in cl]) in A26[cd])
 assert bool(passed)==all(r[2] for r in checks)
assert remaining==[r[1] for r in tests if r[4]] and len(remaining)==32
D=read(P/'C226/proven_picard.m');GG=Quot(D[5]);fc=D[4];fbdegs=D[2];seen=set();left=[];nimages=nw=0
for k in range(8):
 group,ims,tests,lefts=read(P/f'C226/picard{k}.m');iv,rel,hg,cc=group;HQ=Quot(rel,hg);assert GG.index==124501550040000 and HQ.index==350;C={tuple(key):cl for key,cl in ims};nimages+=len(ims)
 for key,cl in ims:
  d=key[0];E,T,Q=src(d,key[1:],F[d]);assert bc(E,T,Q)==[sum(F[d](v)*F[d].gen()**i for i,v in enumerate(key[1:1+d])),sum(F[d](v)*F[d].gen()**i for i,v in enumerate(key[1+d:]))]
 for i,desc,cd,cl,yes in tests:
  assert i not in seen and (i-1)%8==k and desc==remaining[i-1];seen.add(i);rs=[own[desc[1],desc[2]]]+([sm[desc[3],desc[4]]] if desc[0]==2 else []);assert cd==10-sum(r[0] for r in rs);actual=sums([C[tuple([r[0]]+r[2])] for r in rs]);assert GG.contains([a-b for a,b in zip(actual,cl)])
  allowed={HQ.residue([0]*len(cl))}
  if cd>=3:allowed.update(HQ.residue([-x for x in c]) for c in cc)
  assert bool(yes)==(HQ.residue(cl) in allowed)
  if yes:left.append(desc)
 ww=ast.literal_eval(ast.parse((P/f'C226/picard{k}_RR.log.witnesses.py').read_text()).body[0].value)
 for trial,d,result,ev,rv,va,b in ww:
  assert sum(e*fbdegs[j-1] for j,e in ev)==sum(e*fbdegs[j-1] for j,e in rv)==25 and sum(e*fbdegs[j-1] for j,e in va)==1
  actual=[0]*len(result)
  for sign,sp in [(1,ev),(-1,rv),(d,va)]:
   for j,e in sp:actual=[a+sign*e*x for a,x in zip(actual,fc[j-1])]
  assert GG.contains([a-b for a,b in zip(actual,result)]);nw+=1
assert seen==set(range(1,33)) and sorted(left)==sorted([[1,9,52],[2,5,1,5,19],[1,10,54]])
report={'status':'PASS','initial_patterns':len(supports),'after_four_maps':32,'after_full_picard':left,'Picard_group_order':int(GG.index),'rational_image_index':int(HQ.index),'images':nimages,'RR_witness_class_sums':nw,'elapsed_seconds':time.monotonic()-tt}
(P/'C226/independent_sieve_replay.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
