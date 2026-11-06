from pathlib import Path
import json,ast,time
from sage.all import GF,ZZ,matrix,vector,euler_phi
H=Path('work/deg10');P=H/'continuation_20261029b';tt=time.monotonic()
helper=(H/'computations/late_replay_filters_v2.py').read_text().split('# The two low-genus filters')[0];exec(helper)
rankrows=[]
for N,h,g in [(38,1,28),(40,1,25),(50,1,48),(76,39,64),(80,41,57),(100,51,109)]:
 a=read(P/f'winding_N{N}_H{h}.m');assert a[:4]==[N,h,g,g];total=0
 for idx,o,d,p,v,ops in a[4]:
  F=GF(p);assert F.is_prime_field();V=F**d;U=V.subspace([vector(F,v)]);A=[matrix(F,d,d,z) for z in ops]
  while True:
   W=V.subspace(U.basis()+[u*T for u in U.basis() for T in A])
   if W.dimension()==U.dimension():break
   U=W
  assert U.dimension()==d;total+=d*euler_phi(o)
 assert total==g;rankrows.append({'N':N,'H_generator':h,'genus':g,'all_winding_components_full':True})
group,source,small,cover,images,tests,remaining=read(P/'C250/filter25.m');iv,rel,hg,ds,cc,native=group;G25=Quot(rel);H25,A25=allowed(rel,hg,ds,cc);assert G25.index==2503105 and H25.index==11
F50={d:field(d,p) for d,p,c in cover};own={(r[0],r[1]):r for r in source};assert len(own)==137 and not small;cache={}
for kk,e,w,cl in images:
 m,d=kk[:2];E,T,Q=src(d,kk[2:],F50[d]);V=2*Q;assert 50*Q==E(0) and 25*Q!=E(0) and 10*Q!=E(0)
 if m==1:EP=E;U=V
 else:EP,(U,)=iso2(E,[T,25*Q,T+25*Q][m-2],[V])
 assert 25*U==EP(0) and 5*U!=EP(0) and deg(tate(EP,U),d)==e and w==d//e;cache[tuple(kk)]=cl
expected=[[1,d,i] for d,i in own]
for k,desc,cd,checks,passed in tests:
 assert desc==expected[k-1];r=own[desc[1],desc[2]];assert cd==10-r[0]
 for m,cl,yes in checks:
  actual=cache[tuple([m,r[0]]+r[2])];assert G25.contains([a-b for a,b in zip(actual,cl)]) and bool(yes)==(H25.residue([-x for x in cl]) in A25[cd])
 assert bool(passed)==all(x[2] for x in checks)
assert remaining==[r[1] for r in tests if r[4]] and len(remaining)==3
for k in range(3):
 im,tests,left=read(P/f'C250/hecke_shard{k}.m');assert not left and len(tests)==1
 for kk,e,w,xy in im:
  m,d=kk[:2];E,T,Q=src(d,kk[2:],F50[d]);assert m==1;v=tate(E,Q);assert [[int(x) for x in z._vector_()] for z in v]==xy and deg(v,d)==e and w==d//e
 for idx,desc,cd,checks,passed in tests:
  assert idx==k+1 and desc==remaining[k] and not passed and cd==2 and checks==[[1,7,64,False]]
base=json.loads((P/'four_exclusions_replay.json').read_text());assert base['status']=='PASS'
report={'status':'PASS','new_exclusions':[[1,76],[1,80],[1,100],[2,38],[2,40],[2,50]],'C238_supports':[233,35,0],'C240_supports':[781,288,118,0],'C250_supports':[137,3,0],'rank_certificates':rankrows,'C250_target25_rational_group_order':int(G25.index//H25.index),'C250_target25_image_index':int(H25.index),'no_unresolved_case_promoted':True,'elapsed_seconds':time.monotonic()-tt}
(P/'six_exclusions_replay.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
