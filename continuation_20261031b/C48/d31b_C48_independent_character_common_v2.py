from sage.all import *
from pathlib import Path
from collections import Counter
import json,re,time,sys
P=Path('work/deg10/continuation_20261031b/C48');T=time.monotonic();Fp=GF(5);RP=PolynomialRing(Fp,'t');sh,nsh,limit=map(int,sys.argv[1:4])
def rd(n):return json.loads((P/n).read_text().replace('[*','[').replace('*]',']'))
src=rd('new_high_sources.m');pri=rd('high_character_merged.json');primary={r[0]:r for r in pri[1]};fs={d:GF(5**d,'a',modulus=RP(pol))for d,pol in src[0]};text=(P/'inputs/four_characters.m').read_text()
def var(n):return json.loads(re.search(n+r'\s*:=\s*(.*?);',text,re.S)[1].replace('[*','[').replace('*]',']'))
fun=var('FilterFunctions');base=var('FilterBaseValues');comb=var('FilterCombinations');need=[i for i in range(8)if any(row[i]for row in comb)];fun=[[(RP(n),RP(d))for n,d in f]for f in fun];done=0;kept=0;tests=0;digest=[];fallback=0
for ii,(id,d,ky)in enumerate(src[1]):
 if ii%nsh!=sh:continue
 if limit and done>=limit:break
 K=fs[d];a=K.gen();b=RP(ky[:d])(a);c=RP(ky[d:])(a);E=EllipticCurve(K,[1-c,-b,-b,0,0]);Q=E(0,0);S=24*Q;assert S!=E(0)and 2*S==E(0)and 16*Q!=E(0);Z=Q+S;r=S[0];xx=Q[0]+Z[0]-S[0];yy=Q[1]+Z[1]-S[1]+((1-c)*xx-b)/2;A2=-b+(1-c)**2/4;r2=-A2-2*r;delta=xx-r2;bb=delta**3/yy**2;cc=(3*r2+A2)*delta**2/yy**2
 if done<3 or done%1000==0:
  R=PolynomialRing(K,'X');ph=E.isogeny(R.gen()-r);ET=ph.codomain();QT=ph(Q);assert QT.order()==24;assert list(QT.xy())==[Q[0]+Z[0]-S[0],Q[1]+Z[1]-S[1]];ai=ET.a_invariants();assert r2**3+(ai[1]+ai[0]**2/4)*r2*r2+(ai[3]+ai[0]*ai[2]/2)*r2+ai[4]+ai[2]**2/4==0;assert QT[0]==xx and QT[1]+(ai[0]*xx+ai[2])/2==yy;tests+=1
 t=2*bb/(1-cc)-1;q=((1-cc)*t*t-2)/(2*t);ss=t*(q*q-1)/(t*(q*q+1)+2*q);dd=(ss+2)/q;u,v=(ss+dd)/2,(ss-dd)/2;sg={}
 if not primary[id][4]:
  assert primary[id][5]and any(den(v)==0 for i in need for num,den in fun[i]);fallback+=1;done+=1;kept+=1;continue
 for i in need:
  # Horner evaluation in the opposite variable from the native producer.
  val=K(0)
  for num,den in reversed(fun[i]):val=val*u+num(v)/den(v)
  assert val!=0;z=(val/K(base[i]))**((K.order()-1)//2);assert z in[1,-1];sg[i]=int(z==-1)
 chi=[sum(sg[i]*row[i]for i in need)%2 for row in comb];assert chi==primary[id][4];kept+=primary[id][5];done+=1
 if done%1000==0:print('D31B_C48_INDEPENDENT_CHARACTER_PROGRESS',sh,done,kept,tests,time.monotonic()-T,flush=True)
out={'status':'PASS','shard':sh,'shards':nsh,'limit':limit,'source_orbits':done,'admitted_orbits':int(kept),'actual_normalized_Velu_samples':tests,'all_character_vectors_match_primary':True,'conservatively_retained_poles':fallback,'method':'Actual point24P, manual two-point Velu coordinates, independently derived dual kernel and Derickx normalization, opposite-variable Horner evaluation of even-divisor functions and exact finite-field signs','elapsed_seconds':time.monotonic()-T};name=f'independent_character_s{sh}'if not limit else'independent_character_profile';(P/(name+'.json')).write_text(json.dumps(out,indent=2)+'\n');print('D31B_C48_INDEPENDENT_CHARACTER_COMPLETE',sh,done,kept,tests,time.monotonic()-T,flush=True)
