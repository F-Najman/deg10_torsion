#!/usr/bin/env python3
"""Exact modular indices, quotient cusps, local 3-isogeny shape and packet degrees."""
from pathlib import Path
from math import gcd
import ast,json,hashlib
H=Path(__file__).resolve().parent

def cusp(N,a,c):
 c%=N;g=gcd(c,N);return min((a%g,c),((-a)%g,(-c)%N))
def labels(N):return{cusp(N,a,c)for c in range(N)for a in range(gcd(c,N))if gcd(a,gcd(c,N))==1}
def diamond(N,P,u):return cusp(N,pow(u,-1,N)*P[0],u*P[1])
def orbits(values,act):
 left=set(values);out=[]
 while left:
  x=next(iter(left));first=x;part=set()
  while x not in part:part.add(x);x=act(x)
  assert x==first and part<=left;left-=part;out.append(frozenset(part))
 return out
L=labels(69);assert len(L)==88 and all(diamond(69,p,22)!=p for p in L)
quot=orbits(L,lambda p:diamond(69,p,22));assert len(quot)==44 and all(len(o)==2 for o in quot)
Delta=[1,22,47,68];assert all((u*u+1)%69!=0 and(u*u-u+1)%69!=0 for u in Delta)
source_index=2112;quot_index=source_index//2;assert 1+source_index//12-len(L)//2==133 and 1+quot_index//12-len(quot)//2==67
L23=labels(23);assert len(L23)==22 and 1+264//12-len(L23)//2==12 and quot_index//264==4
q=1024;a=[2,59]
for n in range(2,6):a.append(59*a[-1]-q*a[-2])
orders={1:q+1-a[1],5:q**5+1-a[5]}
for r,order in orders.items():
 v=0;n=order
 while n%3==0:v+=1;n//=3
 assert v==1 and pow(q,r,3)==1
p=H/'proof_packet_v2.m';s=p.read_text();assert s.count('Q67_PACKET_V2')==1
packet=ast.literal_eval(s.replace('Q67_PACKET_V2',"'Q67_PACKET_V2'",1).replace('[*','[').replace('*]',']'))
assert packet[0]=='Q67_PACKET_V2' and len(packet)==11
assert len(packet[2])-1==4 and packet[10][:3]==[69,10,2]

def degdiv(rows):
 for place,e in rows:assert place[0]in[0,1]and len(place[2])==2 and isinstance(e,int)
 return sum(place[1]*e for place,e in rows)
assert [degdiv(packet[i])for i in range(4,8)]==[10,10,60,0]
orb=packet[8];assert len(orb)==2
assert sorted([r[1]for r in orb])==[10,50]
assert all(r[1]==r[2]==r[3] and 10*r[4]==r[1]*r[5]for r in orb)
assert sum(r[4]for r in orb)==6 and sum(r[1]*r[5]for r in orb)==60
assert len(packet[9])==4
for base,cover,d,obs,weight in packet[9]:assert base[1]==cover[1]==d and len(base[2])==len(cover[2])==2
out={'status':'PASS','source_level':69,'source_genus':133,'source_modular_index':2112,'source_cusps':88,'quotient_diamond':22,'diamond_action_free_on_all_cusps':True,'quotient_cusps':44,'quotient_modular_index':1056,'quotient_genus':67,'base_X1_23_genus':12,'base_modular_index':264,'relative_degree':4,'3_adic_order_of_finite_E_orders':{'r1':1,'r5':1},'finite_E_orders':orders,'quartic_local_factor_degrees':[1,3],'packet_model_K_degree':len(packet[1])-1,'packet_model_Z_relative_degree':len(packet[2])-1,'packet_prime_generators_per_place':2,'packet_divisor_degrees':[10,10,60,0],'isogeny_orbits':orb,'packet_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'RR_status':'not asserted by this geometry audit'}
(H/'geometry_verified.json').write_text(json.dumps(out,indent=2)+'\n');print('Q67_GEOMETRY_PACKET_AUDIT_COMPLETE',json.dumps(out,sort_keys=True))
