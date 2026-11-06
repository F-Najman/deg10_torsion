#!/usr/bin/env python3
"""Independent GF(2^10) audit of the remaining full C69 orbit input."""
from pathlib import Path
import ast,json,time,hashlib
H=Path(__file__).resolve().parent;t0=time.monotonic()
def read(path):return ast.literal_eval(Path(path).read_text().replace('[*','[').replace('*]',']'))
r=read(H/'selected_input.m');basepath=H.parent/'degree10_base.m';jpath=H.parent/'base_j_functions.m';base=read(basepath);jc=read(jpath)[0]
assert r[:3]==[69,10,2] and base[0:2]==[20,2] and r[3]==base[2][1][0]
def bits(v):return sum(int(c)<<i for i,c in enumerate(v))
modulus=bits(r[4]);assert modulus.bit_length()==11

def prem(a,b):
 while a and a.bit_length()>=b.bit_length():a^=b<<(a.bit_length()-b.bit_length())
 return a

def pgcd(a,b):
 while b:a,b=b,prem(a,b)
 return a

def mul(a,b):
 z=0
 while b:
  if b&1:z^=a
  b>>=1;a<<=1
  if a&1024:a^=modulus
 return z

def power(a,n):
 z=1
 while n:
  if n&1:z=mul(z,a)
  a=mul(a,a);n>>=1
 return z

def inv(a):assert a;return power(a,1022)
def div(a,b):return mul(a,inv(b))
def pol(coeff,x):
 z=0
 for c in reversed(coeff):z=mul(z,x)^c
 return z

def rf(code,x):return div(pol(code[0],x),pol(code[1],x))
assert power(2,1024)==2 and pgcd(power(2,32)^2,modulus)==1 and pgcd(power(2,4)^2,modulus)==1

def jfromkey(key):
 x,y=bits(key[:10]),bits(key[10:20]);assert mul(y,y)^mul(mul(mul(x,x),x)^x^1,y)^mul(x,x)==0
 return rf(jc[0],x)^mul(rf(jc[1],x),y)
j=jfromkey(r[3]);jother=jfromkey(base[2][0][0]);assert pol(r[7],j)==0 and pol(r[7],jother)!=0
v=j;orb=[]
while v not in orb:orb.append(v);v=mul(v,v)
assert v==j and len(orb)==10
A=[bits(x)for x in r[5]];P=tuple(bits(x)for x in r[6]);a1,a2,a3,a4,a6=A

def oncurve(P):
 if P is None:return True
 x,y=P;return mul(y,y)^mul(mul(a1,x),y)^mul(a3,y)^mul(mul(x,x),x)^mul(a2,mul(x,x))^mul(a4,x)^a6==0

def add(P,Q):
 if P is None:return Q
 if Q is None:return P
 x,y=P;u,v=Q
 if x==u:
  if y!=v or(mul(a1,x)^a3)==0:return None
  lam=div(mul(x,x)^a4^mul(a1,y),mul(a1,x)^a3)
 else:lam=div(y^v,x^u)
 nu=y^mul(lam,x);xx=mul(lam,lam)^mul(a1,lam)^a2^x^u;yy=mul(lam^a1,xx)^nu^a3
 out=(xx,yy);assert oncurve(out);return out

def multiple(k,P):
 z=None
 while k:
  if k&1:z=add(z,P)
  P=add(P,P);k>>=1
 return z
assert oncurve(P) and multiple(69,P)is None and multiple(23,P)is not None and multiple(3,P)is not None
b2=mul(a1,a1);b4=mul(a1,a3);b6=mul(a3,a3);b8=mul(b2,a6)^mul(mul(a1,a3),a4)^mul(a2,b6)^mul(a4,a4)
delta=mul(mul(b2,b2),b8)^mul(b6,b6)^mul(mul(b2,b4),b6);c4=mul(b2,b2);assert div(mul(mul(c4,c4),c4),delta)==j
count=1
for x in range(1024):
 k=mul(a1,x)^a3;rhs=mul(mul(x,x),x)^mul(a2,mul(x,x))^mul(a4,x)^a6
 if k==0:count+=1;continue
 z=div(rhs,mul(k,k));trace=0
 for i in range(10):trace^=z;z=mul(z,z)
 assert trace in [0,1]
 if trace==0:count+=2
assert count==966
result={'N':69,'base_case':2,'field_degree':10,'field_polynomial_irreducible':True,'selected_base_key_exact':True,'selected_j_minpoly':r[7],'j_orbit_degree':10,'distinct_from_base_case1_j_orbit':True,'curve_order':966,'point_order':69,'elliptic_curve_equation_verified':True,'j_formula_verified':True,'status':'PASS','elapsed_seconds':time.monotonic()-t0,'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in [H/'selected_input.m',basepath,jpath]}}
(H/'selected_input_verified.json').write_text(json.dumps(result,indent=2)+'\n');print('DIRECT69_SELECTED_INPUT_REPLAY_COMPLETE',json.dumps(result,sort_keys=True))
