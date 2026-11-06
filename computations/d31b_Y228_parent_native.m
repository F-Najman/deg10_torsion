SetColumns(0);SetNthreads(1);SetQuitOnError(true);SetMemoryLimit(8*10^9);Alarm(900);tt:=Cputime();root:="work/deg10/continuation_20261031b/Y228_shared/";load "work/deg10/continuation_20261031b/Y228_shared/fast_primitive_inverse.m";load "work/deg10/continuation_20261031b/Y228_shared/union_cases.m";md:=eval Read(root cat "model_v5_s0.m");k:=GF(3);PK:=PolynomialRing(k);B<t>:=FunctionField(k);R<Y>:=PolynomialRing(B);
function DR(a)return B!(PK!a[1])/(B!(PK!a[2]));end function;
f:=R![DR(a):a in md[1]];assert IsIrreducible(f);F<u>:=FunctionField(f);assert Genus(F)eq 25;OF:=MaximalOrderFinite(F);OI:=MaximalOrderInfinite(F);
function DF(v)return F![DR(a):a in v];end function;coords:=[DF(v):v in md[2]];
function Prime(pd)O:=pd[2]eq 1 select OF else OI;gs:=[DF(v):v in pd[3]];assert #gs eq 2;I:=ideal<O|gs>;assert IsIntegral(I)and IsPrime(I);P:=Place(I);assert Degree(P)eq pd[1]and Ideal(P)eq I;return P;end function;
function Key(P)return Y228PrimitiveJointKey([Evaluate(v,P):v in coords]);end function;
cs:=[Prime(pd):pd in md[3]];assert &and[Degree(p)eq 3:p in cs];fields:=AssociativeArray();for row in UnionFields do fields[row[1]]:=ext<k|PK!row[2]>;end for;
allrows:=[* *];for shard in[0..5]do
 assert eval Read(root cat Sprintf("model_v5_s%o.m",shard))eq md;rows:=eval Read(root cat Sprintf("A5_v5_s%o.m",shard));allrows cat:=rows;
end for;assert #allrows eq 45;out:=[* *];np:=0;nr:=0;
for row in allrows do
 id,cusp,obs,cor,hs:=Explode(row);assert UnionCases[id][2]eq cusp;D:=DivisorGroup(F)!0;bound:=[* *];
 for ob in obs do
  e,mul,ky,ps,pd,hec:=Explode(ob);P:=Prime(ps);V:=Prime(pd);np+:=2;assert Degree(P)eq e and Degree(V)eq e;K:=fields[e];bb:=Evaluate(PK!ky[1..e],K.1);cc:=Evaluate(PK!ky[e+1..2*e],K.1);E:=EllipticCurve([K|0,cc,0,bb*(1-bb-cc),0]);P2:=E![0,0];Q:=E![bb,bb];QD:=5*Q;r:=P2[1];xx:=QD[1]-r;yy:=QD[2];dd:=[xx^3/yy^2,(cc+3*r)*xx^2/yy^2];assert Key(P)eq Y228PrimitiveJointKey([bb,cc])and Key(V)eq Y228PrimitiveJointKey(dd);
  DD:=-Divisor(P)-5*Divisor(V);he:=[* *];assert &+[h[1]:h in hec]eq 6;
  for h in hec do T:=Prime(h[3]);np+:=1;assert Degree(T)*h[2]eq e*h[1];DD+:=h[2]*Divisor(T);Append(~he,[*h[1],h[2],Key(T)*]);end for;
  assert Degree(DD)eq 0;D+:=mul*DD;Append(~bound,[*e,mul,ky,Key(P),Key(V),he*]);
 end for; hh:=[];for v in cor do A:=D+&+[v[j]*Divisor(cs[j]):j in[1..4]];assert Degree(A)eq 0;Append(~hh,#Basis(A:Reduction:=false,Simplification:="None"));nr+:=1;end for;assert hh eq hs;Append(~out,[*id,cusp,bound,cor,hh*]);print "Y228_PARENT_NATIVE_CASE",id,hh,Cputime(tt);
end for;
nei:=eval Read(root cat "m428_neighbor_A5.m");no:=[* *];for row in nei do if row[2]ne 1 then continue;end if;
 K:=fields[10];E:=EllipticCurve([K|Evaluate(PK!v,K.1):v in row[4]]);P2:=E![Evaluate(PK!v,K.1):v in row[5]];Q:=E![Evaluate(PK!v,K.1):v in row[6]];P:=Prime(row[7]);V:=Prime(row[9]);np+:=2;xx:=Q[1]-P2[1];yy:=Q[2];vals:=[xx^3/yy^2,(aInvariants(E)[2]+3*P2[1])*xx^2/yy^2];QD:=5*Q;xx:=QD[1]-P2[1];yy:=QD[2];dia:=[xx^3/yy^2,(aInvariants(E)[2]+3*P2[1])*xx^2/yy^2];assert Key(P)eq Y228PrimitiveJointKey(vals)and Key(V)eq Y228PrimitiveJointKey(dia);assert Degree(P)*row[8]eq 10 and Degree(V)*row[10]eq 50;D:=-row[8]*Divisor(P)-row[10]*Divisor(V);he:=[* *];assert &+[h[1]:h in row[11]]eq 6;
 for h in row[11]do T:=Prime(h[3]);np+:=1;assert h[2]*Degree(T)eq 10*h[1];D+:=h[2]*Divisor(T);Append(~he,[*h[1],h[2],Key(T)*]);end for;assert Degree(D)eq 0;h0:=#Basis(D:Reduction:=false,Simplification:="None");nr+:=1;assert h0 eq row[12]and h0 eq 0;Append(~no,[*row[1],row[2],Key(P),Key(V),he,h0*]);print "Y228_PARENT_NATIVE_NEIGHBOR",row[1],h0,Cputime(tt);
end for;assert #no eq 4;Write(root cat "parent_native_all45_and_neighbors.m",Sprint([*out,no,np,nr*]):Overwrite:=true);print "D31B_Y228_PARENT_NATIVE_COMPLETE",#out,#no,np,nr,Cputime(tt);quit;
