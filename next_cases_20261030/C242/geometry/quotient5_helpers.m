/* Accepted degree-independent projection, inverse place, and source lift. */
function InverseThree(coords,values)
 F:=Parent(coords[1]);K:=Parent(values[1]);k:=PrimeField(K);p:=#k;e:=LCM([Degree(MinimalPolynomial(v,k)):v in values]);found:=false;
 for bits in [[1,0,0],[0,1,0],[0,0,1]] cat [[1,a,b]:a,b in [0..p-1]] cat [[0,1,a]:a in [0..p-1]] do
  gen:=&+[K!bits[i]*values[i]:i in [1..3]];pol:=MinimalPolynomial(gen,k);if Degree(pol) eq e then chosen:=bits;found:=true;break;end if;
 end for;assert found;M:=Matrix(k,[Eltseq(gen^i,k):i in [0..e-1]]);assert Rank(M) eq e;
 cc:=[Eltseq(Solution(M,Vector(k,Eltseq(v,k)))):v in values];assert &and[&+[K!cc[j][i+1]*gen^i:i in [0..e-1]] eq values[j]:j in [1..3]];
 gg:=&+[F!chosen[i]*coords[i]:i in [1..3]];fs:=[Evaluate(pol,gg)] cat [coords[j]-&+[F!cc[j][i+1]*gg^i:i in [0..e-1]]:j in [1..3]];
 O:=MaximalOrderFinite(F);I:=ideal<O|fs>;fa:=Factorization(I);places:=[Place(t[1]):t in fa|t[2] gt 0];
 if #places eq 0 then Oi:=MaximalOrderInfinite(F);Ii:=ideal<Oi|fs>;fi:=Factorization(Ii);places:=[Place(t[1]):t in fi|t[2] gt 0];end if;
 assert #places eq 1 and Degree(places[1]) eq e;assert &and[f eq 0 or Valuation(f,places[1]) gt 0:f in fs];return places[1];
end function;

function Project42(E,P,Q)
 aa:=aInvariants(E);assert aa[1] eq 0 and aa[3] eq 0;
 S:=6*Q;T:=14*Q;R:=21*Q;Z:=P+R;assert 7*S eq E!0 and 3*T eq E!0;
 xs:=S[1];ys:=S[2];m:=(3*xs^2+2*aa[2]*xs+aa[4])/(2*ys);a2p:=3*xs+aa[2]-m^2;a3p:=2*ys;
 kv:=a3p/a2p;bv:=-a2p^3/a3p^2;cv:=1-2*m*a2p/a3p;v:=bv/cv;assert bv eq v^2*(v-1) and cv eq v*(v-1);H:=v^2-v+1;assert H ne 0;
 u:=4*(P[1]-R[1])*(P[1]-Z[1])*(R[1]-Z[1])/(kv^6*v^4*(v-1)^4);
 xx:=(T[1]+aa[2]/3)/(kv^2*H^2);yy:=T[2]/(kv^3*H^3);
 t:=u^2;D:=t^2+13*t+49;N:=t^2+5*t+1;M:=t^4+14*t^3+63*t^2+70*t-7;
 assert yy^2 eq xx^3-N/(48*D)*xx+M/(864*D^2);assert 7*Q ne E!0 and 42*Q eq E!0;
 return [u,xx,yy];
end function;
function Lift42(values)
 d:=Degree(Parent(values[1]));L:=GF(5^(3*d));u,x,y:=Explode([L!z:z in values]);t:=u^2;D:=t^2+13*t+49;N:=t^2+5*t+1;M:=t^4+14*t^3+63*t^2+70*t-7;
 EL:=EllipticCurve([L|0,0,0,-N/(48*D),M/(864*D^2)]);T:=EL![x,y];assert 3*T eq EL!0 and T ne EL!0;
 R<W>:=PolynomialRing(L);s:=t+8;vr:=Roots(W^3-s*W^2+(s-3)*W+1);assert #vr eq 3;v:=vr[1][1];bv:=v^2*(v-1);cv:=v*(v-1);H:=v^2-v+1;assert H ne 0;
 a2v:=-bv+(1-cv)^2/4;S:=EL![a2v/(3*H^2),-bv/(2*H^3)];assert Order(S) eq 7;
 tw:=Roots(W^3-N/(48*D)*W+M/(864*D^2));assert #tw eq 3;P:=EL![tw[1][1],0];R2:=EL![tw[2][1],0];Q:=R2+2*T+6*S;vv:=Project42(EL,P,Q);
 if vv[1] ne u then P,R2:=Explode([R2,P]);Q:=R2+2*T+6*S;vv:=Project42(EL,P,Q);end if;
 assert vv eq [u,x,y] and Order(Q) eq 42 and P ne 21*Q;
 return EL,P,Q;
end function;

