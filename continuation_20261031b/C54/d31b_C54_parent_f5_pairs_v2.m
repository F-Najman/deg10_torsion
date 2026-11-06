load "work/deg10/continuation_20261031b/C54/d31b_C54_parent_f5_setup.m";fib:=eval Read(root cat "class_members_all.m");src:=eval Read(root cat "member_two_isogenies_all.m");assert #fib[1]eq 192 and #src[1]eq 334 and #src[2]eq 192 and #src[3]eq 147;imidx:=AssociativeArray();cl:=AssociativeArray();for row in src[1]do id:=row[1];ky:=row[2];assert not IsDefined(imidx,ky);imidx[ky]:=id;p:=FromKey(ky);cl[id]:=phi(Divisor(p)-Degree(p)*A);if id mod 50 eq 0 then print "D31B_C54_PARENT_F5_TARGET_PROGRESS",id,Cputime(tt);end if;end for;
branches:=AssociativeArray();nker:=0;sourceObs:=[* *];
for row in fib[1]do id,d,ky:=Explode(row);KK:=ext<k|PK!ky[4..4+d]>;gg:=KK.1;bb:=Evaluate(PK!ky[5+d..4+2*d],gg);cc:=Evaluate(PK!ky[5+2*d..4+3*d],gg);EE:=EllipticCurve([KK|1-cc,-bb,-bb,0,0]);QQ:=EE![0,0];assert 27*QQ eq EE!0 and 9*QQ ne EE!0;fac:=Factorization(DivisionPolynomial(EE,2));assert &+[Degree(v[1])*v[2]:v in fac]eq 3 and &and[v[2]eq 1:v in fac];pbs:=[];
 for v in fac do f:=v[1];e:=Degree(f);if e eq 1 then LL:=KK;r0:=-Coefficient(f,0)/Coefficient(f,1);else LL:=ext<KK|f>;r0:=LL.1;end if;EL:=BaseExtend(EE,LL);P:=EL![0,0];aa:=aInvariants(EL);T:=EL![r0,-(aa[1]*r0+aa[3])/2];assert T ne EL!0 and 2*T eq EL!0;RP<z>:=PolynomialRing(LL);ET,ph:=IsogenyFromKernel(EL,z-r0);QT:=ph(P);assert 27*QT eq ET!0 and 9*QT ne ET!0;assert QT[1]eq P[1]+(P+T)[1]-T[1]and QT[2]eq P[2]+(P+T)[2]-T[2];kys:=Key(TateNormalForm_bc(EL,P+T));assert kys[1]eq d*e;ky2:=Key(TateNormalForm_bc(ET,QT));assert IsDefined(imidx,ky2)and d*e mod ky2[1]eq 0;Append(~pbs,[e,imidx[ky2],d*e div ky2[1]]);nker+:=1;
 end for;assert src[2][id][1]eq id and src[2][id][2]eq d;assert Sort(pbs)eq Sort([[v[2],v[3],v[4]]:v in src[2][id][3]]);branches[id]:=pbs;Append(~sourceObs,[*id,d,pbs*]);
end for;assert nker eq 361;
procedure Local(j,left,bs,current,~out)
 if j gt #bs then if left eq 0 then Append(~out,current);end if;return;end if;
 for n in[0..left div bs[j][1]]do w:=current;if n gt 0 then Append(~w,[bs[j][2],n*bs[j][3]]);end if;Local(j+1,left-n*bs[j][1],bs,w,~out);end for;
end procedure;
procedure Global(j,opts,current,~out)
 if j gt #opts then Append(~out,current);return;end if;
 for w in opts[j]do v:=current;for z in w do v[z[1]]+:=z[2];end for;Global(j+1,opts,v,~out);end for;
end procedure;
allowed:=AssociativeArray();for r in dy[8]do x:=Gold!r[7][1..3];y:=cv(Gold!r[8][1..3]);if IsDefined(allowed,x)then Include(~allowed[x],y);else allowed[x]:={y};end if;end for;cusp:=AssociativeArray();cusp[0]:={G!0};cusp[1]:={phi(FunctionFieldDivisor(Divisor(p))-A):p in cs|Degree(p)eq 1};assert #cusp[1]eq 9;expected:=AssociativeArray();for r in src[3]do ky:=[r[2]]cat &cat [[Integers()!z:z in w]:w in r[6]];if IsDefined(expected,ky)then expected[ky]+:=1;else expected[ky]:=1;end if;end for;actual:=AssociativeArray();out:=[* *];cdCounts:=[0,0];
for r in fib[2]do mi,ci,xi,nc,cd:=Explode(r);opts:=[* *];for term in nc do loc:=[];Local(1,term[2],branches[term[1]],[],~loc);Append(~opts,loc);end for;vv:=[];Global(1,opts,[0:i in[1..#src[1]]],~vv);
 for v in vv do sparse:=[[j,v[j]]:j in[1..#v]|v[j]gt 0];ky:=[mi]cat &cat sparse;if IsDefined(actual,ky)then actual[ky]+:=1;else actual[ky]:=1;end if;assert cd in[0,1];cdCounts[cd+1]+:=1;g:=&+[G|v[j]*cl[j]:j in[1..#v]|v[j]ne 0];possible:=(allowed[Gold!xi[1..3]]meet{g+w:w in cusp[cd]});assert #possible eq 0;Append(~out,[*mi,ci,cd,sparse,Eltseq(g),false*]);
 end for;
end for;assert Keys(actual)eq Keys(expected)and &and[actual[ky]eq expected[ky]:ky in Keys(expected)];assert #out eq 147 and cdCounts eq[59,88];Write(root cat "parent_f5_pairs.m",Sprint([*sourceObs,out,[nker,#src[1],#out,0] *]):Overwrite:=true);print "D31B_C54_PARENT_F5_PAIRS_COMPLETE",nker,#src[1],#out,cdCounts,0,Cputime(tt);quit;
