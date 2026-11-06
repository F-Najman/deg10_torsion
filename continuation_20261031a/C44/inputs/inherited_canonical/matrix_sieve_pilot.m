SieveOutput:="work/degree9/cyclic/continuation_20261026/C44_canonical_ring/matrix_pilot.raw.jsonseq";
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(16*10^9);SetSeed(20261044);tt:=Cputime();
load "work/degree9/cyclic/continuation_20261026/C44_canonical_ring/multiplication_data.m";
load "work/degree9/cyclic/continuation_20261026/C44_canonical_ring/sieve_pilot_inputs.m";
assert #TensorBlocks eq 105 and #Std2Exponents eq 105 and #Std4Exponents eq 245 and #Std6Exponents eq 385;
PF:=PolynomialRing(GF(3));
function EvalJets(J)
 dd:=J[1];mm:=J[2];LF<a>:=ext<GF(3)|PF!J[3]>;assert Degree(LF) eq dd;PS<T>:=PowerSeriesRing(LF,mm);pp:=[LF!1];for k in [1..dd-1] do Append(~pp,pp[#pp]*a);end for;zs:=[];
 for row in J[4] do cc:=[];for v in row do assert #v eq dd;Append(~cc,&+[LF!v[k]*pp[k]:k in [1..dd]]);end for;assert #cc eq mm;Append(~zs,PS!cc+O(T^mm));end for;assert #zs eq 36;
 powers:=AssociativeArray();for i in [1..36] do powers[<i,0>]:=PS!1;for k in [1..4] do powers[<i,k>]:=powers[<i,k-1>]*zs[i];end for;end for;
 rows2:=[];rows4:=[];
 for ee in Std2Exponents do vv:=&*[powers[<i,ee[i]>]:i in [1..36]|ee[i] gt 0];Append(~rows2,&cat[Eltseq(Coefficient(vv,k),GF(3)):k in [0..mm-1]]);end for;
 for ee in Std4Exponents do vv:=&*[powers[<i,ee[i]>]:i in [1..36]|ee[i] gt 0];Append(~rows4,&cat[Eltseq(Coefficient(vv,k),GF(3)):k in [0..mm-1]]);end for;
 E2:=Matrix(GF(3),rows2);E4:=Matrix(GF(3),rows4);assert Nrows(E2) eq 105 and Ncols(E2) eq dd*mm and Nrows(E4) eq 245;return E2,E4;
end function;
out:=Open(SieveOutput,"w");fprintf out,"{\"type\":\"header\",\"level\":44,\"p\":3,\"cases\":%o,\"R2\":105,\"R4\":245,\"R6\":385}\n",#JetCases;Flush(out);survivors:=[];
for jj in JetCases do
 k:=jj[1];e2t:=[];e4t:=[];e2u:=[];for J in jj[2] do A,B:=EvalJets(J);Append(~e2t,A);Append(~e4t,B);end for;for J in jj[3] do A,B:=EvalJets(J);Append(~e2u,A);end for;
 ET2:=HorizontalJoin(e2t);ET4:=HorizontalJoin(e4t);EU2:=HorizontalJoin(e2u);assert Ncols(ET2) eq 54 and Ncols(ET4) eq 54 and Ncols(EU2) eq 54;
 assert Rank(ET2) eq 54 and Rank(ET4) eq 54 and Rank(EU2) eq 54;VF:=Nullspace(ET2);VH:=Nullspace(ET4);assert Dimension(VF) eq 51 and Dimension(VH) eq 191;
 rank0:=0;attempts:=0;
 while attempts lt 8 and rank0 lt 105 do
  attempts+:=1;F:=Random(VF);while IsZero(F) do F:=Random(VF);end while;H:=Random(VH);while IsZero(H) do H:=Random(VH);end while;assert IsZero(F*ET2) and IsZero(H*ET4);
  MF:=ZeroMatrix(GF(3),245,385);for i in [1..105] do if F[i] ne 0 then MF+:=F[i]*TensorBlocks[i];end if;end for;assert Rank(MF) eq 245;
  QQ:=Transpose(BasisMatrix(Nullspace(Transpose(MF))));assert Nrows(QQ) eq 385 and Ncols(QQ) eq 140 and IsZero(MF*QQ);
  MH:=Matrix(GF(3),[Eltseq(H*TensorBlocks[i]):i in [1..105]]);AA:=HorizontalJoin(EU2,MH*QQ);assert Nrows(AA) eq 105 and Ncols(AA) eq 194;rank0:=Rank(AA);
 end while;
 if rank0 lt 105 then Append(~survivors,k);end if;
 fprintf out,"{\"type\":\"case\",\"case\":%o,\"attempts\":%o,\"rank\":%o,\"F\":%o,\"H\":%o,\"ET2\":%o,\"ET4\":%o,\"EU2\":%o,\"quotient_projection\":%o,\"matrix\":%o}\n",k,attempts,rank0,[Integers()!c:c in Eltseq(F)],[Integers()!c:c in Eltseq(H)],[Integers()!c:c in Eltseq(ET2)],[Integers()!c:c in Eltseq(ET4)],[Integers()!c:c in Eltseq(EU2)],[Integers()!c:c in Eltseq(QQ)],[Integers()!c:c in Eltseq(AA)];Flush(out);
 printf "NONIC44_CANONICAL_MATRIX_CASE index=%o rank=%o dimension105=true attempts=%o time=%o\n",k,rank0,attempts,Cputime(tt);
end for;
fprintf out,"{\"type\":\"complete\",\"cases\":%o,\"survivors\":%o}\n",#JetCases,survivors;delete out;
printf "NONIC44_CANONICAL_MULTIPLICATION_SIEVE_COMPLETE cases=%o survivors=%o all_rank105=%o time=%o\n",#JetCases,survivors,#survivors eq 0,Cputime(tt);quit;
