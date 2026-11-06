179,0
S,MDPlaceCommonZeros,"Recover the moduli place by factoring the common-zero ideal, avoiding all unrelated zeros of the first coordinate equation",2,0,1,SeqEnum,0,FldFunFracSchElt,1,1,SeqEnum,0,FldElt,2,0,0,0,0,0,0,0,SeqEnum,,0,0,SeqEnum,,PlcCrvElt,Null,Null,Null,Null,Null
S,MDIsogeniesFastPoints,"For each Frobenius-orbit of cyclic subgroups G of order q of E (over the base field K of E), return <E/G, images of points, orbit size e>, with E/G and the images defined over the field of definition L of G",2,0,1,CrvEll,0,FldFin,1,1,SeqEnum,0,PtEll,3,0,0,0,0,0,0,0,RngIntElt,,0,0,SeqEnum,,0,0,CrvEll,,List,Null,Null,Null,Null,Null
S,ModuliPoint3,"Invert the q,t coordinates on the standard model X_1(3,N), with fixed Weil pairing",0,3,0,0,0,0,0,0,0,Rec,,0,0,CrvEll,,0,0,MDCrvMod11,,PlcCrvElt,Null,Null,Null,Null,Null
S,ModuliPoint4,"Invert the q,t coordinates on the standard model X_1(4,N), with fixed Weil pairing",0,3,0,0,0,0,0,0,0,Rec,,0,0,CrvEll,,0,0,MDCrvMod11,,PlcCrvElt,Null,Null,Null,Null,Null
S,ModuliPointFast,"Dispatch to the exact moduli inverse for M=2,3,4",0,3,0,0,0,0,0,0,0,Rec,,0,0,CrvEll,,0,0,MDCrvMod11,,PlcCrvElt,Null,Null,Null,Null,Null
S,DiamondOperator,"Diamond operator on X_1(M,N), M=2,3,4, using the corresponding moduli inverse",0,3,0,0,0,0,0,0,0,PlcCrvElt,,0,0,RngIntElt,,0,0,MDCrvMod11,,PlcCrvElt,Null,Null,Null,Null,Null
S,HeckeOperatorFast,"T_q(x) on X_1(M,N), M=2,3,4; require q=1 modM to preserve the Weil pairing",0,3,0,0,0,0,0,0,0,PlcCrvElt,,0,0,RngIntElt,,0,0,MDCrvMod11,,DivCrvElt,Null,Null,Null,Null,Null
S,HeckeOperatorFast,"T_q(D) on X_1(M,N), for divisors supported on noncuspidal places",0,3,0,0,0,0,0,0,0,DivCrvElt,,0,0,RngIntElt,,0,0,MDCrvMod11,,DivCrvElt,Null,Null,Null,Null,Null
S,HeckeOperatorQuotientFast,"Rational-quotient T_q for M=2,3,4, adjusting the P-generator by q^-1 to fix the Weil pairing",0,3,0,0,0,0,0,0,0,PlcCrvElt,,0,0,RngIntElt,,0,0,MDCrvMod11,,DivCrvElt,Null,Null,Null,Null,Null
S,HeckeOperatorQuotientFast,Rational-quotient T_q on divisors supported on noncuspidal places,0,3,0,0,0,0,0,0,0,DivCrvElt,,0,0,RngIntElt,,0,0,MDCrvMod11,,DivCrvElt,Null,Null,Null,Null,Null
S,PlacesUpToDiamondFast,"One representative per diamond orbit, using simultaneous isomorphism tests on P and Q",1,1,1,SeqEnum,0,PlcCrvElt,2,0,0,0,0,0,0,0,SeqEnum,,0,0,MDCrvMod11,,SeqEnum,Null,Null,Null,Null,Null
