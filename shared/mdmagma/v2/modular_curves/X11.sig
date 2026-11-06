179,0
T,MDCrvMod11,-,1,MDCrvMod
A,MDCrvMod11,10,level,curve,base_ring,N,M,_E,_P,_Q,_coordinates,_classgroup
S,Print,Print X,0,1,0,0,1,0,0,0,0,MDCrvMod11,,Null,Null,Null,Null,Null,Null
S,MDX11,"Create the modular curve X_1(N, M)",0,3,0,0,0,0,0,0,0,Rng,,0,0,RngIntElt,,0,0,RngIntElt,,MDCrvMod11,Null,Null,Null,Null,Null
S,_InitMDCrvMod11,Initialize an MDCrvMod1 object,0,4,0,0,1,0,0,0,0,Rng,,0,0,RngIntElt,,0,0,RngIntElt,,0,0,MDCrvMod11,,Null,Null,Null,Null,Null,Null
S,LevelStructure,"Return the level structure corresponding to the place x on Y1(M,N). To be precise. It returns is a RecFormat r such that r`P is a point of order M and r`Q a point of order N. It raises an error if x is a cusp",0,2,0,0,0,0,0,0,0,PlcCrvElt,,0,0,MDCrvMod11,,Rec,Null,Null,Null,Null,Null
S,ModuliPoint,"Return the place x on X1(M,N), corresponding to the pair (E, levelstructure). It assumes that levelstructure is a RecFormat such that levelstructure`P is a point of order M and levelstructure`Q a point of order N. It only works if M = 2",0,3,0,0,0,0,0,0,0,Rec,,0,0,CrvEll,,0,0,MDCrvMod11,,PlcCrvElt,Null,Null,Null,Null,Null
S,ApplyIsogeny,Return the level structure corresponding to phi(levelstructure),1,1,2,MapSch,0,CrvEll,0,CrvEll,3,0,0,0,0,0,0,0,Rec,,0,0,MapSch,,0,0,MDCrvMod11,,Rec,Null,Null,Null,Null,Null
S,DiamondOperator,"Return the level structure corresponding to d*P, d*Q;",0,3,0,0,0,0,0,0,0,Rec,,0,0,RngIntElt,,0,0,MDCrvMod11,,Rec,Null,Null,Null,Null,Null
S,DegeneracyMap,"If the place x corresponds to (E,P,Q) the place y corresponding to (E,(i*P+Q)*Level(X)/Level(Y)) taking into account multiplicities. I.e. if deg(y) = deg(x) it returns y, if not y deg(x)/deg(y) is returned",0,3,0,0,0,0,0,0,0,PlcCrvElt,,0,0,MDCrvMod1,,0,0,MDCrvMod11,,PlcCrvElt,Null,Null,Null,Null,Null
S,DerickxNormalForm_bc,"Assumes that P is a 2 torsion point Return the b,c of the Derickx normal form of (E,P,Q). I.e.: E = [0,c,0,(1-b-c)*b,0] P = [0,0]; Q = [b,b];",0,3,0,0,0,0,0,0,0,PtEll,,0,0,PtEll,,0,0,CrvEll,,Seq,Null,Null,Null,Null,Null
S,_equation_X11,"Input: m,n - integers such that m divides n base_ring - a ring equation_directory - directory with files X1_m_n.txt containing models zeta_m - a primitive mth root of unity in the base_ring (if unspecified one will be chosen) Output: C - a curve Returns an algebraic model C of the modular curve X_1(m,n) as a curve over base_ring",0,3,0,0,0,0,0,0,0,Any,,0,0,Any,,0,0,Any,,CrvPln,Null,Null,Null,Null,Null
S,Cusps,The cusps on this modular curve,0,1,0,0,0,0,0,0,0,MDCrvMod11,,SeqEnum,Null,Null,Null,Null,Null
