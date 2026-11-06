179,0
T,MDCrvMod1,-,1,MDCrvMod
A,MDCrvMod1,10,level,curve,N,_E,_P,_coordinates,_coordinates_xy,_coordinates_F2F3,_cusp_data,_classgroup
S,Print,Print X,0,1,0,0,1,0,0,0,0,MDCrvMod1,,Null,Null,Null,Null,Null,Null
S,MDX1,Create the modular curve X_1(N),0,2,0,0,0,0,0,0,0,Rng,,0,0,RngIntElt,,MDCrvMod1,Null,Null,Null,Null,Null
S,_InitMDCrvMod1,Initialize an MDCrvMod1 object,0,3,0,0,1,0,0,0,0,Rng,,0,0,RngIntElt,,0,0,MDCrvMod1,,Null,Null,Null,Null,Null,Null
S,LevelStructure,Return the level structure corresponding to the place x on Y1(N). To be precise. It returns is a RecFormat r such that r`P is a point of order N. It raises an error if x is a cusp,0,2,0,0,0,0,0,0,0,PlcCrvElt,,0,0,MDCrvMod1,,Rec,Null,Null,Null,Null,Null
S,ModuliPoint,"Return the place x on X1(N), corresponding to the pair (E, levelstructure). It assumes that levelstructure is a RecFormat such that levelstructure`P is a point of order N",0,3,0,0,0,0,0,0,0,Rec,,0,0,CrvEll,,0,0,MDCrvMod1,,PlcCrvElt,Null,Null,Null,Null,Null
S,TateNormalForm_bc,"Return the b,c of the tate normal form of (E,P) as in equation (2) of http://arxiv.org/pdf/1307.5719v1.pdf",0,2,0,0,0,0,0,0,0,PtEll,,0,0,CrvEll,,Seq,Null,Null,Null,Null,Null
S,TateNormalForm_xy,"return the x,y of the tate normal form of (E,P) as in section 2.1 http://arxiv.org/pdf/1307.5719v1.pdf",0,2,0,0,0,0,0,0,0,PtEll,,0,0,CrvEll,,Seq,Null,Null,Null,Null,Null
S,ApplyIsogeny,Return the level structure corresponding to phi(levelstructure),1,1,2,MapSch,0,CrvEll,0,CrvEll,3,0,0,0,0,0,0,0,Rec,,0,0,MapSch,,0,0,MDCrvMod1,,Rec,Null,Null,Null,Null,Null
S,DiamondOperator,Return the level structure corresponding to d*P; where d is coprime to the level,0,3,0,0,0,0,0,0,0,Rec,,0,0,RngIntElt,,0,0,MDCrvMod1,,Rec,Null,Null,Null,Null,Null
S,DegeneracyMap,Return the level structure corresponding to d*P; where d is a divisor of the level,0,3,0,0,0,0,0,0,0,Rec,,0,0,RngIntElt,,0,0,MDCrvMod1,,Rec,Null,Null,Null,Null,Null
S,DegeneracyMap,"If the place x corresponds to (E,P) the place y corresponding to (E,P*Level(X)/Level(Y)) taking into account multiplicities. I.e. if deg(y) = deg(x) it returns y, if not y deg(x)/deg(y) is returned",0,3,0,0,0,0,0,0,0,PlcCrvElt,,0,0,MDCrvMod1,,0,0,MDCrvMod1,,PlcCrvElt,Null,Null,Null,Null,Null
S,_equation_X1,Input: n - an integer base_ring - a ring equation_directory - directory with files FFFc<n>.txt containing models Output: C - a curve Returns an algebraic model C of the modular curve X_1(n) as a curve over base_ring,0,2,0,0,0,0,0,0,0,Any,,0,0,Any,,CrvPln,Null,Null,Null,Null,Null
S,CongruenceSubgroup,The subgroup of PSL2 corresponding to this modular curve,0,1,0,0,0,0,0,0,0,MDCrvMod1,,GrpPSL2,Null,Null,Null,Null,Null
S,Cusps,The cusps on this modular curve,0,1,0,0,0,0,0,0,0,MDCrvMod1,,SeqEnum,Null,Null,Null,Null,Null
S,CuspSignature,Returns the valuation of F2 and F3 at the place C,0,2,0,0,0,0,0,0,0,PlcCrvElt,,0,0,MDCrvMod1,,SeqEnum,Null,Null,Null,Null,Null
S,CuspOrbitCountQ,Returns the number of Galois orbits of cusps of X_1(n) using 2.8(1) of https://arxiv.org/pdf/2007.13929,0,1,0,0,0,0,0,0,0,MDCrvMod1,,RngIntElt,Null,Null,Null,Null,Null
S,_InitialiseCuspData,"Compute the signatures of all cusps, and create a lookup table to get the cusps corresponding to a signature. It raises an error if the signature doesn't uniquely determine the galois orbit over Q",0,1,0,0,1,0,0,0,0,MDCrvMod1,,Null,Null,Null,Null,Null,Null
S,_CuspData,Return the associative array mapping cusp signatures to cusps,0,1,0,0,0,0,0,0,0,MDCrvMod1,,Assoc,Null,Null,Null,Null,Null
S,CuspsWithSignature,Returns the cusps of the given signature,1,1,1,SeqEnum,0,RngIntElt,2,0,0,0,0,0,0,0,SeqEnum,,0,0,MDCrvMod1,,SeqEnum,Null,Null,Null,Null,Null
S,IsCuspSignatureUnique,Returns true if and only if there is only one galois orbit of cusps per cusp signature,0,1,0,0,0,0,0,0,0,MDCrvMod1,,BoolElt,Null,Null,Null,Null,Null
S,CuspOrbitsQ,The galois orbits over Q of cusps on this modular curve as divisors. Note that when the modular curve is defined over a finite field of positive characteristic one still has an action of the galois group over Q,0,1,0,0,0,0,0,0,0,MDCrvMod1,,SeqEnum,Null,Null,Null,Null,Null
