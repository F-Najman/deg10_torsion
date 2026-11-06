179,1
S,polredabs,Computes a smallest canonical defining polynomial of the etale algebra Q[x]/(f(x)) using pari,0,1,0,0,0,0,0,0,0,SeqEnum,,SeqEnum,Null,Null,Null,Null,Null
S,polredabs,Computes a smallest canonical defining polynomial of the etale algebra Q[x]/(f(x)) using pari,0,1,0,0,0,0,0,0,0,RngUPolElt,,RngUPolElt,Null,Null,Null,Null,Null
S,polredabs,Given a number field returns the same number field defined using a canonical defining polynomial using pari,0,1,0,0,0,0,0,0,0,FldNum,,FldNum,Null,Null,Null,Null,Null
S,polredabs,Given the rational field returns it unchanged,0,1,0,0,0,0,0,0,0,FldRat,,FldRat,Null,Null,Null,Null,Null
S,polredbest,Computes a small (non-canonical) defining polynomial of the etale algebra Q[x]/(f(x)) using pari,0,1,0,0,0,0,0,0,0,SeqEnum,,SeqEnum,Null,Null,Null,Null,Null
S,polredbest,Computes a small (non-canonical) defining polynomial of the etale algebra Q[x]/(f(x)) using pari,0,1,0,0,0,0,0,0,0,RngUPolElt,,RngUPolElt,Null,Null,Null,Null,Null
S,polredbest,Given a number field returns the same number field defined using a small (non-canonical) defining polynomial using pari,0,1,0,0,0,0,0,0,0,FldNum,,FldNum,Null,Null,Null,Null,Null
S,PerfectPowerBase,Returns the least positive integer m such that m^e = n for some positive integer e (m=n if n is not a perfect power),0,1,0,0,0,0,0,0,0,RngIntElt,,RngIntElt,Null,Null,Null,Null,Null
S,IsPolredabsCandidate,Returns true if the polynomial looks like Polredabs can easily handle it,0,1,0,0,0,0,0,0,0,RngUPolElt,,BoolElt,Null,Null,Null,Null,Null
S,IsPolredabsCandidate,Returns true if the polynomial looks like Polredabs can easily handle it,0,1,0,0,0,0,0,0,0,SeqEnum,,BoolElt,Null,Null,Null,Null,Null
S,polredbestify,Call polredbest repeatedly to get a smaller defining polynomial for the etale algebra Q[x]/(f(x)) using pari,0,1,0,0,0,0,0,0,0,RngUPolElt,,RngUPolElt,BoolElt,Null,Null,Null,Null
S,polredbestify,Call polredbest repeatedly to get a smaller defining polynomial for the etale algebra Q[x]/(f(x)) using pari,0,1,0,0,0,0,0,0,0,SeqEnum,,SeqEnum,BoolElt,Null,Null,Null,Null
S,polredbestwithroot,"Returns small polynomial as in Polredbest together with a root, using pari",0,1,0,0,0,0,0,0,0,RngUPolElt,,RngUPolElt,SeqEnum,Null,Null,Null,Null
S,polredabswithroot,"Returns Polredabs(f) together with a root, using pari",0,1,0,0,0,0,0,0,0,RngUPolElt,,RngUPolElt,SeqEnum,Null,Null,Null,Null
S,polredbestifywithroot,"Returns small polynomial as in Polredbestify together with a root, using pari. Will use polredabs if this seems feasible (this is indicated by the third return value)",0,1,0,0,0,0,0,0,0,RngUPolElt,,RngUPolElt,SeqEnum,BoolElt,Null,Null,Null
S,nfisincl,Returns a list of polynomials defining embeddings of the number field K defined by f in the number field L defined by g (each is specified by a polynomial h for which h(L.1) is a generator for the embedding of K in L),0,2,0,0,0,0,0,0,0,RngUPolElt,,0,0,RngUPolElt,,SeqEnum,Null,Null,Null,Null,Null
S,nfisincl,Returns a list of embeddings of the number field K defined by f in the number field L defined by g (each is specified by the coefficient sequence of a polynomial h for which h(L.1) is a generator for the embedding of K in L),0,2,0,0,0,0,0,0,0,SeqEnum,,0,0,SeqEnum,,SeqEnum,Null,Null,Null,Null,Null
