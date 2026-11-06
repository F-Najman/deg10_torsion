179,1
S,ProfileTimes,"Lists vertices in profile graph in order of time. You need to SetProfile(true), run something, then call this (which will SetProfile(false) before dumping)",0,0,0,0,0,0,0,SeqEnum,Null,Null,Null,Null,Null
S,PrintProfile,"Lists vertices in profile graph in order of time. You need to SetProfile(true), run something, then call this (which will SetProfile(false) before dumping)",0,0,0,0,1,0,0,Null,Null,Null,Null,Null,Null
S,Factorization,The prime factorization of the rational number r,0,1,0,0,0,0,0,0,0,FldRatElt,,SeqEnum,Null,Null,Null,Null,Null
S,GSp,The group of symplectic similitudes of degree d over F_q,0,2,0,0,0,0,0,0,0,RngIntElt,,0,0,RngIntElt,,GrpMat,Null,Null,Null,Null,Null
S,Eltseq,Sorted sequence of tuples representing a multiset of integers,1,0,1,SetMulti,0,RngIntElt,1,0,0,0,0,0,0,0,SetMulti,,SeqEnum,Null,Null,Null,Null,Null
S,ReplaceCharacter,Efficiently replace every occurrence of the character c in s with the string d (c must be a single character but d need not be),0,3,0,0,0,0,0,0,0,MonStgElt,,0,0,MonStgElt,,0,0,MonStgElt,,MonStgElt,Null,Null,Null,Null,Null
S,ReplaceString,Greedily replace each occurrence of string c in s with the string d,0,3,0,0,0,0,0,0,0,MonStgElt,,0,0,MonStgElt,,0,0,MonStgElt,,MonStgElt,Null,Null,Null,Null,Null
S,djb2,Returns the b-bit djb2 hash of s. Default value of b is 64,0,1,0,0,0,0,0,0,0,MonStgElt,,RngIntElt,Null,Null,Null,Null,Null
S,PySplit,"Splits using Python semantics (different when #sep > 1, and different when sep at beginning or end)",0,2,0,0,0,0,0,0,0,MonStgElt,,0,0,MonStgElt,,SeqEnum,Null,Null,Null,Null,Null
S,split,"Splits the string s using the delimter d, including empty and trailing elements (equivalent to python r.split(d) in python)",0,2,0,0,0,0,0,0,0,MonStgElt,,0,0,MonStgElt,,SeqEnum,Null,Null,Null,Null,Null
S,getrecs,"Reads a delimited file, returns list of lists of strings (one list per line)",0,1,0,0,0,0,0,0,0,MonStgElt,,SeqEnum,Null,Null,Null,Null,Null
S,putrecs,"Given a list of lists of strings, creates a colon delimited file with one list per line",1,1,1,SeqEnum,1,SeqEnum,0,MonStgElt,2,0,0,1,0,0,0,0,SeqEnum,,0,0,MonStgElt,,Null,Null,Null,Null,Null,Null
S,StringToStrings,"Given a string encoding a list of strings that do not contain commas or whitespace, e.g. ""[cat,dog]"", returns a list of the strings, e.g [ ""cat"", ""dog"" ]",0,1,0,0,0,0,0,0,0,MonStgElt,,SeqEnum,Null,Null,Null,Null,Null
S,sum,Sum of a sequence (empty sum is 0),0,1,0,0,0,0,0,0,0,SeqEnum,,Any,Null,Null,Null,Null,Null
S,sum,Sum of a vector,0,1,0,0,0,0,0,0,0,ModTupRngElt,,Any,Null,Null,Null,Null,Null
S,prod,Product of a sequence (empty product is 1),0,1,0,0,0,0,0,0,0,SeqEnum,,Any,Null,Null,Null,Null,Null
S,prod,Product of a vector,0,1,0,0,0,0,0,0,0,ModTupRngElt,,Any,Null,Null,Null,Null,Null
S,strip,"Strips spaces and carraige returns from string; used to be faster than StripWhiteSpace, now that StripWhiteSpace is faster we just call it",0,1,0,0,0,0,0,0,0,MonStgElt,,MonStgElt,Null,Null,Null,Null,Null
S,sprint,Sprints object X with spaces and carraige returns stripped,0,1,0,0,0,0,0,0,0,Any,,MonStgElt,Null,Null,Null,Null,Null
S,atoi,Converts a string to an integer (alias for StringToInteger),0,1,0,0,0,0,0,0,0,MonStgElt,,RngIntElt,Null,Null,Null,Null,Null
S,itoa,"Converts an integer to a string (equivalent to but slightly slower than IntegerToString, faster than sprint)",0,1,0,0,0,0,0,0,0,RngIntElt,,MonStgElt,Null,Null,Null,Null,Null
S,StringToReal,Converts a decimal string (like 123.456 or 1.23456e40 or 1.23456e-10) to a real number at default precision,0,1,0,0,0,0,0,0,0,MonStgElt,,RngIntElt,Null,Null,Null,Null,Null
S,atof,"Converts a decimal string (like ""123.456"") to a real number at default precision",0,1,0,0,0,0,0,0,0,MonStgElt,,RngIntElt,Null,Null,Null,Null,Null
S,StringsToAssociativeArray,"Converts a list of strings ""a=b"" to an associative array A with string keys and values such that A[a]=b. Ignores strings not of the form ""a=b""",1,0,1,SeqEnum,0,MonStgElt,1,0,0,0,0,0,0,0,SeqEnum,,Assoc,Null,Null,Null,Null,Null
S,atod,"Converts a list of strings ""a=b"" to an associative array A with string keys and values such that A[a]=b (alias for StringsToAssociativeArray)",1,0,1,SeqEnum,0,MonStgElt,1,0,0,0,0,0,0,0,SeqEnum,,Assoc,Null,Null,Null,Null,Null
S,StringToIntegerArray,"Given string representing a sequence of integers, returns the sequence (faster and safer than eval)",0,1,0,0,0,0,0,0,0,MonStgElt,,SeqEnum,Null,Null,Null,Null,Null
S,atoii,Converts a string to a sequence of integers (alias for StringToIntegerArray),0,1,0,0,0,0,0,0,0,MonStgElt,,SeqEnum,Null,Null,Null,Null,Null
S,iitoa,Converts a sequence of integers to a string (faster than sprint),1,0,1,SeqEnum,0,RngIntElt,1,0,0,0,0,0,0,0,SeqEnum,,MonStgElt,Null,Null,Null,Null,Null
S,atoiii,Converts a string to a sequence of sequences of integers,0,1,0,0,0,0,0,0,0,MonStgElt,,SeqEnum,Null,Null,Null,Null,Null
S,atoiiii,Converts a string to a sequence of sequences of sequences of integers,0,1,0,0,0,0,0,0,0,MonStgElt,,SeqEnum,Null,Null,Null,Null,Null
S,StringToRationalArray,"Given string representing a sequence of rational numbers, returns the sequence (faster and safer than eval)",0,1,0,0,0,0,0,0,0,MonStgElt,,SeqEnum,Null,Null,Null,Null,Null
S,StringToRealArray,"Given string representing a sequence of real numbers, returns the sequence (faster and safer than eval)",0,1,0,0,0,0,0,0,0,MonStgElt,,SeqEnum,Null,Null,Null,Null,Null
S,atoff,Converts a string to a sequence of reals (alias for StringToRealArray),0,1,0,0,0,0,0,0,0,MonStgElt,,SeqEnum,Null,Null,Null,Null,Null
S,atofff,Converts a string to a sequence of sequences of real numbers,0,1,0,0,0,0,0,0,0,MonStgElt,,SeqEnum,Null,Null,Null,Null,Null
S,goodp,Returns true if the specified polynomial is square free modulo p (without computing the discrimnant of f),0,2,0,0,0,0,0,0,0,RngIntElt,,0,0,RngUPolElt,,RngIntElt,Null,Null,Null,Null,Null
S,Base26Encode,"Given a nonnegative integer n, returns its encoding in base-26 (a=0,..., z=25, ba=26,...)",0,1,0,0,0,0,0,0,0,RngIntElt,,MonStgElt,Null,Null,Null,Null,Null
S,Base26Decode,"Given string representing a nonnegative integer in base-26 (a=0,..., z=25, ba=26,...) returns the integer",0,1,0,0,0,0,0,0,0,MonStgElt,,RngIntElt,Null,Null,Null,Null,Null
S,PolycyclicPresentation,"Given a sequence of generators in a uniquely represented hashable polycyclic group with composition law m and identity element id, returns sequence of relative orders and a map from group elements to exponent vectors",0,3,0,0,0,0,0,0,0,Any,,0,0,UserProgram,,0,0,SeqEnum,,SeqEnum,Map,Null,Null,Null,Null
S,OrderStats,Multiset of order statistics of elements of the group G,0,1,0,0,0,0,0,0,0,Grp,,SetMulti,Null,Null,Null,Null,Null
S,CyclicGenerators,A list of generators of the distinct cyclic subgroups of the finite abelian group G,0,1,0,0,0,0,0,0,0,GrpAb,,SeqEnum,Null,Null,Null,Null,Null
S,ConjugateIntersection,"Given finite subgroups H1 and H2 of a group G, returns a largest subgroup of G contained in a conjugate of H1 and a conjugate of H2",0,3,0,0,0,0,0,0,0,Grp,,0,0,Grp,,0,0,Grp,,Grp,Null,Null,Null,Null,Null
S,ConjugateCompositum,"Given subgroups H1 and H2 of G, returns a smallest subgroup of G that contains a conjugate of H1 and a conjugate of H2",0,3,0,0,0,0,0,0,0,Grp,,0,0,Grp,,0,0,Grp,,Grp,Null,Null,Null,Null,Null
S,IndexFibers,Given a list of objects S and a function f on S creates an associative array satisfying A[f(s)] = [t:t in S|f(t) eq f(s)],0,2,0,0,0,0,0,0,0,UserProgram,,0,0,SeqEnum,,Assoc,Null,Null,Null,Null,Null
S,IndexFibers,Given a list of objects S and a function f on S creates an associative array satisfying A[f(s)] = [t:t in S|f(t) eq f(s)],0,2,0,0,0,0,0,0,0,UserProgram,,0,0,List,,Assoc,Null,Null,Null,Null,Null
S,IndexFile,"Loads file with colon-delimited columns (or as specified by Delimiter) creating an AssociativeArray with key specified by keys (a column index or list of column indexes) and contents determined by data, which should be either a column index or list of column indexes (empty list means entire record)",0,2,0,0,0,0,0,0,0,Any,,0,0,MonStgElt,,Assoc,Null,Null,Null,Null,Null
S,ReduceToReps,Given a list of objects S and an equivalence relation E on S returns a maximal sublist of inequivalent objects,0,2,0,0,0,0,0,0,0,UserProgram,,0,0,SeqEnum,,SeqEnum,Null,Null,Null,Null,Null
S,Classify,Given a list of objects S and an equivalence relation E on them returns a list of equivalence classes (each of which is a list),0,2,0,0,0,0,0,0,0,UserProgram,,0,0,SeqEnum,,SeqEnum,Null,Null,Null,Null,Null
S,DihedralGroup,Construct the generalized dihedral group dih(G) := G semidirect Z/2Z with Z/2Z acting by inversion,0,1,0,0,0,0,0,0,0,GrpAb,,Grp,Null,Null,Null,Null,Null
S,Quotients,Returns a list of quotients of G (either all non-trivial quotients or those of the specified Order),0,1,0,0,0,0,0,0,0,Grp,,SeqEnum,Null,Null,Null,Null,Null
S,RandomizeForm,Applies a random invertible linear change of variables to the specified homogeneous polynomial (preserves integrality),0,1,0,0,0,0,0,0,0,RngMPolElt,,RngMPolElt,Null,Null,Null,Null,Null
S,RandomizeForms,Applies a random invertible linear change of variables to the specified sequence of homogeneous polynomials (preserves integrality),1,0,1,SeqEnum,0,RngMPolElt,1,0,0,0,0,0,0,0,SeqEnum,,SeqEnum,Null,Null,Null,Null,Null
S,MinimizeGenerators,Given a finite group G tries to reduce the number of generators by sampling random elements. Result is not guaranteed to minimize the number of generators but this is very likely,0,1,0,0,0,0,0,0,0,Grp,,Grp,Null,Null,Null,Null,Null
S,RegularRepresentation,The regular representation of H as a permutation group of degree #H,0,1,0,0,0,0,0,0,0,Grp,,GrpPerm,Null,Null,Null,Null,Null
S,HurwitzClassNumber,"The Hurwitz class number of positive definite binary quadratic forms of discriminant -N with each class C counted with multiplicity 1/#Aut(C), extended by Zagier to H(0)=-1/12 and H(-u^2)=-u/2, with H(-n) = 0 for all nonsquare n>0",0,1,0,0,0,0,0,0,0,RngIntElt,,FldRatElt,Null,Null,Null,Null,Null
S,KroneckerClassNumber,"The sum of the class numbers of the discriminants DD that divide the given imaginary quadratic discriminant D (this is not the same as the Hurwitz class number of -D, in particular, it is always an integer)",0,1,0,0,0,0,0,0,0,RngIntElt,,RngIntElt,Null,Null,Null,Null,Null
S,Log,"Given a,b in (Z/nZ)* returns least nonnegative x such that a^x = b or -1 if no such x exists",0,2,0,0,0,0,0,0,0,RngIntResElt,,0,0,RngIntResElt,,RngIntElt,Null,Null,Null,Null,Null
S,C4C6Invariants,Returns the c4 and c6 invariants of the specified elliptic curve E/Q (assumes E is defined by an integral model),1,0,1,CrvEll,0,FldRat,1,0,0,0,0,0,0,0,CrvEll,,RngInt,RngInt,Null,Null,Null,Null
S,GetFilenames,Return the filenames where such intrinsics are defined,0,1,0,0,0,0,0,0,0,Intrinsic,,SeqEnum,Null,Null,Null,Null,Null
S,WriteStderr,write to stderr,0,1,0,0,1,0,0,0,0,MonStgElt,,Null,Null,Null,Null,Null,Null
S,WriteStderr,write to stderr,0,1,0,0,1,0,0,0,0,Err,,Null,Null,Null,Null,Null,Null
S,NormalizedDixmierOhnoInvariants,"Normalized Dixmier-Ohno invaraints of smooth plane quartic f(x,y,z)=0 defined over Q",0,1,0,0,0,0,0,0,0,RngMPolElt,,SeqEnum,Null,Null,Null,Null,Null
S,NormalizedShiodaInvariants,"Normalized Shioda invariants of genus 3 hyperelliptic curve (the invariants must lie in Q, but the curve need not be defined over Q)",0,1,0,0,0,0,0,0,0,CrvHyp,,SeqEnum,Null,Null,Null,Null,Null
S,NormalizedShiodaInvariants,"Normalized Shioda invariants of genus 3 hyperelliptic curve (the invariants must lie in Q, but the curve need not be defined over Q)",0,2,0,0,0,0,0,0,0,RngUPolElt,,0,0,RngUPolElt,,SeqEnum,Null,Null,Null,Null,Null
S,SPQInvariants,"Normalized Dixmier-Ohno invaraints of smooth plane quartic f(x,y,z)=0 defined over Q",0,1,0,0,0,0,0,0,0,RngMPolElt,,SeqEnum,Null,Null,Null,Null,Null
S,SPQInvariants,"Normalized Dixmier-Ohno invaraints of smooth plane quartic f(x,y,z)=0 defined over Q",0,1,0,0,0,0,0,0,0,MonStgElt,,SeqEnum,Null,Null,Null,Null,Null
S,SPQIsIsomorphic,"Tests isomorphism of smooth plane curves f(x,y,z)=0 by computing a matrix M in GL(3,F) such that f1^M is a scalar multiple of f2",0,2,0,0,0,0,0,0,0,RngMPolElt,,0,0,RngMPolElt,,BoolElt,GrpMatElt,Null,Null,Null,Null
S,MonicQuadraticRoots,Returns the complete list of solutions to x^2+bx+c = 0 in Z/p^eZ for p prime and e ge 1 (does not verify the primality of p),0,4,0,0,0,0,0,0,0,Any,,0,0,RngIntElt,,0,0,RngIntElt,,0,0,RngIntElt,,SeqEnum,Null,Null,Null,Null,Null
S,MonicQuadraticRoots,Returns the complete list of solutions to x^2+bx+c = 0 in Z/mZ,0,3,0,0,0,0,0,0,0,RngIntElt,,0,0,RngIntElt,,0,0,RngIntElt,,SeqEnum,Null,Null,Null,Null,Null
S,ChangeRing,Given f = sum a_i*x^i returns sum pi(a_i)*x^i,0,2,0,0,0,0,0,0,0,Map,,0,0,RngUPolElt,,RngUPolElt,Null,Null,Null,Null,Null
S,PrimePowers,"Ordered list of prime powers q <= B (complexity is O(B log(B) loglog(B)), which is suboptimal but much better than testing individual prime powers)",0,1,0,0,0,0,0,0,0,RngIntElt,,SeqEnum,Null,Null,Null,Null,Null
S,ProperDivisors,Sorted list of postive proper divisors of the integer N,0,1,0,0,0,0,0,0,0,RngIntElt,,SeqEnum,Null,Null,Null,Null,Null
S,PrimesInInterval,"Primes of K with norm in [min,max]",0,3,0,0,0,0,0,0,0,RngIntElt,,0,0,RngIntElt,,0,0,FldNum,,SeqEnum,Null,Null,Null,Null,Null
S,NumberOfRoots,The number of rational roots of the polynomial f,1,0,1,RngUPolElt,0,FldFin,1,0,0,0,0,0,0,0,RngUPolElt,,RngIntElt,Null,Null,Null,Null,Null
S,TracesToLPolynomial,"Given a sequence of Frobenius traces of a genus g curve over Fq,Fq^2,...,Fq^g returns the corresponding L-polynomial",1,0,1,SeqEnum,0,RngIntElt,2,0,0,0,0,0,0,0,RngIntElt,,0,0,SeqEnum,,RngUPolElt,Null,Null,Null,Null,Null
S,PointCountsToLPolynomial,"Given a sequence of point counts of a genus g curve over Fq,Fq^2,...,Fq^g returns the corresponding L-polynomial",1,0,1,SeqEnum,0,RngIntElt,2,0,0,0,0,0,0,0,RngIntElt,,0,0,SeqEnum,,RngUPolElt,Null,Null,Null,Null,Null
S,LPolynomialToTraces,"Returns the sequence of Frobenius traces for a genus g curve over Fq,Fq^2,...,Fq^g corresponding to the givien L-polynomial of degree 2g (or just over Fq,..Fq^d if d is specified)",0,1,0,0,0,0,0,0,0,RngUPolElt,,SeqEnum,RngIntElt,Null,Null,Null,Null
S,LPolynomialToPointCounts,"Returns the sequence of point counrs of a genus g curve over Fq,Fq^2,...,Fq^g corresponding to the givien L-polynomial of degree 2g (or just over Fq,..Fq^d if d is specified)",0,1,0,0,0,0,0,0,0,RngUPolElt,,SeqEnum,RngIntElt,Null,Null,Null,Null
