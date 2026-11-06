#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdint>
#include <cstdlib>
using namespace std;
array<uint64_t,10> columns; vector<array<uint64_t,10>> jets; vector<int> mults; unsigned long long seen=0; bool failed=false;
void leaf(){array<uint64_t,38> basis{};int rank=0;for(uint64_t v:columns){while(v){int p=__builtin_ctzll(v);if(basis[p])v^=basis[p];else{basis[p]=v;++rank;break;}}}++seen;if(rank!=10){failed=true;cerr<<"DEPENDENT rank "<<rank<<" multiplicities ";for(int m:mults)cerr<<m<<" ";cerr<<"\n";}}
void enumerate(int i,int left,int pos){if(failed)return;if(i==(int)jets.size()){if(left==0){if(pos!=10)abort();leaf();}return;}int lo=i+1==(int)jets.size()?left:0;for(int m=lo;m<=left;++m){mults[i]=m;for(int j=0;j<m;++j)columns[pos+j]=jets[i][j];enumerate(i+1,left-m,pos+m);}}
int main(int argc,char**argv){if(argc!=2)return 2;ifstream f(argv[1]);int n,d,m;f>>n>>d>>m;if(n!=18||d!=38||m!=10)return 3;jets.resize(n);mults.resize(n);uint64_t mask=(1ULL<<38)-1;for(auto&row:jets)for(auto&v:row){f>>v;v^=(v<<1)&mask;}array<uint64_t,6>s;for(auto&v:s){f>>v;v^=(v<<1)&mask;}if(!f)return 4;seen=0;enumerate(0,10,0);if(failed||seen!=8436285ULL)return 5;auto rational=seen;seen=0;for(int j=0;j<6;++j)columns[j]=s[j];enumerate(0,4,6);if(failed||seen!=5985ULL)return 6;cout<<"{\"status\":\"PASS\",\"rational_cusp_patterns\":"<<rational<<",\"sextic_plus_cusp_patterns\":"<<seen<<",\"rank\":10,\"rows\":38,\"method\":\"Full Gaussian rank at every leaf after independent invertible row-basis change\"}\n";}
