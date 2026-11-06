#include <array>
#include <fstream>
#include <vector>
#include <iostream>
#include <cstdint>
#include <stdexcept>
struct V {std::array<uint64_t,4>x{};};
std::array<V,10> selected;std::vector<std::vector<V>> rats;uint64_t total;
int highest(V const&v){for(int j=3;j>=0;--j)if(v.x[j])return 64*j+63-__builtin_clzll(v.x[j]);return -1;}
int rank10(){std::array<V,256>b{};int r=0;for(V v:selected){int p;while((p=highest(v))>=0){if(highest(b[p])<0){b[p]=v;++r;break;}for(int j=0;j<4;++j)v.x[j]^=b[p].x[j];}}return r;}
void visit(int i,int rem,int pos){if(i==(int)rats.size()-1){for(int j=0;j<rem;++j)selected[pos+j]=rats[i][j];if(pos+rem!=10||rank10()!=10)throw std::runtime_error("rank below ten");++total;return;}for(int m=0;m<=rem;++m){if(m)selected[pos+m-1]=rats[i][m-1];visit(i+1,rem-m,pos+m);}}
V readv(std::ifstream&in){V v;for(auto&x:v.x)in>>x;return v;}
int main(int argc,char**argv){if(argc!=2)return 2;std::ifstream in(argv[1]);int nf;in>>nf;for(int f=0;f<nf;++f){int fixed,rem,nrat;in>>fixed>>rem>>nrat;if(fixed+rem!=10)return 3;for(int j=0;j<fixed;++j)selected[j]=readv(in);rats.clear();for(int j=0;j<nrat;++j){std::vector<V>a;for(int k=0;k<rem;++k)a.push_back(readv(in));rats.push_back(a);}if(!in)return 4;total=0;visit(0,rem,fixed);std::cout<<f<<" "<<total<<"\n";}std::cout<<"PARENT_MIXED_FULL_RANK_COMPLETE "<<nf<<"\n";return 0;}
