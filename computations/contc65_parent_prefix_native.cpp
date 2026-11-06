#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdint>
#include <stdexcept>
struct V {uint64_t lo,hi;};
struct G {int cost,cap;std::vector<V> c;};
std::vector<G> groups;std::array<V,10> selected;uint64_t count;int top;
int rank10(){std::array<V,128>b{};int rank=0;for(V v:selected){while(v.lo||v.hi){int p=v.hi?127-__builtin_clzll(v.hi):63-__builtin_clzll(v.lo);if(!b[p].lo&&!b[p].hi){b[p]=v;++rank;break;}v.lo^=b[p].lo;v.hi^=b[p].hi;}}return rank;}
int capacity(int at){int s=0;for(int j=at;j<(int)groups.size();++j)s+=groups[j].cost*(groups[j].cost==1?std::min(top,5):2);return s;}
void visit(int at,int left,int pos){if(at==(int)groups.size()){if(left)return;if(pos!=10||rank10()!=10)throw std::runtime_error("dependent restriction prefix");++count;return;}const auto&g=groups[at];int cap=g.cost==1?std::min(top,5):2;int lim=std::min(cap,left/g.cost);for(int m=0;m<=lim;++m){if(m&&pos+m*g.cost<=10)for(int k=(m-1)*g.cost;k<m*g.cost;++k)selected[pos+k]=g.c[k];int rem=left-m*g.cost;if(rem<=capacity(at+1))visit(at+1,rem,pos+m*g.cost);}}
int main(int argc,char**argv){if(argc!=2)return 2;std::ifstream in(argv[1]);int ng,nw;in>>ng>>nw;if(ng!=30||nw!=106)return 3;for(int i=0;i<ng;++i){G g;int nv;in>>g.cost>>g.cap>>nv;for(int j=0;j<nv;++j){V v;in>>v.lo>>v.hi;g.c.push_back(v);}groups.push_back(g);}if(!in)return 4;uint64_t total=0;for(top=1;top<=10;++top){for(int j=0;j<top;++j)selected[j]=groups[0].c[j];count=0;visit(1,10-top,top);total+=count;std::cout<<top<<" "<<count<<"\n";}std::cout<<"C65_PARENT_FULL_PREFIX_NATIVE_COMPLETE "<<total<<"\n";return 0;}
