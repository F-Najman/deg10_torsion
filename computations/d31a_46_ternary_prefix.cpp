#include <algorithm>
#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdint>
#include <stdexcept>
struct V{uint64_t one=0,two=0;};constexpr uint64_t mask=(1ull<<45)-1;
V add(V a,V b){uint64_t a0=mask&~(a.one|a.two),b0=mask&~(b.one|b.two);return{(a0&b.one)|(a.one&b0)|(a.two&b.two),(a0&b.two)|(a.two&b0)|(a.one&b.one)};}
int rank10(std::array<V,10>const&vv){std::array<V,45>basis{};int r=0;for(V v:vv){while(v.one|v.two){int p=63-__builtin_clzll(v.one|v.two);uint64_t bit=1ull<<p;if(!(basis[p].one|basis[p].two)){if(v.two&bit)std::swap(v.one,v.two);basis[p]=v;++r;break;}V b=basis[p];if(v.one&bit)std::swap(b.one,b.two);v=add(v,b);}}return r;}
std::vector<std::vector<V>> flags;std::array<V,10>selected;std::array<int,22>mult;int top;uint64_t count;
void visit(int at,int left,int pos){if(at==22){if(left)return;if(pos!=10)throw std::runtime_error("wrong width");if(rank10(selected)!=10){std::cout<<"DEPENDENCY";for(auto m:mult)std::cout<<" "<<m;std::cout<<std::endl;throw std::runtime_error("rank below10");}++count;return;}int cap=std::min({top,5,left});if(left>cap+5*(21-at))return;for(int m=0;m<=cap;++m){mult[at]=m;if(m)selected[pos+m-1]=flags[at][m-1];visit(at+1,left-m,pos+m);}}
int main(int argc,char**argv){if(argc!=2)return 2;for(int a=0;a<3;++a)for(int b=0;b<3;++b){V x{a==1?1ull:0ull,a==2?1ull:0ull},y{b==1?1ull:0ull,b==2?1ull:0ull};auto z=add(x,y);int n=(z.one&1)?1:(z.two&1)?2:0;if(n!=(a+b)%3)return 3;}
std::ifstream in(argv[1]);int baseA,baseB;in>>baseA>>baseB;std::vector<std::vector<V>> raw(22);for(int i=0;i<22;++i){int n;in>>n;for(int j=0;j<n;++j){V v;in>>v.one>>v.two;raw[i].push_back(v);}}if(!in)return 4;
for(int base:{baseA,baseB}){flags.clear();flags.push_back(raw[base]);for(int i=0;i<22;++i)if(i!=base)flags.push_back(raw[i]);for(top=1;top<=10;++top){for(int j=0;j<top;++j)selected[j]=flags[0][j];mult[0]=top;count=0;visit(1,10-top,top);std::cout<<base<<" "<<top<<" "<<count<<std::endl;}}std::cout<<"C46_TERNARY_ALL_PREFIX_RANKS_COMPLETE"<<std::endl;return 0;}
