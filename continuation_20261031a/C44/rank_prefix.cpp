#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
using namespace std;
struct V{uint64_t a=0,b=0;};
V plusv(V x,V y){uint64_t xx=x.a|x.b,yy=y.a|y.b;return {(x.a&~yy)|(y.a&~xx)|(x.b&y.b),(x.b&~yy)|(y.b&~xx)|(x.a&y.a)};}
array<V,64> piv;int n,dim,target;vector<int> ds,mults;vector<vector<V>> blocks;uint64_t leaves=0,bad=0;ofstream out;
int insertv(V v){while(v.a|v.b){int j=63-__builtin_clzll(v.a|v.b);if(!(piv[j].a|piv[j].b)){if(v.b&(1ULL<<j))swap(v.a,v.b);piv[j]=v;return j;}V q=piv[j];if(v.a&(1ULL<<j))swap(q.a,q.b);v=plusv(v,q);}return -1;}
void dfs(int i,int rem,int rank){if(rem==0){++leaves;if(rank!=target){++bad;out<<"[";for(int j=0;j<n;++j){if(j)out<<",";out<<(j<i?mults[j]:0);}out<<"]\n";}return;}if(i==n)return;vector<int> added;for(int m=0;m<=rem/ds[i];++m){if(m){for(int c=(m-1)*ds[i];c<m*ds[i];++c){int j=insertv(blocks[i][c]);if(j>=0)added.push_back(j);}}mults[i]=m;dfs(i+1,rem-m*ds[i],rank+added.size());}for(int j:added)piv[j]={0,0};mults[i]=0;}
int main(int argc,char**argv){if(argc!=3)return 2;ifstream in(argv[1]);in>>n>>dim>>target;if(!in||dim>64)return 3;ds.resize(n);mults.resize(n);blocks.resize(n);for(int i=0;i<n;++i){int cols;in>>ds[i]>>cols;blocks[i].resize(cols);for(auto &v:blocks[i])in>>v.a>>v.b;}if(!in)return 4;out.open(argv[2]);dfs(0,target,0);out.close();cout<<"D31A_C44_PREFIX_RANK_COMPLETE patterns="<<leaves<<" bad="<<bad<<" dimension="<<dim<<"\n";return leaves==4302898?0:5;}
