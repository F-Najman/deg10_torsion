#include <array>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
struct Col { uint64_t a=0,b=0; };
Col add(Col x,Col y){return {x.a^y.a,x.b^y.b};}
Col mul(Col x,int c){if(c==0)return {};if(c==1)return x;if(c==2)return{x.b,x.a^x.b};return{x.a^x.b,x.a};}
struct Basis {array<Col,36> v{};int rank=0;void push(Col x){while(x.a||x.b){int p=__builtin_ctzll(x.a|x.b);int c=((x.a>>p)&1)+2*((x.b>>p)&1);if(v[p].a||v[p].b)x=add(x,mul(v[p],c));else{v[p]=mul(x,c==2?3:(c==3?2:1));rank++;break;}}}};
array<array<Col,5>,40> jets;array<int,40> partner{},cnt{};uint64_t total=0,bad=0;array<uint64_t,3> dist{};
void rec(int start,int left,Basis B,int pairs){if(left==0){total++;if(B.rank!=5-pairs){bad++;if(bad<=20){cerr<<"FAIL rank "<<B.rank<<" pairs "<<pairs<<" divisor ";for(int i=0;i<40;i++)if(cnt[i])cerr<<i<<":"<<cnt[i]<<" ";cerr<<"\n";}}else dist[pairs]++;return;}
 for(int i=start;i<40;i++){Basis T=B;for(int m=1;m<=left;m++){T.push(jets[i][m-1]);cnt[i]=m;int extra=partner[i]<i?min(m,cnt[partner[i]]):0;rec(i+1,left-m,T,pairs+extra);}cnt[i]=0;}}
int main(int argc,char**argv){if(argc!=3)return 2;ifstream in(argv[1]);int n,r;in>>n>>r;if(n!=40||r!=36)return 3;for(int i=0;i<40;i++){in>>partner[i];for(int j=0;j<36;j++)for(int m=0;m<5;m++){int c;in>>c;if(c<0||c>3)return 4;if(c&1)jets[i][m].a|=uint64_t(1)<<j;if(c&2)jets[i][m].b|=uint64_t(1)<<j;}}
 for(int i=0;i<40;i++)if(partner[i]<0||partner[i]>=40||partner[partner[i]]!=i||partner[i]==i)return 5;rec(0,5,Basis{},0);if(total!=1086008)return 6;ofstream out(argv[2]);out<<"{\"status\":\""<<(bad?"FAIL":"PASS")<<"\",\"divisors\":"<<total<<",\"bad\":"<<bad<<",\"rank_distribution_by_pulledback_degree\":["<<dist[0]<<","<<dist[1]<<","<<dist[2]<<"]}\n";cout<<"D31B_M333_PREFIX_GF4_COMPLETE "<<total<<" bad="<<bad<<" ranks="<<dist[0]<<","<<dist[1]<<","<<dist[2]<<"\n";return bad?1:0;}
