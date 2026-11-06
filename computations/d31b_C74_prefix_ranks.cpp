#include <bits/stdc++.h>
using namespace std;using U=uint64_t;
vector<vector<U>> groups; vector<U> initial;uint64_t visited=0,failures=0;int target;
bool insert(vector<U>& b,U v){while(v){int p=63-__builtin_clzll(v);if(b[p])v^=b[p];else{b[p]=v;return true;}}return false;}
void dfs(int i,int left,const vector<U>& B,int rank){if(i==(int)groups.size()){if(left==0){visited++;if(rank!=target){failures++;if(failures<10)cerr<<"BAD rank "<<rank<<"\n";}}return;}vector<U>b=B;int r=rank;for(int n=0;n<=left && n<=(int)groups[i].size();n++){dfs(i+1,left-n,b,r);if(n<left && n<(int)groups[i].size())r+=insert(b,groups[i][n]);}}
int main(int argc,char**argv){ifstream f(argv[1]);int g,cap,k,total;f>>g>>cap>>k>>total;groups.assign(g,{});for(auto&v:groups)for(int j=0;j<cap;j++){U x;f>>x;v.push_back(x);}vector<U>b(64);int rank=0;for(int j=0;j<k;j++){U x;f>>x;rank+=insert(b,x);}target=total+k;dfs(0,total,b,rank);cout<<"D31B_C74_PREFIX_RANK_COMPLETE "<<visited<<" "<<failures<<" "<<target<<"\n";return failures?1:0;}
