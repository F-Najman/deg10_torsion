#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using Mat=std::vector<std::vector<int>>;
int add9[9][9],sub9[9][9],mul9[9][9],inv9[9];
void require(bool ok){if(!ok)throw std::runtime_error("independent finite-field check failed");}
int mod3(int a){a%=3;return a<0?a+3:a;}
int reduce(Mat& m,std::vector<int>* piv=nullptr){int row=0,n=m.size(),w=m[0].size();for(int j=0;j<w&&row<n;j++){int a=row;while(a<n&&!m[a][j])a++;if(a==n)continue;std::swap(m[a],m[row]);int inv=inv9[m[row][j]];for(int k=0;k<w;k++)m[row][k]=mul9[m[row][k]][inv];for(int i=0;i<n;i++)if(i!=row&&m[i][j]){int b=m[i][j];for(int k=0;k<w;k++)m[i][k]=sub9[m[i][k]][mul9[b][m[row][k]]];}if(piv)piv->push_back(j);row++;}return row;}
uint32_t u32(std::ifstream& in){uint32_t v;in.read((char*)&v,4);require(bool(in));return v;}
uint16_t u16(std::ifstream& in){uint16_t v;in.read((char*)&v,2);require(bool(in));return v;}
int byte(std::ifstream& in){unsigned char b;in.read((char*)&b,1);require(bool(in)&&b<9);return b;}
int main(int argc,char**argv){try{require(argc==2);std::ifstream in(argv[1],std::ios::binary);int c0=u32(in),c1=u32(in),count=u32(in);require(count==5512&&c0<3&&c1<3);for(int a=0;a<9;a++)for(int b=0;b<9;b++){int a0=a%3,a1=a/3,b0=b%3,b1=b/3;add9[a][b]=mod3(a0+b0)+3*mod3(a1+b1);sub9[a][b]=mod3(a0-b0)+3*mod3(a1-b1);mul9[a][b]=mod3(a0*b0-c0*a1*b1)+3*mod3(a0*b1+a1*b0-c1*a1*b1);}for(int a=1;a<9;a++){inv9[a]=0;for(int b=1;b<9;b++)if(mul9[a][b]==1)inv9[a]=b;require(inv9[a]!=0);}
int rat[24][4][25],quad[108][25][2],cubic[320][25][3];for(auto& p:rat)for(auto& j:p)for(int& a:j)a=byte(in);for(auto& p:quad)for(auto& row:p)for(int& a:row)a=byte(in);for(auto& p:cubic)for(auto& row:p)for(int& a:row)a=byte(in);
uint64_t triples=0;
for(int rec=0;rec<count;rec++){int type=u16(in),a=u16(in),b=u16(in),c=u16(in);Mat mt(22,std::vector<int>(24));for(auto& row:mt)for(int& z:row)z=byte(in);std::array<int,24> multiplicity{};Mat bt(3,std::vector<int>(25));if(type==1){require(a<24&&b<24&&c<24);multiplicity[a]++;multiplicity[b]++;multiplicity[c]++;int col=0;for(int p=0;p<24;p++)for(int j=0;j<multiplicity[p];j++)for(int k=0;k<25;k++)bt[col+j][k]=rat[p][j][k];for(int p=0,cc=0;p<24;p++){for(int j=0;j<multiplicity[p];j++){for(int k=0;k<25;k++)bt[cc][k]=rat[p][j][k];cc++;}}}else if(type==2){require(a<108&&b<24);multiplicity[b]=1;for(int k=0;k<25;k++){bt[0][k]=quad[a][k][0];bt[1][k]=quad[a][k][1];bt[2][k]=rat[b][0][k];}}else{require(type==3&&a<320);for(int k=0;k<25;k++)for(int j=0;j<3;j++)bt[j][k]=cubic[a][k][j];}
std::vector<int> piv;require(reduce(bt,&piv)==3);Mat expected;for(int f=0;f<25;f++){bool pivot=false;for(int p:piv)if(p==f)pivot=true;if(pivot)continue;std::vector<int> row(24);for(int j=0;j<24;j++){int z=rat[j][multiplicity[j]][f];for(int p=0;p<3;p++)z=sub9[z][mul9[bt[p][f]][rat[j][multiplicity[j]][piv[p]]]];row[j]=z;}expected.push_back(row);}require(expected.size()==22);Mat mcopy=mt,ecopy=expected,stack=mt;stack.insert(stack.end(),expected.begin(),expected.end());int rank=reduce(mcopy);require(rank==reduce(ecopy)&&rank==reduce(stack));
for(int i=0;i<22;i++){std::array<int,22> u{};int pi=-1;for(int t=0;t<22;t++)if(mt[t][i]){pi=t;break;}require(pi>=0);for(int t=0;t<22;t++)u[t]=mul9[mt[t][i]][inv9[mt[pi][i]]];for(int j=i+1;j<23;j++){std::array<int,22> v{};int pj=-1;for(int t=0;t<22;t++){v[t]=sub9[mt[t][j]][mul9[mt[pi][j]][u[t]]];if(pj<0&&v[t])pj=t;}require(pj>=0);int iv=inv9[v[pj]];for(int t=0;t<22;t++)v[t]=mul9[v[t]][iv];for(int k=j+1;k<24;k++){int fu=mt[pi][k],fv=sub9[mt[pj][k]][mul9[fu][u[pj]]];bool nonzero=false;for(int t=0;t<22;t++)if(sub9[sub9[mt[t][k]][mul9[fu][u[t]]]][mul9[fv][v[t]]]){nonzero=true;break;}require(nonzero);triples++;}}}
}
require(triples==11156288);require(in.peek()==std::char_traits<char>::eof());std::cout<<"GONALITY420_INDEPENDENT_REPLAY_PASS complements="<<count<<" triples="<<triples<<" kernels_and_projection_spaces_verified=true\n";return 0;}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
