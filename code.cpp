//aravind ravishankar
#include<bits/stdc++.h>
#include<NTL/ZZ.h>
#include<NTL/ZZX.h>
#include<NTL/ZZ_pX.h>
#include<NTL/ZZXFactoring.h> 
#include<NTL/ZZ_pEX.h> 

using namespace NTL;
using namespace std;

typedef pair<int,int>   II;
typedef vector< II >      VII;
typedef vector<int>     VI;
typedef vector< VI > 	VVI;
typedef long long int 	LL;

#define PB push_back
#define MP make_pair
#define F first
#define S second
#define SZ(a) (int)(a.size())
#define ALL(a) a.begin(),a.end()
#define SET(a,b) memset(a,b,sizeof(a))

#define si(n) scanf("%d",&n)
#define dout(n) printf("%d\n",n)
#define sll(n) scanf("%lld",&n)
#define lldout(n) printf("%lld\n",n)
#define fast_io ios_base::sync_with_stdio(false);cin.tie(NULL)

#define TRACE

#ifdef TRACE
#define trace(...) __f(#__VA_ARGS__, __VA_ARGS__)
template <typename Arg1>
void __f(const char* name, Arg1&& arg1){
	cerr << name << " : " << arg1 << std::endl;
}
template <typename Arg1, typename... Args>
void __f(const char* names, Arg1&& arg1, Args&&... args){
	const char* comma = strchr(names + 1, ',');cerr.write(names, comma - names) << " : " << arg1<<" | ";__f(comma+1, args...);
}
#else
#define trace(...)
#endif

//FILE *fin = freopen("in","r",stdin);
//FILE *fout = freopen("out","w",stdout);

bool check_if_power(ZZ n){
  int logn = (log(n)/log(ZZ(2)))+1;
  for(int i=2;i<=logn;i++){
    ZZ l=ZZ(2),r=ZZ(n);
    while(l<=r){
      ZZ mid = (l+r)/ZZ(2);
      ZZ tmp = power(mid,i);
      if(tmp>n)
        r=mid-1;
      else if(tmp<n)
        l=mid+1;
      else return true;
    }
  }
  return false;
}
int findsmallestr(ZZ n){
  double logn = (log(n)/log(ZZ(2)));
  int log2n=(int)(logn*logn);
  ZZ ret(2);
  int rett=2;
  while(1){
    ZZ tmp(1);
    int fl=0;
    for(int i=1;i<=log2n;i++){
      tmp=MulMod(tmp,n,ret);
      if(tmp==ZZ(1)){
        fl=1;
        break;
      }
    }
    if(fl==0)
      break;
    ret++;
    rett++;
  }
  return rett;
}
int phi(ZZ n){
  int ret=0;
  for(ZZ i(0);i<=n;i++){
    ret+=(GCD(i,n)==ZZ(1));
  }
  return ret;
}
int main()
{
  ZZ n;
  cin>>n;
  if(check_if_power(n)){
    cout<<"Not Prime\n";
    return 0;
  }
  int r = findsmallestr(n);
  for(int i=1;i<=r;i++){
    if(GCD(n,ZZ(r))>1&&GCD(n,ZZ(r))<n){
      cout<<"Not Prime\n";
      return 0;
    }
  }
  if(n<=r){
    cout<<"Prime\n";
    return 0;
  }
  int phir=phi(ZZ(r));
  double tmp=sqrt(phir);
  tmp=tmp*log(n);
  int lim=floor(tmp);
  ZZ_pX lhs,rhs,mod;
  ZZ_pXModulus Mod;
  ZZ_p::init(n); 
  mod=ZZ_pX(r, to_ZZ_p(1)); 
  sub(mod,mod,1); 
  build(Mod,mod); 
  for(int a=1;a<=lim;a++){
    PowerXPlusAMod(lhs,to_ZZ_p(a),n,Mod); 
    PowerXMod(rhs,n,Mod); 
    if(!(IsZero((lhs - rhs - a) % Mod))) {
      cout<<"Not Prime\n";
      return 0;
    }
  }
  cout<<"Prime\n";
	return 0;
}
