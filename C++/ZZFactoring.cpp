// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/vec_ZZ.h>
#include<NTL/pair.h>
using namespace NTL;

void factor(Vec<Pair<long, long> >& f, long n)
// f = prime factorization of |n|
//     vector of (prime, exponent) pair
//     in increasing order of primes
{
	long i(0),j,p,m;
	if(n<0) n=-n;
	if(n==0 || n==1) {
		f.SetLength(0);
		return;
	}
	m = SqrRoot(n);
	PrimeSeq ps;
	while((p = ps.next()) <= m) {
		if(p==0) Error("too large n");
		for(j=0; n%p==0; j++) n/=p;
		if(j==0) continue;
		f.SetLength(i+1);
		f[i].a = p;
		f[i].b = j;
		if(n==1) return;
		i++;
		m = SqrRoot(n);
	}
	f.SetLength(i+1);
	f[i].a = n;
	f[i].b = 1;
}

long CubeFree(long& f, long &g, long m)
// m = fg^2 such that g is largest
// return 0 if m is divisible by cube
// else return 1
// assume m>0
{
    long i,c(1);
    Vec<Pair<long, long> > p;
    factor(p,m);
    f = 1;
    for(i=0; i<p.length(); i++) {
        if(p[i].b&1) f *= p[i].a;
        if(p[i].b>2) c=0;
    }
    g = SqrRoot(m/f);
    return c;
}

long RootInt(long n, long k)
// return n^{1/k}; assume n>0, k>0
{
    long k1(k-1),x,y(long(round(pow(n,1./k))));
    do {
        x = y;
        y = (k1*x + n/power_long(x,k1))/k;
    } while(x>y);
    return x;
}
