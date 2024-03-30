// uses NTL
//   http://www.shoup.net/ntl

#include "ZZ3P.h"
#include<exception>
using namespace NTL;

mat_ZZ ZZ3P::W; // integral basis
ZZ ZZ3P::d; // denominator of integral basis
ZZ ZZ3P::d3; // d^3
ZZX ZZ3P::T; // modulus polynomial (x^3 - m)
Vec<mat_ZZ> ZZ3P::M; // multiplication table of W
Vec<mat_ZZ> ZZ3P::N; // multiplication table of conj(W)
Mat<RR> ZZ3P::L; // w, real(conj(w)), imag(conj(w))

long CubeFree(long&, long&, long);
long RootInt(long, long);

void ZZ3P::init(long m, long cubefree)
// set m of cubic field Q(m^{1/3})
// if cubefree, raise exception if m is not cubefree
// assume m!=0 and m is not cube
// ref: H. Cohen "A Course in Computational
//      Algebraic Number Theory" theorem 6.4.13
{
	long i,j,f,g,fg;

    m = abs(m);
    if(m==0) throw std::runtime_error("m==0");
    i = RootInt(m,3);// i = m^{1/3}
    if(i*i*i == m)
        throw std::runtime_error("m is cube");
    i = CubeFree(f,g,m);// m = fg^2 (g is largest)
    if(cubefree & !i)
        throw std::runtime_error("m is not cubefree");

    SetCoeff(T,3);
    SetCoeff(T,0,-m);// T = x^3 - m

    // integral basis
    // w_i = \sum_{j=0}^2 W[i][j] m^{j/3} (i=0,1,2)
    W.SetDims(3,3); clear(W);

    // multiplication tables
    // w_i * w_j = \sum_{k=0}^2 M[i][j][k] * w_k (j<=i)
    // w_i' * w_i'' = \sum_{k=0}^2 N[i][i][k] * w_k
    // w_i' * w_j'' + w_i'' * w_j'
    //      = \sum_{k=0}^2 N[i][j][k] * w_k (j<i)
    M.SetLength(3);
    N.SetLength(3);

    for(i=0; i<3; i++) {
        M[i].SetDims(i+1, 3);
        N[i].SetDims(i+1, 3);
        clear(M[i]);
        clear(N[i]);
        set(M[i][0][i]);
    }

    fg = f*g;
    i = m%9;
    if(i!=1 && i!=8) {
        // w_0 = 1, w_1 = m^{1/3}, w_2 = m^{2/3}/g
        d = g;
        M[1][1][2] = N[1][1][2] = g;
        M[2][1][0] = fg;
        M[2][2][1] = N[2][2][1] = f;
        N[2][1][0] = -fg;
    }
    else {
        // w_0 = 1, w_1 = m^{1/3},
        // w_2 = (b + am^{1/3} + m^{2/3}/g)/3
        long a(fg%3), b(g%3), ag,bg,abg,afg;
        if(a==2) a = -1;
        if(b==2) b = -1;
        ag = a*g;
        bg = b*g;
        abg = a*bg;
        afg = a*fg;
        d = 3*g;
        W[2][0] = bg;
        W[2][1] = ag;

        M[1][1][0] = -bg;
        M[1][1][1] = -ag;
        M[1][1][2] = d;
        M[2][1][0] = (fg - abg)/3;
        M[2][1][1] = (b-g)/3;
        M[2][1][2] = ag;
        M[2][2][0] = (2*afg - bg - 1)/9;
        M[2][2][1] = (f - ag)/9;
        M[2][2][2] = (2*b + g)/3;

        N[2][0][0] = b;
        N[1][1] = M[1][1];
        N[2][1][0] = -(fg + 2*abg)/3;
        N[2][1][1] = -(b + 2*g)/3;
        N[2][1][2] = 2*ag;
        N[2][2][0] = (2 - afg - bg)/9;
        N[2][2][1] = (f - ag)/9;
        N[2][2][2] = (g-b)/3;
    }
    W[0][0] = W[1][1] = d;
    set(W[2][2]);
    set(N[0][0][0]);
    N[1][0][1] = N[2][0][2] = -1;
    power(d3,d,3);

    // lattice corrdinates
    Mat<RR> K;
    RR t(m),b(3),a;
    inv(a,b); pow(t,t,a);
    SqrRoot(a,b); a/=2; negate(b,a);

    // K[i][0] = real value of m^{i/3}
    // K[i][1] = real part of (m^{1/3}\rho)^i
    // K[i][2] = imaginary part of (m^{1/3}\rho)^i
    //   where \rho = (-1 + \sqrt{-3})/2
    K.SetDims(3,3);
    set(K[0][0]);
    set(K[0][1]);
    sqr(K[2][0], K[1][0] = t);
    div(K[1][1], K[1][0], -2);
    div(K[2][1], K[2][0], -2);
    mul(K[1][2], K[1][0], a);
    mul(K[2][2], K[2][0], b);
    conv(L,W);
    conv(t,d); inv(t,t);
    L *= K;
    L *= t;
    // L[i][0] = real value of w_i
    // L[i][1] = real part of w_i'
    // L[i][2] = imaginary part of w_i'
}

void set(ZZ3P& a) { set(a[0]); clear(a[1]); clear(a[2]); }
// a = 1

void LatticePoint(Vec<RR>& b, const ZZ3P& a)
// b[0] = real value of a
// b[1] = real part of a'
// b[2] = imaginary part of a'
//   where a' is conjugate of a as obtained by
//   replacing m^{1/3} with m^{1/3}\rho
{
    conv(b,a);
    b *= ZZ3P::L;
}

void norm(NTL::ZZ& n, const ZZ3P& a)
// n = norm of (a) = a * a' ** a''
//     where a'' is complex congugate of a'
{
    ZZX p;
    mul(p.rep, a, ZZ3P::W);
    p.normalize();
    NormMod(n, p, ZZ3P::T);
    n /= ZZ3P::d3;
}

void content(ZZ& c, const ZZ3P& a)
// c = gcd(a[0], a[1], a[2])
{
    ZZX b; conv(b,a);
    content(c,b);
    abs(c,c);
}

ZZ3P& operator*=(ZZ3P& b, const ZZ3P& a) { mul(b,b,a); return b; }
// b = a*b

void mul(ZZ3P& c, const ZZ3P& a, const ZZ3P& b)
// c = a*b; &c==&a or &c==&b is allowed
{
    if(&c==&a) { mul(c, ZZ3P(a), b); return; }
    if(&c==&b) { mul(c, a, ZZ3P(b)); return; }
    long i,j,k;
    ZZ t;
    clear(c);
    for(i=0; i<3; i++) {
        for(j=0; j<=i; j++) {
            mul(t, a[i], b[j]);
            if(j<i) MulAddTo(t, a[j], b[i]);
            for(k=0; k<3; k++)
                MulAddTo(c[k], t, ZZ3P::M[i][j][k]);
        }
    }
}

void sqr(ZZ3P& b, const ZZ3P& a)
// b = a*a; &b==&a is allowed
{
    if(&b==&a) { sqr(b, ZZ3P(a)); return; }
    long i,j,k;
    ZZ t;
    clear(b);
    for(i=0; i<3; i++) {
        for(j=0; j<=i; j++) {
            if(i==j) sqr(t, a[i]);
            else { mul(t, a[i], a[j]); t <<= 1; }
            for(k=0; k<3; k++)
                MulAddTo(b[k], t, ZZ3P::M[i][j][k]);
        }
    }
}

void inv(ZZ& d, ZZ3P& b, const ZZ3P& a)
// b/d = 1/a = a'*a''/norm(a); &b==&a is allowed
// d>0, d and content(b) are coprime
{
    if(&b==&a) { inv(d, b, ZZ3P(a)); return; }
    long i,j,k;
    ZZ t;
    clear(b);
    for(i=0; i<3; i++) {
        for(j=0; j<=i; j++) {
            mul(t, a[i], a[j]);
            for(k=0; k<3; k++)
                MulAddTo(b[k], t, ZZ3P::N[i][j][k]);
        }
    }
    clear(d);
    for(i=0; i<3; i++) {
        for(j=0; j<=i; j++) {
            mul(t, a[i], b[j]);
            if(j<i) MulAddTo(t, a[j], b[i]);
            MulAddTo(d, t, ZZ3P::M[i][j][0]);
        }
    }
    if(sign(d) < 0) {
        negate(d,d);
        negate(b,b);
    }
    GCD(t,d,b[0]); if(IsOne(t)) return;
    GCD(t,t,b[1]); if(IsOne(t)) return;
    GCD(t,t,b[2]); if(IsOne(t)) return;
    for(i=0; i<3; i++) b[i] /= t;
    d /= t;// minimize denominator
}

long divide(ZZ3P& q, const ZZ3P& a, const NTL::ZZ& b)
// if b divides a, set q = b/a and return 1
// else return 0 (q is unchanged)
{
    ZZX t; conv(t,a);
    if(!divide(t,t,b)) return 0;
    VectorCopy(q,t,3);
    return 1;
}

void PolyFromZZ3P(ZZ& d, ZZX& b, const ZZ3P& a)
// b/d = polynomial representation of a in m^{1/3}
// d and content(b) are coprime
{
    ZZ c;
    mul(b.rep, a, ZZ3P::W);
    b.normalize();
    content(c,b); abs(c,c);
    d = ZZ3P::d;
    if(!IsOne(c)) GCD(c,c,d);
    if(!IsOne(c)) { b/=c; d/=c; }
}