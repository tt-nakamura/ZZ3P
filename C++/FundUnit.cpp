// references:
//  H. Wada, "Application of computers to number theory"
//    chapter 4, section 2 (in Japanese)
//  S. A. Hambleton and H. C. Williams
//    "Cubic Fields with Geometry" sections 7.9, 8.8
//  B. N. Delone and D. K. Faddeev, "The Theory of 
//    Irrationalities of the Third Degree" section 59

// uses NTL
//   http://www.shoup.net/ntl

#include "ZZ3P.h"
using namespace NTL;

void cfrac(Mat<RR>& B, Mat<ZZ>& S, long k)
// continued fraction expansion in Voronoi algorithm
// B = basis lattice in punctured plane
// S = euclid coefficients (s_i, t_i, s_{i+1}, t_{i+1})
// k = direction of expansion (0 for x, 1 for z)
{
    ZZ q;
    RR t;
    Vec<ZZ> s;
    Vec<RR> b;
    div(t, B[1][k], B[0][k]);
    negate(t,t);
    floor(t,t);
    conv(q,t);
    swap(B[0], B[1]);
    swap(S[0], S[1]);
    mul(b, B[1], t); B[0] += b;
    mul(s, S[1], q); S[0] += s;
}

void VoronoiBasis(ZZ3P& e1, ZZ3P& e2, const ZZ& n)
// Voronoi basis of lattice spanned by (1, e1/n, e2/n)
// (e1,e2) are replaced by new basis
//   such that e1 is relative minumum to 1
// ref: Wada, p45, fig4
{
    static Mat<RR> P(INIT_SIZE,5,3);
    static Mat<RR> A(INIT_SIZE,2,3);
    static Mat<RR> B(INIT_SIZE,2,2);
    static Mat<ZZ> E(INIT_SIZE,2,3);
    static Mat<ZZ> S(INIT_SIZE,2,2);

    long i,j;
    ZZ q;
    RR c,r,s,t;

    LatticePoint(A[0], e1);
    LatticePoint(A[1], e2);
    conv(c,n); inv(c,c); A *= c;
    // puncture
    sub(B[0][0], A[0][0], A[0][1]);
    sub(B[1][0], A[1][0], A[1][1]);
    B[0][1] = A[0][2];
    B[1][1] = A[1][2];
    ident(S,2);
    // cfrac in x-direction
    while(sign(B[0][1]) != sign(B[1][1]) ||
          (B[0][1] < 0.5 && B[0][1] > -0.5) ||
          (B[1][1] < 0.5 && B[1][1] > -0.5))
        cfrac(B,S,0);

    // make B[0][0] > 0 and B[1][0] > 0
    if(sign(B[1][0]) < 0) {
        negate(B[1], B[1]);
        negate(S[1], S[1]);
    }
    else {
        negate(B[0] ,B[0]);
        negate(S[0], S[0]);
    }
    // cfrac in z-direction
    while(B[0][1] > 0.5 || B[0][1] < -0.5)
        cfrac(B,S,1);
    // five puncture theorem
    conv(B,S);
    mul(P[0], B[0], A);
    mul(P[1], B[1], A);
    sub(P[2], P[0], P[1]);
    add(P[3], P[0], P[1]);
    mul(P[4], P[0], 2); P[4] += P[1];
    for(i=0; i<5; i++) {
        floor(t, P[i][1]);
        P[i][1] -= t;
        sqr(r, P[i][1]);
        sqr(s, P[i][2]); r += s;// cylinder radius
        mul(s, P[i][1], 2);
        if(r<s) t++;
        else if(r>=1) continue;// outside of cylinder
        sub(s, P[i][0], t);
        if(i==0 || s<c) { c=s; j=i; conv(q,t); }
    }
    if(j) swap(S[0], S[1]);
    if(j==2) sub(S[0], S[1], S[0]);
    else if(j==3) S[0] += S[1];
    else if(j==4)
        for(i=0; i<2; i++)
            MulAddTo(S[0][i], S[1][i], 2);
    E[0] = e1;
    E[1] = e2;
    mul(e1, S[0], E);
    mul(e2, S[1], E);
    MulSubFrom(e1[0], n, q);
}

void ZZ3P::FundUnit(ZZ3P& u)
// u = fundamental unit of pure cubic field
// assume ZZ3P::init(m) has been executed
// ref: Hambleton and Williams, algorithm 7.1
{
    ZZ d(1),n(1),m,k;
    ZZ3P e1,e2;

    set(u);
    set(e1[1]);
    set(e2[2]);
    for(;;) {
        VoronoiBasis(e1,e2,n);
        d *= n;
        u *= e1;
        if(divide(u,u,d)) {
            norm(k,u);
            if(IsOne(k)) break;
            if(k==-1) { negate(u,u); break; }
            set(d);
        }
        inv(m,e1,e1);// basis1 = (n/m)e1
        e2 *= e1;// basis2 = (1/m)e2
        content(k,e2); GCD(k,k,m); GCD(k,k,n);
        if(!IsOne(k)) {
            n /= k;
            m /= k;
            divide(e2,e2,k);
        }
        e1 *= n;
        n = m;
    }
}