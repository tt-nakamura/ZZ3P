// uses NTL
//   http://www.shoup.net/ntl

#ifndef __ZZ3P_h__
#define __ZZ3P_h__

#include<NTL/ZZX.h>
#include<NTL/mat_ZZ.h>
#include<NTL/mat_RR.h>

// integers pure cubic field
struct ZZ3P : NTL::vec_ZZ {
	static NTL::mat_ZZ W;// integral basis
	static NTL::ZZ d,d3;// denominator of integral basis
	static NTL::ZZX T;// modulus polynomial of cubic field
	static NTL::Vec<NTL::mat_ZZ> M,N;// multiplication table
	static NTL::Mat<NTL::RR> L;// basis lattice coordinates
    static void init(long m, long cubefree=1);// set m
    static void FundUnit(ZZ3P& u);// Fundamental Unit
    ZZ3P() { SetLength(3); }
};

void set(ZZ3P& a);// a=1

void LatticePoint(NTL::Vec<NTL::RR>& b, const ZZ3P& a);
// b = lattice coordinates of a

void norm(NTL::ZZ& n, const ZZ3P& a);// n = norm(a)
void content(NTL::ZZ& c, const ZZ3P& a);// c = gcd(a0,a1,a2)

ZZ3P& operator*=(ZZ3P& b, const ZZ3P& a);// b *= a

void mul(ZZ3P& c, const ZZ3P& a, const ZZ3P& b);// c = a*b
void sqr(ZZ3P& b, const ZZ3P& a);// b = a*a

void inv(NTL::ZZ& d, ZZ3P& b, const ZZ3P& a);// b/d = 1/a

long divide(ZZ3P& q, const ZZ3P& a, const NTL::ZZ& b);
// if b divides a, set q = b/a and return 1
// else return 0 (q is unchanged)

void PolyFromZZ3P(NTL::ZZ& d, NTL::ZZX& b, const ZZ3P& a);
// b/d = polynomial representation of a in m^{1/3}

#endif // __ZZ3P_h__