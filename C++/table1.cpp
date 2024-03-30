// reference:
//  H. Wada, ``A Table of Fundamental Units of Purely Cubic Fields"
//   Proceedings of the Japan Academy 46 (1970) 1135

#include<fstream>
#include<exception>
#include<cmath>
#include "ZZ3P.h"
using namespace NTL;

main() {
    double l[3], l10(log(10)), L(85);
    std::ofstream fs("table1.txt");
    long m,i;
    ZZ n;
    ZZ3P u;
    ZZX f;
    for(m=251; m<323; m++) {
        try { ZZ3P::init(m); }
        catch(std::exception) { continue; }
        ZZ3P::FundUnit(u);
        PolyFromZZ3P(n,f,u);
        for(i=0; i<3; i++) l[i] = log(f[i])/l10;
        fs << '$' << m << "$ & $";// latex table
        if(!IsOne(n)) fs << '(';
        fs << f[0];
        if(l[0]+l[1] > L) fs << "$ \\\\\n & $ + ";
        else { fs << " + "; l[1] += l[0]; }
        fs << f[1] << "\\theta";
        if(l[1]+l[2] > L) fs << "$ \\\\\n & $ + ";
        else fs << " + ";
        fs << f[2] << "\\theta^2";
        if(!IsOne(n)) fs << ")/" << n;
        fs << "$ \\\\\n";
    }
}
