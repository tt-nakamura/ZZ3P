# references:
#  H. Wada, "Application of computers to number theory"
#    chapter 4, section 2 (in Japanese)
#  S. A. Hambleton and H. C. Williams
#    "Cubic Fields with Geometry" sections 7.9, 8.8
#  B. N. Delone and D. K. Faddeev, "The Theory of 
#    Irrationalities of the Third Degree" section 59

import numpy as np
from ZZ3P import ZZ3P
from ZZlib import content

def cfrac(B,S,k):
    """ continued fraction expansion in Voronoi algorithm
    B = basis lattice in punctured plane
    S = euclid coefficients (s_i, t_i, s_{i+1}, t_{i+1})
    k = direction of expansion (0 for x, 1 for z)
    """
    q = np.floor(-B[1,k]/B[0,k])
    B[1],B[0] = B[0], B[1] + B[0]*q
    S[1],S[0] = S[0], S[1] + S[0]*int(q)

def VoronoiBasis(e0,e1,e2):
    """ Voronoi basis of lattice spanned by (e0,e1,e2)
    e0: int; e1,e2: ZZ3P object
    return (f1,f2) such that f1 is relative minumum to e0
      and (e0,f1,f2) span the same lattice. 
    ref: Wada, p45, fig4
    """
    A = np.vstack([e1.LatticePoint(),
                   e2.LatticePoint()])/e0
    B = np.c_[A[:,0] - A[:,1], A[:,2]] # puncture
    S = np.eye(2, dtype='object')

    # cfrac in x-direction
    while(B[0,1]*B[1,1] < 0 or
          abs(B[0,1]) < 0.5 or
          abs(B[1,1]) < 0.5): cfrac(B,S,0)

    # make B[0,0] > 0 and B[1,0] > 0
    if B[1,0]<0: B[1],S[1] = -B[1],-S[1]
    else:        B[0],S[0] = -B[0],-S[0]

    # cfrac in z-direction
    while abs(B[0,1]) > 0.5: cfrac(B,S,1)

    # five puncture theorem
    S = [S[0], S[1], S[0] - S[1],
         S[0] + S[1], S[1] + S[0]*2]
    P = np.dot(S,A).astype('float').T
    q = np.floor(P[1])
    P[1] -= q
    r = P[1]**2 + P[2]**2 # cylinder radius
    s = P[1]*2
    q[r<s] += 1
    q[(r>s)&(r>1)] = np.nan # outside of cylinder
    i = np.nanargmin(P[0] - q)

    S = [S[i], S[int(i==0)]]
    e1,e2 = np.dot(S, [e1.v, e2.v])
    e1[0] -= e0*int(q[i])
    return ZZ3P(e1), ZZ3P(e2)

def FundUnit(m):
    """ fundamental unit of pure cubic field Q(m^{1/3})
    ref: Hambleton and Williams, algorithm 7.1
    """
    ZZ3P.init(m)
    d,e0 = 1,1
    u  = ZZ3P([1,0,0])
    e1 = ZZ3P([0,1,0])
    e2 = ZZ3P([0,0,1])

    while True:
        e1,e2 = VoronoiBasis(e0,e1,e2)
        u,d = u*e1, d*e0
        k = content(u.v, d)
        u,d = u//k, d//k
        if d == 1:
            n = u.norm()
            if n == 1: return u
            if n == -1: return -u
        e1,n = e1.inv()
        e2 *= e1
        k = content(e0, n, e2.v)
        e0,e2,n = e0//k, e2//k, n//k
        e0,e1 = n, e1*e0
