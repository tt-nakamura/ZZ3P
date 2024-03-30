import numpy as np
from sympy.ntheory.factor_ import core
from ZZlib import SqrRoot,RootInt,content
from fractions import gcd

class ZZ3P:
    """ integers in pure cubic field Q(m^{1/3}) """
    @staticmethod
    def init(m, CubeFree=True):
        """ set m if Q(m^{1/3})
        if cubeFree, raise exception if m is not cubefree
        raise exception if m is cube
        assume m>0
        ref: H. Cohen "A Course in Computational
             Algebraic Number Theory" theorem 6.4.13
        """
        if RootInt(m,3)**3 == m:
            raise RuntimeError('m is cube')
        # m = fg^2 where g is largest
        f = core(m)
        g = SqrRoot(m//f)
        if CubeFree and (gcd(f,g)>1 or core(g)!=g):
            raise RuntimeError('m is not cubefree')

        # integral basis
        # w_i = \sum_{j=0}^2 W[i,j] m^{j/3} (i=0,1,2)
        W = np.zeros((3,3), dtype='object')
        # multiplication tables
        # w_i * w_j = \sum_{k=0}^2 M[i,j,k] * w_k
        # w_i' * w_i'' = \sum_{k=0}^2 N[i,i,k] * w_k
        # w_i' * w_j'' + w_i'' * w_j'
        #      = \sum_{k=0}^2 N[i,j,k] * w_k (j<i)
        M = np.zeros((3,3,3), dtype='object')
        N = np.zeros((3,3,3), dtype='object')

        M[0] = np.eye(3, dtype='object')
        M[1:,0] = M[0,1:]

        fg = f*g
        k = m%9
        if k!=1 and k!=8:
            # w_0 = 1, w_1 = m^{1/3}, w_2 = m^{2/3}/g
            d = g

            M[1,1,2] = g
            M[1,2,0] = fg
            M[2,1,0] = fg
            M[2,2,1] = f

            N[1,1,2] = g
            N[2,1,0] = -fg
            N[2,2,1] = f

        else:
            # w_0 = 1, w_1 = m^{1/3},
            # w_2 = (b + am^{1/3} + m^{2/3}/g)/3
            a,b = fg%3,g%3
            if a==2: a=-1
            if b==2: b=-1
            d,ag,bg = 3*g,a*g,b*g
            abg,afg = a*bg,a*fg
            W[2,0] = bg
            W[2,1] = ag         

            M[1,1,:2] = -W[2,:2]
            M[1,1,2] = d
            M[1,2,0] = (fg - abg)//3
            M[1,2,1] = (b-g)//3
            M[1,2,2] = ag
            M[2,1] = M[1,2]
            M[2,2,0] = (2*afg - bg - 1)//9
            M[2,2,1] = (f - ag)//9
            M[2,2,2] = b - M[1,2,1]

            N[2,0,0] = b
            N[1,1] = M[1,1]
            N[2,1,0] = 2*M[1,2,0] - fg
            N[2,1,1] = -M[1,2,1] - g
            N[2,1,2] = 2*ag
            N[2,2,0] = (2 - afg - bg)//9
            N[2,2,1] = M[2,2,1]
            N[2,2,2] = -M[1,2,1]
           
        W[0,0] = d
        W[1,1] = d
        W[2,2] = 1
        N[0,0,0] = 1
        N[1,0,1] = -1
        N[2,0,2] = -1

        t = m**(1/3)
        t2 = t**2
        s3 = np.sqrt(3)/2
        L = [[1, 1, 0],
             [t, -t/2, t*s3],
             [t2, -t2/2, -t2*s3]]
        L = np.dot(W,L)/d;
        # L[i,0] = real value of w_i
        # L[i,1] = real part of w_i'
        # L[i,2] = imaginary part of w_i'

        ZZ3P.d = d
        ZZ3P.W = W
        ZZ3P.M = M
        ZZ3P.N = N
        ZZ3P.L = L
                  
    def __init__(a,v):
        a.v = np.asarray(v, dtype='object')

    def __repr__(a):
        return str(a.v)

    def __getitem__(a,i):
        return a.v[i]

    def __neg__(a):
        return ZZ3P(-a.v)

    def __mul__(a,b): # b is either int or ZZ3P
        if not isinstance(b, ZZ3P): v = a.v*b
        else: v = np.dot(a.v, np.dot(b.v, a.M))
        return ZZ3P(v)

    def __floordiv__(a,b):
        return ZZ3P(a.v//b)

    def inv(a):
        """ return b,d such that b/d = 1/a = a'a''/N(a)
        d>0, d and content(b) are coprime
        """
        v = np.dot(a.v, np.dot(a.v, a.N))
        n = np.dot(a.v, np.dot(v, a.M[:,:,0]))
        if n<0: v,n = -v,-n
        c = content(v,n)
        return ZZ3P(v//c), n//c
        
    def norm(a): # a * a' * a''
        v = np.dot(a.v, np.dot(a.v, a.N))
        n = np.dot(a.v, np.dot(v, a.M[:,:,0]))
        return n

    def LatticePoint(a):
        """
        a[0] = real value of a
        a[1] = real part of a'
        a[2] = imaginary part of a'
        """ 
        return np.dot(a.v, a.L).astype('float')

    def content(a): # gcd(a[0],a[1],a[2])
        return content(a.v)

    def to_poly(a):
        """
        return b,d such that
        b/d = polynomial representation of a in m^{1/3}
        d>0, d and content(b) are coprime
        """
        v = np.dot(a.v, a.W)
        m = content(v, a.d)
        return ZZ3P(v//m), a.d//m
