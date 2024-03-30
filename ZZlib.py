import numpy as np
from fractions import gcd

def SqrRoot(n):
    """ floor(sqrt(n)) for huge integer n """
    x,y = n+1,n
    while x>y: x,y = y, (y + n//y)//2
    return x

def RootInt(n,k):
    """ floor((n)++(1/k)) for huge integer n """
    x,y,k1 = n+1,n,k-1
    while x>y: x,y = y, (k1*y + n//y**k1)//k
    return x

def content(*a):
    """ greatest common divisor of a's """
    a = np.r_[a]
    c = a[0]
    for a in a[1:]: c = gcd(c,a)
    return abs(c)
