# reference:
#  H. Wada, ``A Table of Fundamental Units of Purely Cubic Fields"
#   Proceedings of the Japan Academy 46 (1970) 1135

from FundUnit import FundUnit

for m in range(251,323):
    try: e = FundUnit(m)
    except: continue
    e,d = e.to_poly()
    if d>1: e = str(e) + '/' + str(d)
    print(m, e)
