from sympy import bernoulli
import utils as ut


#Dedekind-Rademacher sum
def testdederad(v):
    final = 0
    def b1(x):
        return bernoulli(1,x)
    for k in [0,1,2]:
        final += ut.Dedekind(v,k,[b1]*3)
    derived = (v[0]**2+v[1]**2+v[2]**2)/(12*v[0]*v[1]*v[2]) - 1/4
    print(v)
    #print('realval = ', final)
    #print('derived value =', derived)
    print('ratio =', final/derived)
    return 0

N=5
for n in range(N):
    v = ut.randvs(2)
    testdederad(v)