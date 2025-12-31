
from sympy import bernoulli
import utils as ut
import zetas as zet 
import dimtwo as dt
import math


def lemma51(k,v,q):
    normq = sum([q[j] for j in [0,1,2] if j!= k])
    derived = dt.Zvalued(v,q,k)*ut.zeta(normq)
    final = zet.Zperp(k,v,q)
    print('realval = ', final)
    print('derived value =', derived)
    if derived == 0:
        print('DIVISON BY 0',  final,'----',derived)
    else:
        print('ratio =', final/derived) 
    return 0


def lemma52(v,q,k,l):
    U = [[1,-1,-1],[1,-1,1],[1,1,-1]]
    ul = U[l]
    final = zet.Yperpu(v, k, q, ul)
    normq = sum([q[j] for j in [0,1,2] if j!= k])
    Qk = dt.Qvaluekl(v,q,k,l)
    derived = Qk*ut.zeta(normq)
    print('realval = ', final)
    print('derived value =', derived)
    if derived == 0:
        print('DIVISON BY 0',  final,'----',derived)
    else:
        print('ratio =', final/derived) 
    return 0

def coro53(v,q):
    fs = [lambda x, s=qj: bernoulli(s, x) for qj in q]
    final = 0
    for k in range(len(v)):
        if q[k] == 1:
            final += ut.Dedekind(v,k,fs)
    derived = 0
    for k in range(len(v)):
        if q[k] == 1:
            derived += v[k]*(dt.Zvalued(v,q,k)/2 + dt.Qvalued(v,q,k))
    derived *= dt.constquno(q)
    print(derived)
    delta = dt.delta2d(q)
    derived += -delta
    print('realval = ', final)
    print('derived value =', derived)
    if derived == 0:
        print('DIVISON BY 0',  final,'----',derived)
    else:
        print('ratio =', final/derived) 
    return 0


def remark54(v,q):
    Q = [1,q,q]
    func = [lambda x, s=qj : bernoulli(s,x) for qj in Q]
    final = ut.Dedekind(v,0,func)
    derived = (-1)**q*(v[1]*v[2])**(-q) + dt.Qvalued(v,Q,0)
    derived *= -bernoulli(2*q, 0) * v[0]*(math.factorial(q)**2) / math.factorial(2*q)
    derived += -bernoulli(q,0)**2 
    print('realval = ', final)
    print('derived value =', derived)
    if derived == 0:
        print('DIVISON BY 0',  final,'----',derived)
    else:
        print('ratio =', final/derived) 
    return 0

v=[2,3,7]
k=0
l = 1
q = [1,1,1]

lemma51(k,v,q)
lemma52(v,q,k,l)
coro53(v,q)
remark54(v,4)