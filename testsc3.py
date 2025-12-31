import scipy.integrate as integrate
import matplotlib.pyplot as plt
import itertools as iter
from sympy import bernoulli
import utils as ut
import zetas as zet 
import cmath
import math


def testexamp321(v,a0,a1,a2):
    a = [a0,a1,a2]
    func = [lambda x , s = aj : ut.frac(x)-s for aj in a]
    final = 0
    for k in range(0,3):
        final += ut.Dedekind(v,k,func)
    derived = -1 + (a0+a1+a2) - a1*a2 - a1*a0 - a2*a0
    derived += v[0]*(1/2-a1)*(1/2-a2) + v[1]*(1/2-a0)*(1/2-a2) + v[2]*(1/2-a1)*(1/2-a0)
    derived += (v[0]**2 + v[1]**2+v[2]**2)/(12*v[0]*v[1]*v[2])
    print('realval = ', final)
    print('derived value =', derived)
    print('ratio =', final/derived)   
    return 0

def testexamp322(v,a):
    def fa(x):
        return ut.frac(x)-a
    final = 0
    functions = [fa,fa,fa]
    for k in range(0,3):
        final += ut.Dedekind(v,k,functions)
    derived = (v[0]**2 + v[1]**2+v[2]**2)/(12*v[0]*v[1]*v[2]) + (1/2 - a)*(1/2-a)*(v[0]+v[1]+v[2]) -1 + 3*a-3*a**2
    print('realval = ', final)
    print('derived value =', derived)
    print('ratio =', final/derived)  
    return 0


 
def testexamp351(v,k):
    def fourexp(x):
        return cmath.exp(2*math.pi*1j*x)
    r = len(v)-1
    functions = [fourexp]*(r+1)
    final = ut.Dedekind(v,k,functions)
    derived = -1
    sumvk = sum(v) - v[k]
    if sumvk%v[k] == 0:
        derived += v[k]
    print('realval = ', final)
    print('derived value =', derived)
    print('ratio =', final/derived)  
    return 0


def testexamp352(v,k):
    def normcos(x):
        return math.cos(2*math.pi*x)
    r = len(v)-1
    functions = [normcos]*(r+1)
    final = ut.Dedekind(v,k,functions)
    derived = -1
    betak = len(ut.U(v,k))
    if betak != 0:
        derived += betak*v[k]*(2**(-r))
    print('realval = ', final)
    print('derived value =', derived)
    print('ratio =', final/derived)  
    return 0

def testexamp353(v,k):
    def normsin(x):
        return math.sin(2*math.pi*x)
    r = len(v)-1
    dif = (1+(-1)**r)*(2*1j)**(-r)
    functions = [normsin]*(r+1)
    final = ut.Dedekind(v,k,functions)
    Uvk = ut.U(v,k)
    barU = []
    for u in Uvk:
        if u[0] == 1:
            barU.append(u)
    barbetak = 0
    for u in barU:
        barbetak += ut.prodlist(u)
    derived = dif*barbetak*v[k]
    print('realval = ', final)
    print('derived value =', derived)
    if derived == 0:
        print('DIVISON BY 0',  final,'----',derived)
    else:
        print('ratio =', final/derived)  
    return 0


def testeq11(v,q):
    funcq = [lambda x , s = qj : ut.frac(x)**s for qj in q]
    recq = 0 
    for k in range(0,len(v)):
        recq += ut.Dedekind(v, k, funcq)
    SS = ut.S(q)
    final = 0
    for s in SS:
        print(s)
        A = 1
        funcbq = [lambda x , t = sj : bernoulli(t,x) for sj in s]
        for j in range(len(s)):
            A *= math.comb(q[j]+1, s[j])
        ded = 0
        for k in range(len(s)):
            if s[k] == 1:
                ded += ut.Dedekind(v, k, funcbq)        
        final += A*ded
    Q = 1
    for qj in q:
        Q *= 1/(qj+1)
    final *= Q
    print('recq = ', recq)
    print('derived value =', final)    
    return 0

def testcoro36(v,q):
    final = zet.bernrec(v, q)
    derived = zet.bernreczeta(v, q) - zet.deltabernrec(q)
    print('realval = ', final)
    print('derived value =', derived)
    print('ratio =', final/derived)     
    return 0

def testcoro37(v,q):
    check = 1 
    for k in range(len(q)):
        if sum([q[j] for j in range(len(q)) if j != k])%2 == 0:
            check = 0
    if check == 0:
        print('bad q, no all subsums are odd')
        return 0
    if check == 1:
        derived = - zet.deltabernrec(q)
    final = zet.bernrec(v, q)
    print('realval = ', final)
    print('derived value =', derived)
    if derived == 0:
        print('DIVISON BY 0',  final,'----',derived)
    else:
        print('ratio =', final/derived) 
    return 0

v=[2,3,7]
q = [2,2,1]
#k=1
#a0,a1,a2 = [1/2,1/3,1/5]
#testexamp321(v,a0,a1,a2)
#a=1/3
#testexamp322(v,a)
#testexamp351(v,k)
#testexamp352(v,k)
#testexamp353(v,k)
#testeq11(v,q)
testcoro36(v,q)
#testcoro37(v,q)