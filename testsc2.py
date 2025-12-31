import utils as ut
import conesfunc as con


def testexample26(v,k):
    r = len(v)-1
    def power(v,k):
        suma = 0
        for i in range(1,v[k]):
            prod = 1
            for j in range(0,r+1):
                if j != k:
                    prod = prod*((i*v[j])%v[k])
            suma = suma + prod
        return suma
    final = ut.Dedekind(v, k, [ut.frac]*(r+1))*(v[k]**r)
    derived = power(v,k)
    #print('realval = ', final)
    #print('derived value =', derived)
    if derived == 0:
        print('DIVISON BY 0',  final,'----',derived)
    else:
        print('ratio =', final/derived) 
    return 0

def testeq4(L,q):
    normq = sum(q)
    val1 = con.zetacone(L,q)
    val2 = con.zetaconegcd(L,q)
    val3 = ut.zeta(normq)
    final = val1
    derived =val2*val3
    print('realval = ', final)
    print('derived value =', derived)
    print('ratio =', final/derived)
    return 0

def testremark27(v,q):
    L = [v]
    normq = sum(q)
    final = con.zetacone(L,q)
    derived = ut.zeta(normq)
    for i in range(len(v)):
        derived *= v[i]**(-q[i])
    print('realval = ', final)
    print('derived value =', derived)
    print('ratio =', final/derived)
    return 0

def testsexam28(r,q):
    if r != len(q):
        return 'error'
    can = []
    for j in range(0,r):
        e = [0]*r
        e[j]=1 
        can.append(e)
    final = con.zetacone(can,q)
    derived = 1 
    for j in range(len(q)):
        derived *= ut.zeta(q[j])
    print('realval = ', final)
    print('derived value =', derived)
    print('ratio =', final/derived)
    return 0

v=[2,3,5]
q = [2,2,4]

for k in range(0,len(v)):
    testexample26(v,k)

#L = [[1,0,-1],[0,1,-1]]
#testeq4(L,q)

#testremark27(v,q)

#r=3
#testsexam28(r,q)