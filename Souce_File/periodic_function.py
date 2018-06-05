from sage.all import *


def inverse(f, n):
    # calcul pour une serie formelle de la forme F=1+XG de son inverse modulo N
    if n == 1:
        return f.parent().one()
    else:
        G = inverse(f, (n + 1) // 2)
        return G + (1 - G.multiplication_trunc(f, n)).multiplication_trunc(G, n)


def arctan_series(f,n) :
    #f(0) doit être différent de 0
    g = f.derivative()
    h = f.parent().one()+f.multiplication_trunc(f,n)
    return (g.multiplication_trunc(inverse(h,n-1),n-1)).integral()


def tan_series(f,n) :
    #f(0) doit être différent de 0
    s = f
    for k in range(log(n,2)) :
        t = (f.parent().one()+s.multiplication_trunc(s,n))
        s += (f-arctan_series(s,n)).multiplication_trunc(t,n)
    return s

def sin_series(f,n) :
    t = tan(f/2,n)
    return (2*t).multiplication_trunc(inverse(f.parent().one()+t.multiplication_trunc(t,n),n),n)

def cos_series(f,n) :
    t = tan(f/2,n)
    t = t.multiplication_trunc(t,n)
    one = f.parent().one()
    return (one-t).multiplication_trunc(inverse(one+t,n),n)
