from sage.all import *


def inverse(f, n):
    # calcul pour une serie formelle de la forme F=1+XG de son inverse modulo N
    if n == 1:
        return f.parent().one()
    else:
        G = inverse(f, (n + 1) // 2)
        return G + (1 - G.multiplication_trunc(f, n)).multiplication_trunc(G, n)


def log1(f, n):
    # calcul pour une serie formelle de la forme F=XG de log(1+F) modulo N
    D = f.derivative()
    I = inverse(f.parent().one() + f, n)
    return (D.multiplication_trunc(I, n - 1)).integral()


def inverse2(f, n):

    # calcul l'inverse une serie formelle F avec F(0)!=0
    if n == 1:
        return f.parent()(1 / f)
    else:
        G = inverse(f, (n + 1) // 2)
        return G + (1 - G.multiplication_trunc(f, n)).multiplication_trunc(G, n)


def log2(f, n):
    a = f[0]
    return log(a) + log1((f / a - 1), n)


def exp(f, n):
    """
    fonctionne aussi avec r=n et en mettant effectuant dans la boucle les calculs à la précision k+2

    """
    one = f.parent().one()
    s = one
    r = int(log(n,2))
    for k in range(r):
        t = f - log1(s-one,n)
        s += s.multiplication_trunc(t,n)
    return s

def arctan(f,n) :
    #f(0) doit être différent de 0
    g = f.derivative()
    h = f.parent().one()+f.multiplication_trunc(f,n)
    return (g.multiplication_trunc(inverse(h,n-1),n-1)).integral()
    

def tan(f,n) :
    #f(0) doit être différent de 0
    s = f
    for k in range(log(n,2)) :
        t = (f.parent().one()+s.multiplication_trunc(s,n))
        s += (f-arctan(s,n)).multiplication_trunc(t,n)
    return s

def sin(f,n) :
    t = tan(f/2,n)
    return (2*t).multiplication_trunc(inverse(f.parent().one()+t.multiplication_trunc(t,n),n),n)

def cos(f,n) :
    t = tan(f/2,n)
    t = t.multiplication_trunc(t,n)
    one = f.parent().one()
    return (one-t).multiplication_trunc(inverse(one+t,n),n)
