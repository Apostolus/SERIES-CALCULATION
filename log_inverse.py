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
    I = inverse(1 + f, n)
    return (D.multiplication_trunc(I, n - 1)).integral()


def inverse2(f, n):
    # calcul pour une serie formelle de la forme F=1+XG de son inverse modulo N
    if n == 1:
        return f.parent()(1 / f)
    else:
        G = inverse(f, (n + 1) // 2)
        return G + (1 - G.multiplication_trunc(f, n)).multiplication_trunc(G, n)


def log2(f, n):
    a = f[0]
    return log(a) + log1((f / a - 1), n)


def exp(f, n):
    s = f.parent().one()
    for k in range(int(log(n))):
        s += s.multiplication_trunc(f - log2(s, n), n)
    return s
