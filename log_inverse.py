from sage.all import *


def inverse(f, n):
    # calcul pour une serie formelle de la forme F=1+XG de son inverse modulo N
    if n == 1:
        return f.parent().one()
    else:
        G = inverse(f, (n + 1) // 2)
        return G + (1 - G.multiplication_trunc(f, n)).multiplication_trunc(G, n)


def log1(f, n):
	"""
	EXAMPLES ::
	sage: R=PolynomialRing(QQ,'x')
	sage: x=R.gen()
	sage: R=x
	sage: log1(R,10)
	1/9*x^9 - 1/8*x^8 + 1/7*x^7 - 1/6*x^6 + 1/5*x^5 - 1/4*x^4 + 1/3*x^3 - 1/2*x^2 + x

	sage: R=x
	sage: exp(R,10)
	sage: 1/362880*x^9 + 1/40320*x^8 + 1/5040*x^7 + 1/720*x^6 + 1/120*x^5 + 1/24*x^4 + 1/6*x^3 + 1/2*x^2 + x + 1
	
	sage: 
	sage: 
	sage: 
	
	"""
    # calcul pour une serie formelle de la forme F=XG de log(1+F) modulo N
    D = f.derivative()
    I = inverse(1 + f, n)
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
    s = f.parent().one()
    for k in range(n):
        s += s.multiplication_trunc(f - log2(s, n), n)
    return s

def fsin(f,n):
	t=exp(f*CC.0,n)
	return -0.5*CC.0*(t-inverse2(t,n))

def fcos(f,n):
	t=exp(f*CC.0,n)
	return 0.5*(t+inverse2(t,n))

def ftan(f,n):
	return fsin(f,n)*inverse2(fcos(f,n),n)

def farctan(f,n):
	return CC.0/2*inverse2(log2((CC.0+f).multiplication_trunc(CC.0-f,n),n),n)
