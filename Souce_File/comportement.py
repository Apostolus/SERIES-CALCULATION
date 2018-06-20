from sage.all import *
from series import *
from random import *
from timeit import *
from time import *
from math import *

Q=PolynomialRing(QQ,'x')
x=Q.gen()
R=PolynomialRing(RR,'y')
y=R.gen()
C=PolynomialRing(CC,'z')
z=C.gen()
AQ=PolynomialRing(Q,'a')
a=AQ.gen()


def poly_random(switch,nbTermes,maxInt):
    """
    C'est une fonction qui revoie un polynome random avec des restrictions donnees:
    switch = le code de mode: 1->rationnel, 2->reel, 3->imaginaire, 4->anneau de l'anneau quotient
    nbTermes = le nombre de termes dans un polynome
    maxInt = la valeur absolue maximale de numerateur et denominateur pour le mode rationnel
             ou celle de coefficiences pour le mode reel
             ou celle de parties reelle et imaginaire de coefficience pour le mode imaginaire
    """
    P = 0

    if(switch == 1):
        for k in range(nbTermes):
            P = P*x + int(2*maxInt*random()-maxInt)/int(2*maxInt*random()-maxInt+1)
    elif(switch == 2):
        for k in range(nbTermes):
            P = P*y + 2*maxInt*random()-maxInt
    elif switch == 3: 
        for k in range(nbTermes):
            P = P*z + 2*maxInt*random()-maxInt + (2*maxInt*random()-maxInt)*i
    else:
        print("Mode non-connu")
        sys.exit(1)

    return P

def comp_inverses_un(switch,nbTermes,maxInt,precision):
    P = poly_random(switch,nbTermes,maxInt)
    P = P - P(0) + P.parent().one()

    return inverse_un_series(P,precision)


def test():
    print("Test de inverse_un_series, 10 termes, maxInt=80, et 12 fois precision")
    timeit('comp_inverses_un(1,10,80,12)',number = 20)
