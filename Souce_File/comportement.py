from sage.all import *
from series import *
from random import *
from time import *
from math import *
import matplotlib.pyplot as plt


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
            P = P*x + (int(2*maxInt*random())-maxInt)/(int(maxInt*random())+1)
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

def tempsMoyen(nb,fun,*paras):
    """
    C'est une fonction qui teste le temps moyen d'effectuer la fonction 'fun(*paras)' avec 'nb' fois des iterations
    """
    moyen = 0
    for k in range(nb):
        debut = time()
        fun(*paras)
        fin = time()
        moyen = moyen + fin - debut

    return moyen/nb

def comp_inverse_un(switch,nbTermes,maxInt,precision):
    """
    C'est une fonction qui change le polynome random pour que la fonction 'inverse_un_series' puis le prendre comme un parametre, et effectuera cette fonction
    """
    P = poly_random(switch,nbTermes,maxInt)
    P = P - P(0) + P.parent().one()
    return inverse_un_series(P,precision)


def test():
    """
    C'est une fonction qui teste le comportement de toutes les fonctions demandees dans le fichier series.py
    """
    maxInt = 100
    nb = 100
    fois = 5
    moyen = 0
    file_abs = ""

    #---------------------------------inverse_un------------------------------------#
    print("Test de inverse_un_series, maxInt=%d\n" %(maxInt))

    ##----------------variable = nombre de termes-------------------##
    print("\tinverse_un TERMES")
    file_abs = "/usr/share/SageMath/2I013/data/inverse_un_terme.txt"
    f = open(file_abs,"w")
    for k in range(1,nb+1):
        print("\t\t%d termes, precision fixee = 12, on effectuera %d fois des iterations" %(k,fois))
        moyen = tempsMoyen(fois,comp_inverse_un,1,k,maxInt,12)
        f.write("%d %.15f\n" %(k,moyen))
        print("\t\t  Temps moyen = %.15f seconds" %(moyen))
    f.close()
    print("\n")

    ##---------------variable = niveau de precision-----------------##
    print("\tinverse_un PRECISION")
    file_abs = "/usr/share/SageMath/2I013/data/inverse_un_precision.txt"
    f = open(file_abs,"w")
    for k in range(1,nb+1):
        print("\t\t%d fois de precision, nombre de termes fixe = 20, on effectuera %d fois des iterations" %(k,fois))
        moyen = tempsMoyen(fois,comp_inverse_un,1,20,maxInt,k)
        f.write("%d %.15f\n" %(k,moyen))
        print("\t\t  Temps moyen = %.15f seconds" %(moyen))
    f.close()
    print("\n")
