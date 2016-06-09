import sys
sys.path.append('../ccepy')
import unittest
# import doctest

from hypothesis import given, assume  # , settings
from hypothesis.strategies import integers, lists, sampled_from

# from ccepy import cuerpos_finitos
from ccepy.curvas_elipticas import curva_eliptica_sobre_Fq
from ccepy.cuerpos_finitos import Fq

if __name__ == '__main__':
    # E = curva_eliptica_sobre_Q(0, 4)
    # print(E.coeficientes)
    # print("Discriminante:", E.discriminante)
    # P = E(0, -2)
    # print("P =", P)
    # print("-P =", -P)
    # for i in range(3):
    #     print(i, "P =", i * P)
    #
    # E = curva_eliptica_sobre_Q(-18, -72)
    # print(E.coeficientes)
    # print("Discriminante:", E.discriminante)
    # P = E(6, 6)
    # print("P =", P)
    # print("-P =", -P)
    # for i in range(5):
    #     print(i, "P =", i * P)

    # E = curva_eliptica_sobre_Fq(2, 3, 97, 1)
    # print(E.coeficientes)
    # print("Discriminante:", E.discriminante)
    # F97 = Fq(97, 1)
    # P = E(0, 10)
    # print("P =", P)
    # print("-P =", -P)
    # for i in range(5):
    #     print(i, "P =", i * P)
    # Q = E(3, 6)
    # print("Q =", Q)
    # print("P + Q =", P + Q)

    E = curva_eliptica_sobre_Fq(1, 1, 5, 2)
    print(E.coeficientes)
    print("Discriminante:", E.discriminante)
    F25 = Fq(5, 2)
    P = E(F25.cero(), F25.uno())
    print("P =", P)
    print("-P =", -P)
    for i in range(5):
        print(i, "P =", i * P)
