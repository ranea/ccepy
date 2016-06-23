import sys
sys.path.append('../ccepy')
import unittest
import doctest

from hypothesis import given, assume
from hypothesis.strategies import integers, sampled_from

from ccepy import curvas_elipticas  # para cargar los docstring
from ccepy.curvas_elipticas import curva_eliptica_sobre_Fq
from ccepy.listado_curvas_elipticas import curvas_eliptipcas_sobre_Fq_famosas
from ccepy.listado_curvas_elipticas import procesar_parametros_curva_eliptica


# TODO: ver que hacer con esto
# ParametrosCurvaElipticaSobreF2m = namedtuple('CurvaElípticaF2m', ['nombre', 'm', 'pol_irreducible', 'a', 'b', 'x1', 'y1'])
#
# # TODO: ref guide to elliptic curve
# curvas_eliptipcas_sobre_F2m_famosas = [
#     ParametrosCurvaElipticaSobreF2m(
#         "B-163",
#         163,
#         "z^163 + z^7 + z^6 + z^3 + 1",
#         1,
#         0x000000020A601907B8C953CA1481EB10512F78744A3205FD,
#         0x00000003F0EBA16286A2D57EA0991168D4994637E8343E36,
#         0x00000000D51FBC6C71A0094FA2CDD545B11C5C0C797324F1
#     ),
# ]


class TestCurvaElipticaFq(unittest.TestCase):
    """Conjuto de test para PuntosFqRacionales"""
    @given(sampled_from(curvas_eliptipcas_sobre_Fq_famosas), integers(), integers(), integers())
    def test_propiedades_grupo(self, ce, k1, k2, k3):
        assume(k1 >= 0)
        assume(k2 >= 0)
        assume(k3 >= 0)

        E, generador, orden = procesar_parametros_curva_eliptica(ce)
        # print(ce.nombre)

        P = generador * k1
        Q = generador * k2
        R = generador * k3
        O = E.elemento_neutro()

        # print("\tP: {0}\n\tQ: {1}\n\tR: {2}".format(P, Q, R))
        assert P + O == P == O + P
        assert P + (-P) == O
        assert P + Q == Q + P
        assert P + (Q + R) == (P + Q) + R

    @given(sampled_from(curvas_eliptipcas_sobre_Fq_famosas), integers(), integers(min_value=-100, max_value=100))
    def test_multiplicacion_por_duplicacion(self, ce, k, e):
        assume(k >= 0)
        assume(e)

        E, generador, orden = procesar_parametros_curva_eliptica(ce)
        # print(ce.nombre)

        P = generador * k

        multiplicacion = E.elemento_neutro()
        if e < 0:
            for i in range(-e):
                multiplicacion += -P
        elif e > 0:
            for i in range(e):
                multiplicacion += P

        # print("\tP: {0}\n\tP': {1}".format(P, multiplicacion))
        assert P * e == multiplicacion


# TODO: decidir si poner o no
# class TestCurvaElipticaF2m(unittest.TestCase):
#     """Conjuto de test para PuntosF2mRacionales"""
#     @given(sampled_from(curvas_eliptipcas_sobre_F2m_famosas))
#     def setUp(self, parametros_curva_eliptica):
#         nombre = parametros_curva_eliptica.nombre
#         m = parametros_curva_eliptica.m
#
#         rep_pol_irreducible = parametros_curva_eliptica.pol_irreducible
#         monomios = rep_pol_irreducible.split(" + ")
#         pol_irreducible = PolinomioZp([0], p=2)
#         for monomio in monomios:
#             if monomio == "1":
#                 pol_irreducible += PolinomioZp([1], p=2)
#                 continue
#
#             grado = int(monomio.split("^")[1])
#             pol_irreducible += PolinomioZp.monomio(1, grado, 2)
#
#         a = parametros_curva_eliptica.a
#
#         F2m = Fq(2, m, pol_irreducible)
#
#         b = parametros_curva_eliptica.b
#         b = F2m([int(i) for i in reversed(bin(b)[2:])])
#         x1 = parametros_curva_eliptica.x1
#         x1 = F2m([int(i) for i in reversed(bin(x1)[2:])])
#         y1 = parametros_curva_eliptica.y1
#         y1 = F2m([int(i) for i in reversed(bin(y1)[2:])])
#
#         self.nombre, self.m, self.pol_irreducible = nombre, m, pol_irreducible
#         self.a, self.b, self.x1, self.y1 = a, b, x1, y1
#
#     @given(integers(min_value=0, max_value=3), integers(min_value=0, max_value=3), integers(min_value=0, max_value=3))
#     def test_propiedades_grupo(self, k1, k2, k3):
#         assume(k1 >= 0)
#         assume(k2 >= 0)
#         assume(k3 >= 0)
#
#         nombre, m, pol_irreducible = self.nombre, self.m, self.pol_irreducible
#         a, b, x1, y1 = self.a, self.b, self.x1, self.y1
#
#         E = curva_eliptica_sobre_F2m(a, b, m, pol_irreducible)
#         generador = E(x1, y1)
#         print(nombre)
#
#         P = generador * k1
#         Q = generador * k2
#         R = generador * k3
#         O = E.elemento_neutro()
#
#         print("\tP: {0}\n\tQ: {1}\n\tR: {2}".format(P, Q, R))
#         assert P + O == P == O + P
#         assert P + (-P) == O
#         assert P + Q == Q + P
#         assert P + (Q + R) == (P + Q) + R
#
#     @given(integers(min_value=-10, max_value=10), integers(min_value=-10, max_value=10))
#     def test_multiplicacion_por_duplicacion(self, k, e):
#         assume(k >= 0)
#         assume(e)
#
#         nombre, m, pol_irreducible = self.nombre, self.m, self.pol_irreducible
#         a, b, x1, y1 = self.a, self.b, self.x1, self.y1
#
#         E = curva_eliptica_sobre_F2m(a, b, m, pol_irreducible)
#         generador = E(x1, y1)
#         print(nombre)
#
#         P = generador * k
#
#         multiplicacion = E.elemento_neutro()
#         if e < 0:
#             for i in range(-e):
#                 multiplicacion += -P
#         elif e > 0:
#             for i in range(e):
#                 multiplicacion += P
#
#         print("\tP: {0}\n\tP': {1}".format(P, multiplicacion))
#         assert P * e == multiplicacion


def load_tests(loader, tests, ignore):
    """Añade los ejemplos insertados en los docstring."""
    tests.addTests(doctest.DocTestSuite(curvas_elipticas))
    return tests


if __name__ == '__main__':
    unittest.main()
