import sys
sys.path.append('../ccepy')
import unittest
import doctest

from hypothesis import given, assume
from hypothesis.strategies import integers, sampled_from

from ccepy import curvas_elipticas  # para cargar los docstring
from ccepy.listado_curvas_elipticas import curvas_eliptipcas_sobre_Fq_famosas
from ccepy.listado_curvas_elipticas import procesar_parametros_curva_eliptica


class TestCurvaElipticaFq(unittest.TestCase):
    """Conjuto de test para PuntosFqRacionales"""
    @given(sampled_from(curvas_eliptipcas_sobre_Fq_famosas), integers(), integers(), integers())
    def test_propiedades_grupo(self, ce, k1, k2, k3):
        assume(k1 >= 0)
        assume(k2 >= 0)
        assume(k3 >= 0)

        E, generador, orden = procesar_parametros_curva_eliptica(ce)

        P = generador * k1
        Q = generador * k2
        R = generador * k3
        O = E.elemento_neutro()

        assert P + O == P == O + P
        assert P + (-P) == O
        assert P + Q == Q + P
        assert P + (Q + R) == (P + Q) + R

    @given(sampled_from(curvas_eliptipcas_sobre_Fq_famosas), integers(), integers(min_value=-100, max_value=100))
    def test_multiplicacion_por_duplicacion(self, ce, k, e):
        assume(k >= 0)
        assume(e)

        E, generador, orden = procesar_parametros_curva_eliptica(ce)

        P = generador * k

        multiplicacion = E.elemento_neutro()
        if e < 0:
            for i in range(-e):
                multiplicacion += -P
        elif e > 0:
            for i in range(e):
                multiplicacion += P

        assert P * e == multiplicacion


def load_tests(loader, tests, ignore):
    """Añade los ejemplos insertados en los docstring."""
    tests.addTests(doctest.DocTestSuite(curvas_elipticas))
    return tests


if __name__ == '__main__':
    unittest.main()
