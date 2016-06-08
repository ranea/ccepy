import sys
sys.path.append('../ccepy')
import unittest
# import doctest

from hypothesis import given, assume  # , settings
from hypothesis.strategies import integers, lists, sampled_from

# from ccepy import cuerpos_finitos
from ccepy.cuerpos_finitos import Fq

# TODO: decidir que hacer con los primos grandes
# secuencia A000040
primos = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
    # 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131,
    # 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197,
    # 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271
]


# TODO: rellenar los docstring?
class TestElementoFpn(unittest.TestCase):
    """Conjuto de test para ElementoFpn"""
    # TODO: quitar debug
    # @settings(perform_health_check=False)
    @given(sampled_from(primos), integers(min_value=2, max_value=3), lists(integers()),
        lists(integers()), lists(integers()))
    def test_propiedades_cuerpo(self, p, n, l1, l2, l3):
        assume(n)
        assume(l1)
        assume(l2)
        assume(l3)
        # print("p: {0}, n: {1}".format(p, n))
        Fpn = Fq(p, n)
        x = Fpn(l1)
        y = Fpn(l2)
        t = Fpn(l3)
        cero = Fpn.cero()
        uno = Fpn.uno()
        # print("{0}\n\tx: {1}\n\ty: {2}\n\tt: {3}".format(Fpn.__name__, x, y, t))
        # propieadades de anillo conmutativo
        assert x + (y + t) == (x + y) + t
        assert x + y == y + x
        assert x + cero == x == cero + x
        assert x + (-x) == cero
        assert x * (y * t) == (x * y) * t
        assert x * uno == x
        assert x * (y + t) == (x * y) + (x * t)
        assert (x + y) * t == (x * t) + (y * t)
        assert x * y == y * x
        # propiedades de cuerpo
        if x != Fpn.cero():
            assert x * x.inverso() == uno

    # TODO: quitar debug
    # @settings(perform_health_check=False)
    @given(sampled_from(primos), integers(min_value=2, max_value=3), lists(integers()),
        integers(min_value=-10, max_value=10), integers(min_value=-10, max_value=10))
    def test_propiedades_potencias(self, p, n, l1, e, f):
        assume(n)
        assume(l1)
        assume(e)
        assume(f)
        Fpn = Fq(p, n)
        x = Fpn(l1)
        # print("{0}\te: {1}\tf: {2}\tx: {3}".format(Fpn.__name__, e, f, x))
        assume(not (x == Fpn.cero() and (e < 0 or f < 0)))  # para evitar hacer 0 ** (-1)
        assert x ** e * x ** f == x ** (e + f)
        assert (x ** e) ** f == x ** (e * f)
        if x != Fpn.cero():
            assert x ** e / x ** f == x ** (e - f)


# def load_tests(loader, tests, ignore):
#     tests.addTests(doctest.DocTestSuite(cuerpos_finitos))
#     return tests

if __name__ == '__main__':
    unittest.main()
