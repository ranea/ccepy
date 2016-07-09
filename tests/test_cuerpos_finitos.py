import sys
sys.path.append('../ccepy')
import unittest

from hypothesis import given, assume
from hypothesis.strategies import integers, lists, sampled_from

from ccepy.cuerpos_finitos import Fq


# parte de la secuencia A000040 de OEIS
primos = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
]


class TestElementoFpn(unittest.TestCase):
    """Conjuto de test para ElementoFpn"""
    @given(sampled_from(primos), integers(min_value=2, max_value=3), lists(integers()),
        lists(integers()), lists(integers()))
    def test_propiedades_cuerpo(self, p, n, l1, l2, l3):
        assume(n)
        assume(l1)
        assume(l2)
        assume(l3)
        Fpn = Fq(p, n)
        x = Fpn(l1)
        y = Fpn(l2)
        t = Fpn(l3)
        cero = Fpn.cero()
        uno = Fpn.uno()
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

    @given(sampled_from(primos), integers(min_value=2, max_value=3), lists(integers()),
        integers(min_value=-10, max_value=10), integers(min_value=-10, max_value=10))
    def test_propiedades_potencias(self, p, n, l1, e, f):
        assume(n)
        assume(l1)
        assume(e)
        assume(f)
        Fpn = Fq(p, n)
        x = Fpn(l1)
        assume(not (x == Fpn.cero() and (e < 0 or f < 0)))  # para evitar hacer 0 ** (-1)
        assert x ** e * x ** f == x ** (e + f)
        assert (x ** e) ** f == x ** (e * f)
        if x != Fpn.cero():
            assert x ** e / x ** f == x ** (e - f)


if __name__ == '__main__':
    unittest.main()
