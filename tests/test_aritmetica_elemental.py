import sys
sys.path.append('../ccepy')
import unittest
from math import gcd
import doctest

from hypothesis import given, assume
from hypothesis.strategies import integers, lists, sampled_from

from ccepy import aritmetica_elemental  # para cargar los docstring
from ccepy.aritmetica_elemental import PolinomioZp, Zp, alg_euclides, alg_euclides_polinomios

# secuencia A000040 de OEIS
primos = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
    61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131,
    137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197,
    199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271
]


# http://mathworld.wolfram.com/IrreduciblePolynomial.html
irreducibles_Z2_hasta_grado_4 = [
    PolinomioZp([1, 1], p=2),
    PolinomioZp([0, 1], p=2),
    PolinomioZp([1, 1, 1], p=2),
    PolinomioZp([1, 1, 0, 1], p=2),
    PolinomioZp([1, 0, 1, 1], p=2),
    PolinomioZp([1, 1, 0, 0, 1], p=2),
    PolinomioZp([1, 1, 1, 1, 1], p=2),
    PolinomioZp([1, 0, 0, 1, 1], p=2)
]


class TestEnteroModuloP(unittest.TestCase):
    """Conjuto de test para EnteroModuloP"""
    @given(sampled_from(primos), integers(), integers(), integers())
    def test_propiedades_cuerpo(self, q, n, m, k):
        assume(n)  # para que sea un entero
        assume(m)
        assume(k)
        Zq = Zp(q)
        x = Zq(n)
        y = Zq(m)
        t = Zq(k)
        cero = Zq(0)
        uno = Zq(1)
        # propiedades de anillo conmutativo
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
        if x != 0:
            assert x * x.inverso() == uno

    @given(sampled_from(primos), integers(), integers(), integers())
    def test_propiedades_potencias(self, q, n, e, f):
        assume(n)
        assume(e)
        assume(f)
        Zq = Zp(q)
        x = Zq(n)
        assume(not(x == 0 and (e < 0 or f < 0)))  # para evitar hacer 0 ** (-1)
        assert x ** e * x ** f == x ** (e + f)
        assert (x ** e) ** f == x ** (e * f)
        if x != 0:
            assert x ** e / x ** f == x ** (e - f)


class TestPolinomioZp(unittest.TestCase):
    """Conjuto de test para la clase PolinomioZp"""
    @given(lists(integers()), lists(integers()), lists(integers()), sampled_from(primos))
    def test_propiedades_anillo(self, l1, l2, l3, primo):
        assume(l1)  # para que no sea None ni []
        assume(l2)
        assume(l3)
        p = PolinomioZp(l1, primo)
        q = PolinomioZp(l2, primo)
        r = PolinomioZp(l3, primo)
        cero = PolinomioZp([0], primo)
        uno = PolinomioZp([1], primo)
        assert p + (q + r) == (p + q) + r
        assert p + q == q + p
        assert p + cero == p == cero + p
        assert p + (-p) == cero
        assert p * (q * r) == (p * q) * r
        assert p * uno == p
        assert p * (q + r) == (p * q) + (p * r)
        assert (p + q) * r == (p * r) + (q * r)
        assert p * q == q * p

    @given(lists(integers()), lists(integers()), sampled_from(primos))
    def test_composicion_suma_resta(self, l1, l2, primo):
        assume(l1)
        assume(l2)
        p = PolinomioZp(l1, primo)
        q = PolinomioZp(l2, primo)
        assert p == (p - q) + q

    @given(lists(integers()), lists(integers()), sampled_from(primos))
    def test_composicion_multiplicacion_division(self, l1, l2, primo):
        assume(l1)
        assume(l2)
        p = PolinomioZp(l1, primo)
        q = PolinomioZp(l2, primo)
        cero = PolinomioZp([0], primo)
        assume(q != cero)
        assert p == (p * q) / q

    @given(lists(integers()), integers(), sampled_from(primos))
    def test_multiplicacion_escalares(self, l1, n, primo):
        assume(l1)
        assume(n)
        p = PolinomioZp(l1, primo)
        assert p * n == n * p

    @given(lists(integers()), integers(min_value=0, max_value=3), integers(min_value=0, max_value=3), sampled_from(primos))
    def test_propiedades_potencias(self, l1, e, f, primo):
        assume(l1)
        assume(e)
        assume(f)
        p = PolinomioZp(l1, primo)
        assert p ** e * p ** f == p ** (e + f)
        assert (p ** e) ** f == p ** (e * f)

    @given(lists(integers()), lists(integers()), sampled_from(primos))
    def test_grado(self, l1, l2, primo):
        assume(l1)
        assume(l2)
        p = PolinomioZp(l1, primo)
        q = PolinomioZp(l2, primo)
        assert ((p * q).grado()) == p.grado() + q.grado()
        assert (p + q).grado() <= max(p.grado(), q.grado())

    @given(integers(min_value=1, max_value=4))
    def test_polinomios_irreducibles_Z2(self, n):
        assume(n)
        f = PolinomioZp.genera_irreducible(grado=n, p=2)
        assert f in irreducibles_Z2_hasta_grado_4


class TestAlgoritmoExtendidoEuclides(unittest.TestCase):
    """Conjuto de test para los algoritmos extendidos de euclides"""
    @given(integers(min_value=1), integers(min_value=1))
    def test_alg_euclides(self, a, b):
        x, y, d = alg_euclides(a, b)
        assert x * a + y * b == d
        assert d == gcd(a, b)
        if (a // b) * b != a and (b // a) * a != b:  # a no divide a b y viceversa
            assert x < b // d
            assert y < a // d

    @given(lists(integers()), lists(integers()), sampled_from(primos))
    def test_alg_euclides_polinomios(self, l1, l2, primo):
        assume(l1)
        assume(l2)
        g = PolinomioZp(l1, primo)
        h = PolinomioZp(l2, primo)
        cero = PolinomioZp([0], primo)
        assume(g != cero)
        s, t, d = alg_euclides_polinomios(g, h, p=primo)
        assert s * g + t * h == d
        assert g % d == 0 and h % d == 0  # vemos si el gcd divide a ambos
        if h != cero:
            assert s.grado() <= h.grado() and t.grado() <= g.grado()


def load_tests(loader, tests, ignore):
    """Añade los ejemplos insertados en los docstring."""
    tests.addTests(doctest.DocTestSuite(aritmetica_elemental))
    return tests

if __name__ == '__main__':
    unittest.main()
