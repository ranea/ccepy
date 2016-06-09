"""Curvas elípticas.
...
"""
import copy
from fractions import Fraction
from collections import namedtuple
from abc import ABCMeta, abstractmethod

EcuacionWeierstrass = namedtuple('Coeficientes', ['a', 'b'])

from ccepy.cuerpos_finitos import Fq


class PuntosRacionales(metaclass=ABCMeta):
    """..."""
    __slots__ = ('_x', '_y')  # para mejorar la eficiencia si hay muchos objetos

    @abstractmethod
    def __init__(self, x, y):
        """Inicializa self._x y self._y."""

    def es_elemento_neutro(self):
        return self._x is None or self._y is None

    # TODO: explicar en algún lado el uso de property
    @property
    def x(self):
        if self.es_elemento_neutro():
            raise AttributeError("El elemento neutro no tiene componente x")
        else:
            return self._x

    @property
    def y(self):
        if self.es_elemento_neutro():
            raise AttributeError("El elemento neutro no tiene componente y")
        else:
            return self._y

    @classmethod
    def elemento_neutro(cls):
        return cls(None, None)

    @abstractmethod
    def __eq__(self, other):
        return

    def __ne__(self, other):
        return not self.__eq__(other)

    @abstractmethod
    def __add__(self, other):
        return

    @abstractmethod
    def __neg__(self):
        return

    def __sub__(self, other):
        return self + (-other)

    @abstractmethod
    def __mul__(self, other):
        return

    @abstractmethod
    def __rmul__(self, other):
        return

    def __str__(self):
        if self.es_elemento_neutro():
            return "Elemento neutro"
        else:
            return "({0},{1})".format(self.x, self.y)


def curva_eliptica_sobre_Q(a, b):
    class PuntosQRacionales(PuntosRacionales):
        """..."""
        # coeficientes (a, b) de la ecuación y^2 = x^3 + a*x + b
        coeficientes = None
        discriminante = None

        @classmethod
        def contiene(cls, x, y):
            a, b = cls.coeficientes
            lado_izquierdo_ecuacion = y**2
            lado_derecho_ecuacion = x**3 + a * x + b
            return lado_izquierdo_ecuacion == lado_derecho_ecuacion

        def __init__(self, x, y):
            if x is None or y is None:
                self._x = None
                self._y = None
            else:
                if PuntosQRacionales.contiene(x, y):
                    self._x = Fraction(x)  # para aceptar también int
                    self._y = Fraction(y)
                else:
                    raise ValueError("El punto ({0}, {1}) no pertenece a la" +
                                    " curva".format(x, y))

        def __eq__(self, other):
            if self.es_elemento_neutro():
                return other.es_elemento_neutro()
            elif other.es_elemento_neutro():
                return False

            return self.x == other.x and self.y == other.y

        def __add__(self, other):
            if self.es_elemento_neutro():
                return other
            elif other.es_elemento_neutro():
                return self

            x1, y1 = self.x, self.y
            x2, y2 = other.x, other.y
            a = PuntosQRacionales.coeficientes.a

            if self == other:
                if y1 == 0:
                    return PuntosQRacionales.elemento_neutro()
                else:
                    m = (3 * x1**2 + a) / (2 * y1)
                    x3 = m**2 - 2 * x1
                    y3 = m * (x1 - x3) - y1
                    return PuntosQRacionales(x3, y3)
            elif x1 == x2:
                # y1 != y2
                return PuntosQRacionales.elemento_neutro()
            else:
                m = (y2 - y1) / (x2 - x1)
                x3 = m**2 - x1 - x2
                y3 = m * (x1 - x3) - y1
                return PuntosQRacionales(x3, y3)

        def __neg__(self):
            return PuntosQRacionales(self.x, -self.y)

        def __mul__(self, entero):
            producto = PuntosQRacionales.elemento_neutro()  # una copia
            for i in range(entero):
                producto += self
            return producto

        __rmul__ = __mul__

    discriminante = 4 * a**3 + 27 * b**2
    if discriminante == 0:
        raise ValueError("El discriminant, 4a^3 + 27b^2, no puede ser cero.")

    PuntosQRacionales.discriminante = discriminante
    PuntosQRacionales.coeficientes = EcuacionWeierstrass(a, b)
    return PuntosQRacionales


def curva_eliptica_sobre_Fq(a, b, p, n=1, pol_irreducible=None):
    class PuntosFqRacionales(PuntosRacionales):
        """..."""
        # coeficientes (a, b) de la ecuación y^2 = x^3 + a*x + b
        coeficientes = None
        discriminante = None
        Fq = None

        @classmethod
        def contiene(cls, x, y):
            a, b = cls.coeficientes
            lado_izquierdo_ecuacion = y**2
            lado_derecho_ecuacion = x**3 + a * x + b
            return lado_izquierdo_ecuacion == lado_derecho_ecuacion

        def __init__(self, x, y):
            if x is None or y is None:
                self._x = None
                self._y = None
            else:
                self._x = PuntosFqRacionales.Fq(x)
                self._y = PuntosFqRacionales.Fq(y)
                if not PuntosFqRacionales.contiene(self._x, self._y):
                    raise ValueError("El punto ({0}, {1}) no pertenece a la" +
                                    " curva".format(x, y))

        def __eq__(self, other):
            if self.es_elemento_neutro():
                return other.es_elemento_neutro()
            elif other.es_elemento_neutro():
                return False

            return self.x == other.x and self.y == other.y

        def __add__(self, other):
            if self.es_elemento_neutro():
                return other
            elif other.es_elemento_neutro():
                return self

            x1, y1 = self.x, self.y
            x2, y2 = other.x, other.y
            a = PuntosFqRacionales.coeficientes.a
            E = PuntosFqRacionales

            if self == other:
                if y1 == E.Fq.cero():
                    return PuntosFqRacionales.elemento_neutro()
                else:
                    # TODO: pensar como solucionar
                    Fq_3 = E.Fq.uno() + E.Fq.uno() + E.Fq.uno()
                    Fq_2 = E.Fq.uno() + E.Fq.uno()
                    m = (Fq_3 * x1**2 + a) / (Fq_2 * y1)
                    x3 = m**2 - Fq_2 * x1
                    y3 = m * (x1 - x3) - y1
                    return PuntosFqRacionales(x3, y3)
            elif x1 == x2:
                # y1 != y2
                return PuntosFqRacionales.elemento_neutro()
            else:
                m = (y2 - y1) / (x2 - x1)
                x3 = m**2 - x1 - x2
                y3 = m * (x1 - x3) - y1
                return PuntosFqRacionales(x3, y3)

        def __neg__(self):
            return PuntosFqRacionales(self.x, -self.y)

        # añadir referencia (adaptación del handbook)
        @classmethod
        def _multiplicacion_por_duplicacion(cls, punto, k):
            rep_binaria_k = "".join(reversed(bin(k)[2:]))
            t = len(rep_binaria_k) - 1

            Q = PuntosFqRacionales.elemento_neutro()
            if k == 0:
                return Q
            P = copy.deepcopy(punto)
            if rep_binaria_k[0] == "1":
                Q = copy.deepcopy(P)
            for i in range(1, t + 1):
                P = P + P  # duplicar
                if rep_binaria_k[i] == "1":
                    Q = P + Q  # suar
            return Q

        def __mul__(self, entero):
            if self.es_elemento_neutro():
                return self
            else:
                return PuntosFqRacionales._multiplicacion_por_duplicacion(self, entero)

        __rmul__ = __mul__

    if p == 2 or p == 3:
        raise ValueError("p no puede ser ni 2 ni 3.")

    F_q = Fq(p, n, pol_irreducible)
    if n == 1:
        A = F_q(a)
        B = F_q(b)
        discriminante = 4 * A**3 + 27 * B**2
    else:
        # n > 2
        A = F_q([a])
        B = F_q([b])
        discriminante = F_q([4]) * A**3 + F_q([27]) * B**2
    if discriminante == F_q.cero():
        raise ValueError("El discriminant, 4a^3 + 27b^2, no puede ser cero.")

    PuntosFqRacionales.discriminante = discriminante
    PuntosFqRacionales.coeficientes = EcuacionWeierstrass(A, B)
    PuntosFqRacionales.Fq = F_q
    return PuntosFqRacionales
