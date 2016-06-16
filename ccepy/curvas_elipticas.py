"""Aritmética con curvas elípticas.

Este módulo permite operar con el grupo de puntos de una curva elíptica.

Para utilizar las funciones y las clases de este módulo, debe importarlo
previamente: ::

    # reemplace ... por la función/clase que desea utilizar
    from ccepy.curvas_elipticas import ...

Para operar con puntos de una curva elíptica, use las funciones de la forma
``curva_eliptica_sobre_*`` y los operadores aritméticos habituales.

    >>> E = curva_eliptica_sobre_Fq(a=2, b=3, p=97)  # y^2 = x^3 + 2x + 3 sobre F97
    >>> E.coeficientes
    Coeficientes(a=2, b=3)
    >>> P = E(0, 10)
    >>> P
    (0,10)
    >>> Q = E(3, 6)
    >>> Q
    (3,6)
    >>> P + Q
    (85,71)
    >>> -P
    (0,87)
    >>> 3 * P
    (23,24)
"""
import copy
from fractions import Fraction
from collections import namedtuple
from abc import ABCMeta, abstractmethod

EcuacionWeierstrass = namedtuple('Coeficientes', ['a', 'b'])

from ccepy.cuerpos_finitos import Fq


class PuntoRacional(metaclass=ABCMeta):
    """..."""
    __slots__ = ('_x', '_y')  # para mejorar la eficiencia si hay muchos objetos

    @classmethod
    @abstractmethod
    def contiene(cls, x, y):
        """Comprueba si (x,y) esta en la curva."""
        return

    @abstractmethod
    def __init__(self, x, y):
        """Inicializa self._x y self._y."""
        return

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
    class PuntoQRacional(PuntoRacional):
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
                if PuntoQRacional.contiene(x, y):
                    self._x = Fraction(x)  # para aceptar también int
                    self._y = Fraction(y)
                else:
                    raise ValueError("El punto ({0}, {1})".format(x, y) +
                                    " no pertenece a la curva.")

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
            a = PuntoQRacional.coeficientes.a

            if self == other:
                if y1 == 0:
                    return PuntoQRacional.elemento_neutro()
                else:
                    m = (3 * x1**2 + a) / (2 * y1)
                    x3 = m**2 - 2 * x1
                    y3 = m * (x1 - x3) - y1
                    return PuntoQRacional(x3, y3)
            elif x1 == x2:
                # y1 != y2
                return PuntoQRacional.elemento_neutro()
            else:
                m = (y2 - y1) / (x2 - x1)
                x3 = m**2 - x1 - x2
                y3 = m * (x1 - x3) - y1
                return PuntoQRacional(x3, y3)

        def __neg__(self):
            if self.es_elemento_neutro():
                return self
            else:
                return PuntoQRacional(self.x, -self.y)

        def __mul__(self, entero):
            producto = PuntoQRacional.elemento_neutro()  # una copia
            for i in range(entero):
                producto += self
            return producto

        __rmul__ = __mul__

    discriminante = 4 * a**3 + 27 * b**2
    if discriminante == 0:
        raise ValueError("El discriminant, 4a^3 + 27b^2, no puede ser cero.")

    PuntoQRacional.discriminante = discriminante
    PuntoQRacional.coeficientes = EcuacionWeierstrass(a, b)
    return PuntoQRacional


def curva_eliptica_sobre_Fq(a, b, p, n=1, pol_irreducible=None):
    """Devuelve el constructor de puntos de una curva elíptica sobre
    un cuerpo finito de característica distinta de 2 y 3.

        >>> E = curva_eliptica_sobre_Fq(1, 1, 5, 2)  # y^2 = x^3 + x + 1 sobre F25
        >>> E
        <class 'ccepy.curvas_elipticas.curva_eliptica_sobre_Fq.<locals>.PuntoFqRacional'>
        >>> E(0, 1)
        ({[0, 0]; 25},{[1, 0]; 25})

    Los dos primeros argumentos (``a``, ``b``) son los coeficientes de la ecuación
    de Weierstrass simplificada: :math:`y^2 = x^3 + a x + b`. Estos valores
    pueden ser bien de tipo ``int`` o bien de tipo ``EnteroModuloP`` o ``ElementoFpn``
    según sea ``n`` uno o mayor que uno respectivamente.

    Los tres últimos argumentos (``p``, ``n``, ``pol_irreducible``) definen el cuerpo
    finito de p**n elementos sobre el que se define la curva eliptipca.

    Args:
        a : el coeficiente que acompaña a x en la ecuación de Weierstrass
        b : el término independiente de la ecuación de Weierstrass
        p (int): un número primo.
        n (Optional[int]): un número natural.
        pol_irreducible (Optional[PolinomioZp]): un polinomio de grado
            *n* irreducible.

    Return:
        PuntoFqRacional: la clase que representa los puntos de la curva elíptica.
    """
    class PuntoFqRacional(PuntoRacional):
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
                self._x = PuntoFqRacional.Fq(x)
                self._y = PuntoFqRacional.Fq(y)
                if not PuntoFqRacional.contiene(self._x, self._y):
                    raise ValueError("El punto ({0}, {1})".format(x, y) +
                                    " no pertenece a la curva.")

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
            a = PuntoFqRacional.coeficientes.a
            Fq = PuntoFqRacional.Fq

            if self == other:
                if y1 == Fq(0):
                    return PuntoFqRacional.elemento_neutro()
                else:
                    m = (Fq(3) * x1**2 + a) / (Fq(2) * y1)
                    x3 = m**2 - Fq(2) * x1
                    y3 = m * (x1 - x3) - y1
                    return PuntoFqRacional(x3, y3)
            elif x1 == x2:
                # y1 != y2
                return PuntoFqRacional.elemento_neutro()
            else:
                m = (y2 - y1) / (x2 - x1)
                x3 = m**2 - x1 - x2
                y3 = m * (x1 - x3) - y1
                return PuntoFqRacional(x3, y3)

        def __neg__(self):
            if self.es_elemento_neutro():
                return self
            else:
                return PuntoFqRacional(self.x, -self.y)

        # añadir referencia (adaptación del handbook)
        @classmethod
        def _multiplicacion_por_duplicacion(cls, punto, k):
            rep_binaria_k = "".join(reversed(bin(k)[2:]))
            t = len(rep_binaria_k) - 1

            Q = PuntoFqRacional.elemento_neutro()
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
            elif entero < 0:
                return PuntoFqRacional._multiplicacion_por_duplicacion(-self, -entero)
            else:
                return PuntoFqRacional._multiplicacion_por_duplicacion(self, entero)

        __rmul__ = __mul__

    if p == 2 or p == 3:
        raise ValueError("p no puede ser ni 2 ni 3.")

    F_q = Fq(p, n, pol_irreducible)
    A = F_q(a)
    B = F_q(b)
    discriminante = F_q(4) * A**3 + F_q(27) * B**2
    if discriminante == F_q.cero():
        raise ValueError("El discriminant, 4a^3 + 27b^2, no puede ser cero.")

    PuntoFqRacional.discriminante = discriminante
    PuntoFqRacional.coeficientes = EcuacionWeierstrass(A, B)
    PuntoFqRacional.Fq = F_q
    return PuntoFqRacional


# TODO: remover clase tras realizar documentación
class PuntoFqRacional(PuntoRacional):
    """Representa un punto de una curva elíptica sobre un cuerpo finito de
    característica distinta de 2 y 3.

        >>> E = curva_eliptica_sobre_Fq(1, 1, 5, 2)  # y^2 = x^3 + x + 1 sobre F25
        >>> E.coeficientes
        >>> F25 = Fq(5, 2)
        >>> P = E(F25.cero(), F25.uno())
        >>> P
        ({[0, 0]; 25},{[1, 0]; 25})
        >>> Q = E(F25([4, 0]), F25([2, 0]))
        >>> Q
        ({[4, 0]; 25},{[2, 0]; 25})
        >>> P + Q
        ({[2, 0]; 25},{[1, 0]; 25})
        >>> -P
        ({[0, 0]; 25},{[4, 0]; 25})
        >>> 4 * P
        ({[3, 0]; 25},{[4, 0]; 25})

    Soporta los operadores ``+``, ``-``, ``*`` con su significado habitual.

    Args:
        x: ...
        y: ...

    Attributes:
        coeficientes: ...
        discriminante: ...
        Fq: ...
    """
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
            self._x = PuntoFqRacional.Fq(x)
            self._y = PuntoFqRacional.Fq(y)
            if not PuntoFqRacional.contiene(self._x, self._y):
                raise ValueError("El punto ({0}, {1})".format(x, y) +
                                " no pertenece a la curva.")

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
        a = PuntoFqRacional.coeficientes.a
        Fq = PuntoFqRacional.Fq

        if self == other:
            if y1 == Fq(0):
                return PuntoFqRacional.elemento_neutro()
            else:
                m = (Fq(3) * x1**2 + a) / (Fq(2) * y1)
                x3 = m**2 - Fq(2) * x1
                y3 = m * (x1 - x3) - y1
                return PuntoFqRacional(x3, y3)
        elif x1 == x2:
            # y1 != y2
            return PuntoFqRacional.elemento_neutro()
        else:
            m = (y2 - y1) / (x2 - x1)
            x3 = m**2 - x1 - x2
            y3 = m * (x1 - x3) - y1
            return PuntoFqRacional(x3, y3)

    def __neg__(self):
        if self.es_elemento_neutro():
            return self
        else:
            return PuntoFqRacional(self.x, -self.y)

    # añadir referencia (3.27 guide to elliptic curve)
    @classmethod
    def _multiplicacion_por_duplicacion(cls, punto, k):
        rep_binaria_k = "".join(bin(k)[2:])  # (k_t, k_{t-1},..., k_0)
        Q = PuntoFqRacional.elemento_neutro()

        for k_i in rep_binaria_k:
            Q = Q + Q  # duplicar
            if k_i == "1":
                Q = Q + P  # sumar

        return Q

    def __mul__(self, entero):
        if self.es_elemento_neutro():
            return self
        elif entero < 0:
            return PuntoFqRacional._multiplicacion_por_duplicacion(-self, -entero)
        else:
            return PuntoFqRacional._multiplicacion_por_duplicacion(self, entero)

    __rmul__ = __mul__
