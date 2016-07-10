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

from ccepy.cuerpos_finitos import Fq, PolinomioZp  # PolinomioZp para los test


class PuntoRacional(metaclass=ABCMeta):
    """Clase abstracta que representa un punto racional de una curva elíptica.
    El resto de las clases ``Punto*Racional`` heredan de esta.

    Esta clase no se puede instanciar. Sirve como punto de partida para
    crear curvas elípticas sobre nuevos cuerpos.
    """
    __slots__ = ('_x', '_y')  # para mejorar la eficiencia si hay muchos objetos

    @classmethod
    @abstractmethod
    def contiene(cls, x, y):
        """Comprueba si (x,y) esta en la curva."""
        return

    @abstractmethod
    def __init__(self, x, y):
        # Debe inicializar self._x y self._y.
        return

    def es_elemento_neutro(self):
        """Comprueba si es el elemento neutro (el punto del infinito).

        Returns:
            bool: verdadero o falso.
        """
        return self._x is None or self._y is None

    @property
    def x(self):
        """La componente x del punto si no es el elemento neutro. Es un
        atributo de solo lectura."""
        if self.es_elemento_neutro():
            raise AttributeError("El elemento neutro no tiene componente x")
        else:
            return self._x

    @property
    def y(self):
        """La componente y del punto si no es el elemento neutro. Es un
        atributo de solo lectura."""
        if self.es_elemento_neutro():
            raise AttributeError("El elemento neutro no tiene componente y")
        else:
            return self._y

    @classmethod
    def elemento_neutro(cls):
        """Devuelve el elemento neutro.

        Returns:
            El elemento neutro."""
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

    __repr__ = __str__


def curva_eliptica_sobre_Fq(a, b, p, n=1, pol_irreducible=None):
    """Devuelve el constructor de puntos de una curva elíptica sobre
    un cuerpo finito de q elementos de característica distinta de 2 y 3.

        >>> E = curva_eliptica_sobre_Fq(1, 1, 5, 2)  # y^2 = x^3 + x + 1 sobre F25
        >>> E
        <class 'ccepy.curvas_elipticas.curva_eliptica_sobre_Fq.<locals>.PuntoFqRacional'>
        >>> E(0, 1)
        ({[0, 0]; 25},{[1, 0]; 25})

    Los dos primeros argumentos (``a``, ``b``) son los coeficientes de la ecuación
    de Weierstrass simplificada: :math:`y^2 = x^3 + a x + b`. Estos valores
    pueden ser bien de tipo :py:class:`int` o bien de tipo :class:`.EnteroModuloP` o
    :class:`.ElementoFq` según sea ``n`` uno o mayor que uno respectivamente.

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
    # Copiar la clase fuera de la función para que aparezca en la documentación
    class PuntoFqRacional(PuntoRacional):
        """Representa un punto de una curva elíptica sobre un cuerpo finito de
        q elementos de característica distinta de 2 y 3.

            >>> E = curva_eliptica_sobre_Fq(1, 1, 5, 2)  # y^2 = x^3 + x + 1 sobre F25
            >>> F25 = Fq(5, 2)
            >>> P = E(F25.cero(), F25.uno())
            >>> type(P)
            <class 'ccepy.curvas_elipticas.curva_eliptica_sobre_Fq.<locals>.PuntoFqRacional'>
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

        La curva elíptica está definida por la ecuación de Weierstrass
        simplificada :math:`y^2 = x^3 + a x + b`.

        Soporta los operadores ``+``, ``-``, ``*`` con su significado habitual.

        Los parámetros ``x``, ``y`` deben ser del tipo :class:`.EnteroModuloP` o
        :class:`.ElementoFq` según el cuerpo finito tenga un número primo de
        elementos o un número potencia de un primo de elementos. Para construir
        elementos de estos tipos, utilice :func:`.Fq`.

        Args:
            x: un elemento del cuerpo finito de q elementos.
            y: un elemento del cuerpo finito de q elementos.

        Los elementos de ``coeficientes`` y ``discriminante`` serán del tipo
        :class:`.EnteroModuloP` o :class:`.ElementoFq` según el cuerpo finito
        tenga un número primo de elementos o un número potencia de un primo de
        elementos.

        Attributes:
            coeficientes (Tuple): los coeficientes (a, b) de la ecuación de Weierstrass. (atributo de clase)
            discriminante: El discriminate de la curva elíptica. (atributo de clase)
            Fq: El constructor de elementos del cuerpo finito de q elementos. (atributo de clase)
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

        @classmethod
        def _multiplicacion_por_duplicacion(cls, punto, k):
            """Realiza la multiplicación k * punto mediante el método de
            multiplicación por duplicación."""
            rep_binaria_k = "".join(bin(k)[2:])  # (k_t, k_{t-1},..., k_0)
            Q = PuntoFqRacional.elemento_neutro()
            P = punto

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


def curva_eliptica_sobre_F2m(a, b, m, pol_irreducible=None):
    """Devuelve el constructor de puntos de una curva elíptica sobre
    el cuerpo finito de 2**m elementos.

        >>> pol_irreducible = PolinomioZp([1, 1, 0, 0, 1], p=2)
        >>> F16 = Fq(2, 4, pol_irreducible)
        >>> a = F16([0, 0, 0, 1])
        >>> b = F16([1, 0, 0, 1])
        >>> E = curva_eliptica_sobre_F2m(a, b, 4, pol_irreducible)
        >>> E
        <class 'ccepy.curvas_elipticas.curva_eliptica_sobre_F2m.<locals>.PuntoF2mRacional'>
        >>> E(F16.uno(), F16.uno())
        ({[1, 0, 0, 0]; 16},{[1, 0, 0, 0]; 16})

    Los dos primeros argumentos (``a``, ``b``) son los coeficientes de la ecuación
    de Weierstrass simplificada: :math:`y^2 + x y = x^3 + a x^2 + b`. Estos valores
    pueden ser bien de tipo :py:class:`int` o bien de tipo :class:`.EnteroModuloP` o
    :class:`.ElementoFq` según sea ``n`` uno o mayor que uno respectivamente.

    Los dos últimos argumentos (``m``, ``pol_irreducible``) definen el cuerpo
    finito de 2**m elementos sobre el que se define la curva eliptipca.

    Args:
        a : el coeficiente que acompaña a x^2 en la ecuación de Weierstrass
        b : el término independiente de la ecuación de Weierstrass
        m ([int]): un número natural.
        pol_irreducible (Optional[PolinomioZp]): un polinomio de grado
            *m* irreducible.

    Return:
        PuntoF2mRacional: la clase que representa los puntos de la curva elíptica.
    """
    # Copiar la clase fuera de la función para que aparezca en la documentación
    class PuntoF2mRacional(PuntoRacional):
        """Representa un punto de una curva elíptica sobre el cuerpo finito de
        2**m elementos.

            >>> pol_irreducible = PolinomioZp([1, 1, 0, 0, 1], p=2)
            >>> F16 = Fq(2, 4, pol_irreducible)
            >>> a = F16([0, 0, 0, 1])
            >>> b = F16([1, 0, 0, 1])
            >>> E = curva_eliptica_sobre_F2m(a, b, 4, pol_irreducible)
            >>> E.coeficientes
            Coeficientes(a={[0, 0, 0, 1]; 16}, b={[1, 0, 0, 1]; 16})
            >>> P = E(F16([0, 1, 0, 0]), F16([1, 1, 1, 1]))
            >>> type(P)
            <class 'ccepy.curvas_elipticas.curva_eliptica_sobre_F2m.<locals>.PuntoF2mRacional'>
            >>> P
            ({[0, 1, 0, 0]; 16},{[1, 1, 1, 1]; 16})
            >>> Q = E(F16([0, 0, 1, 1]), F16([0, 0, 1, 1]))
            >>> Q
            ({[0, 0, 1, 1]; 16},{[0, 0, 1, 1]; 16})
            >>> P + Q
            ({[1, 0, 0, 0]; 16},{[1, 0, 0, 0]; 16})
            >>> -P
            ({[0, 1, 0, 0]; 16},{[1, 0, 1, 1]; 16})
            >>> 2 * P
            ({[1, 1, 0, 1]; 16},{[0, 1, 0, 0]; 16})

        La curva elíptica está definida por la ecuación de Weierstrass
        simplificada :math:`y^2 + x y = x^3 + a x^2 + b`.

        Soporta los operadores ``+``, ``-``, ``*`` con su significado habitual.

        Los parámetros ``x``, ``y`` deben ser del tipo :class:`.EnteroModuloP` o
        :class:`.ElementoFq` según m sea uno o mayor que uno. Para construir
        elementos de estos tipos, utilice :func:`.Fq`.

        Args:
            x: un elemento del cuerpo finito de 2**m elementos.
            y: un elemento del cuerpo finito de 2**m elementos.

        Los elementos de ``coeficientes`` y ``discriminante`` serán del tipo
        :class:`.EnteroModuloP` o :class:`.ElementoFq` según m sea uno o
        mayor que uno.

        Attributes:
            coeficientes (Tuple): los coeficientes (a, b) de la ecuación de Weierstrass. (atributo de clase)
            discriminante: El discriminate de la curva elíptica. (atributo de clase)
            Fq: El constructor de elementos del cuerpo finito de 2**m elementos. (atributo de clase)
        """
        # coeficientes (a, b) de la ecuación y^2 + x y = x^3 + a x^2 + b
        coeficientes = None
        discriminante = None
        F2m = None

        @classmethod
        def contiene(cls, x, y):
            """Comprueba si el punto (x, y) pertenece a la curva elíptica.

            Args:
                a : un elemento del cuerpo finito de 2**m elementos.
                b : un elemento del cuerpo finito de 2**m elementos.

            Returns:
                bool: verdadero o falso.
            """
            a, b = cls.coeficientes
            lado_izquierdo_ecuacion = y**2 + x * y
            lado_derecho_ecuacion = x**3 + a * x**2 + b
            return lado_izquierdo_ecuacion == lado_derecho_ecuacion

        def __init__(self, x, y):
            if x is None or y is None:
                self._x = None
                self._y = None
            else:
                self._x = PuntoF2mRacional.F2m(x)
                self._y = PuntoF2mRacional.F2m(y)
                if not PuntoF2mRacional.contiene(self._x, self._y):
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
            a = PuntoF2mRacional.coeficientes.a
            F2m = PuntoF2mRacional.F2m

            if self == other:
                if x1 == F2m(0):
                    # P = Q, P = -P, calculamos 2P
                    return PuntoF2mRacional.elemento_neutro()
                else:
                    # P = Q, P != -P, calculamos 2P
                    m = x1 + y1 / x1
                    x3 = m**2 + m + a
                    y3 = x1**2 + (m + 1) * x3
                    return PuntoF2mRacional(x3, y3)
            elif x1 == x2:
                # (y1 != y2) | P != Q, P = -Q, calculamos P - P
                return PuntoF2mRacional.elemento_neutro()
            else:
                # P != +-Q, calculamos P + Q
                m = (y1 + y2) / (x1 + x2)
                x3 = m**2 + m + x1 + x2 + a
                y3 = m * (x1 + x3) + x3 + y1
                return PuntoF2mRacional(x3, y3)

        def __neg__(self):
            if self.es_elemento_neutro():
                return self
            else:
                return PuntoF2mRacional(self.x, self.x + self.y)

        @classmethod
        def _multiplicacion_por_duplicacion(cls, punto, k):
            rep_binaria_k = "".join(bin(k)[2:])  # (k_t, k_{t-1},..., k_0)
            Q = PuntoF2mRacional.elemento_neutro()
            P = punto

            for k_i in rep_binaria_k:
                Q = Q + Q  # duplicar
                if k_i == "1":
                    Q = Q + P  # sumar

            return Q

        def __mul__(self, entero):
            if self.es_elemento_neutro():
                return self
            elif entero < 0:
                return PuntoF2mRacional._multiplicacion_por_duplicacion(-self, -entero)
            else:
                return PuntoF2mRacional._multiplicacion_por_duplicacion(self, entero)

        __rmul__ = __mul__

    F2m = Fq(2, m, pol_irreducible)
    A = F2m(a)
    B = F2m(b)
    discriminante = B
    if discriminante == F2m.cero():
        raise ValueError("El discriminant, b, no puede ser cero.")

    PuntoF2mRacional.discriminante = discriminante
    PuntoF2mRacional.coeficientes = EcuacionWeierstrass(A, B)
    PuntoF2mRacional.F2m = F2m
    return PuntoF2mRacional


def curva_eliptica_sobre_Q(a, b):
    """Devuelve el constructor de puntos de una curva elíptica sobre los
    números racionales.

        >>> E = curva_eliptica_sobre_Q(0, 4)  # y^2 = x^3 + 4 sobre Q
        >>> E
        <class 'ccepy.curvas_elipticas.curva_eliptica_sobre_Q.<locals>.PuntoQRacional'>
        >>> E(0, -2)
        (0,-2)

    Los argumentos (``a``, ``b``) son los coeficientes de la ecuación
    de Weierstrass simplificada: :math:`y^2 = x^3 + a x + b`. Estos valores
    pueden ser bien de tipo :py:class:`int` o bien de tipo :py:class:`fractions.Fraction`.

    Args:
        a : el coeficiente que acompaña a x en la ecuación de Weierstrass
        b : el término independiente de la ecuación de Weierstrass

    Return:
        PuntoQRacional: la clase que representa los puntos de la curva elíptica.
    """
    # Copiar la clase fuera de la función para que aparezca en la documentación
    class PuntoQRacional(PuntoRacional):
        """Representa un punto de una curva elíptica sobre los
        números racionales.

            >>> E = curva_eliptica_sobre_Q(-18, -72)  # y^2 = x^3 - 18 * x - 72 sobre Q
            >>> P = E(6, 6)
            >>> type(P)
            <class 'ccepy.curvas_elipticas.curva_eliptica_sobre_Q.<locals>.PuntoQRacional'>
            >>> P
            (6,6)
            >>> Q = E(Fraction(177, 4), Fraction(-2343, 8))
            >>> Q
            (177/4,-2343/8)
            >>> P + Q
            (28102/2601,4183750/132651)
            >>> -P
            (6,-6)
            >>> 4 * P
            (111795513/9759376,-1067078260371/30488290624)

        La curva elíptica está definida por la ecuación de Weierstrass
        simplificada :math:`y^2 = x^3 + a x + b`.

        Soporta los operadores ``+``, ``-``, ``*`` con su significado habitual.

        Los parámetros ``x``, ``y`` deben ser del tipo :py:class:`int` o
        :py:class:`fractions.Fraction`.

        Args:
            x: un número racional.
            y: un número racional.

        Los elementos de ``coeficientes`` y ``discriminante`` serán del tipo
        :py:class:`int` o :py:class:`fractions.Fraction` según se haya llamado
        a :func:`.curva_eliptica_sobre_Q`.

        Attributes:
            coeficientes (Tuple): los coeficientes (a, b) de la ecuación de Weierstrass. (atributo de clase)
            discriminante: El discriminate de la curva elíptica. (atributo de clase)
        """
        coeficientes = None
        discriminante = None

        @classmethod
        def contiene(cls, x, y):
            """Comprueba si el punto (x, y) pertenece a la curva elíptica.

            Args:
                a : un número racional.
                b : un número racional.

            Returns:
                bool: verdadero o falso.
            """
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
            producto = PuntoQRacional.elemento_neutro()
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
