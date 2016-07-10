"""Aritmética elemental con enteros y polinomios.

Este módulo permite operar con enteros módulo un primo p y polinomios
cuyos coeficientes sean enteros módulo un primo p.

Para utilizar las funciones y las clases de este módulo, debe importarlo
previamente: ::

    # reemplace ... por la función/clase que desea utilizar
    from ccepy.aritmetica_elemental import ...

Para operar con enteros módulo un primo p, use la función :func:`Zp` y
los operadores aritméticos habituales.

    >>> Z7 = Zp(7)
    >>> n, m = Z7(2), Z7(6)
    >>> n
    2
    >>> m
    6
    >>> n + m
    1
    >>> n * m
    5
    >>> m ** (-1)
    6

Para operar con polinomios con coeficientes enteros módulo un primo p,
use la función :func:`PolinomioZp` y los operadores aritméticos
habituales.

    >>> f = PolinomioZp([0, 0, 1], p=2)
    >>> f
    X^2
    >>> g = PolinomioZp([1, 1], p=2)
    >>> g
    X + 1
    >>> f + g
    X^2 + X + 1
    >>> f * g
    X^3 + X^2
    >>> f ** 3
    X^6
"""
from itertools import zip_longest
import functools
import math
import copy
import random


def alg_euclides(a, b):
    """Calcula el algoritmo extendido de Euclides para enteros.

    Esto es, los (x, y, d) tal que :math:`a x + b y = d`, siendo d el máximo
    común divisor de (a, b).

        >>> alg_euclides(54, 24)
        (1, -2, 6)

    Args:
        a (int): un número positivo.
        b (int): otro número positivo.

    Returns:
        List[int]: la lista [x, y, d].
    """
    if b > a:
        y, x, d = alg_euclides(b, a)  # Intercambiamos x, y
    else:
        if b == 0:
            d, x, y = a, 1, 0
        else:
            x2, x1, y2, y1 = 1, 0, 0, 1
            while b > 0:
                q, r = divmod(a, b)
                x, y = x2 - q * x1, y2 - q * y1
                a, b = b, r
                x2, x1, y2, y1 = x1, x, y1, y
            d, x, y = a, x2, y2
    return x, y, d


def alg_euclides_polinomios(g, h, p):
    """Calcula el algoritmo extendido de Euclides para polinomios.

    Esto es, los (s, t, d) tal que :math:`s g + t h = d`,
    siendo d el máximo común divisor mónico de (s, g).

        >>> f = PolinomioZp([0, 0, 0, 1], p=2)
        >>> f
        X^3
        >>> g = PolinomioZp([1, 0, 1, 1], p=2)
        >>> g
        X^3 + X^2 + 1
        >>> alg_euclides_polinomios(f, g, p=2)
        (X^2 + X + 1, X^2 + 1, 1)

    Args:
        g (PolinomioZp): un polinomio no nulo con coeficientes enteros
            módulo p.
        h (PolinomioZp): otro polinomio con coeficientes enteros módulo p.
        p (int): el primo p.

    Returns:
        List[PolinomioZp]: la lista [s, t, d].
    """
    gg = PolinomioZp(g.coeficientes, p)
    hh = PolinomioZp(h.coeficientes, p)
    cero = PolinomioZp([0], p)
    uno = PolinomioZp([1], p)
    if hh == cero:
        d, s, t = gg, uno, cero
    else:
        s2, s1, t2, t1 = uno, cero, cero, uno
        while hh != cero:
            q, r = divmod(gg, hh)
            s, t = s2 - q * s1, t2 - q * t1
            gg, hh = hh, r
            s2, s1, t2, t1 = s1, s, t1, t
        d, s, t = gg, s2, t2
    if d.coeficiente_lider() != 1:
        Z_p = Zp(p)
        inverso = Z_p(d.coeficiente_lider()).inverso()
        k = PolinomioZp([inverso], p)
        s, t, d = s * k, t * k, d * k
    return s, t, d


@functools.lru_cache()
def Zp(p):
    """Devuelve el constructor de enteros módulo un primo p.

        >>> Z2 = Zp(2)
        >>> Z2
        <class 'ccepy.aritmetica_elemental.Zp.<locals>.EnteroModuloP'>
        >>> Z2(11)  # 11 mod 2
        1

    Args:
        p (int): un número primo.

    Return:
        EnteroModuloP: la clase que representa los enteros módulo un primo p.
    """
    # Copiar la clase fuera de la función para que aparezca en la documentación
    class EnteroModuloP(int):
        """Representa un entero módulo un primo p.

            >>> Z7 = Zp(7)
            >>> n, m = Z7(2), Z7(6)
            >>> type(n)
            <class 'ccepy.aritmetica_elemental.Zp.<locals>.EnteroModuloP'>
            >>> n
            2
            >>> m
            6
            >>> n + m
            1
            >>> n * m
            5
            >>> m ** (-1)
            6

        Soporta los operadores ``+``, ``-``, ``*``, ``/`` y ``**`` con su
        significado habitual.

        Los operandos pueden ser ambos de tipo :class:`EnteroModuloP` o bien uno
        de tipo :class:`EnteroModuloP` y otro de tipo :py:class:`int`. En ambos
        casos el resultado será de tipo :class:`EnteroModuloP`.

        Args:
            entero (int): el valor del entero.

        Attributes:
            p (int): el primo p (*atributo de clase*).
        """
        p = None

        @classmethod
        def cero(cls):
            """Devuelve el cero.

            Return:
                EnteroModuloP: el cero.
            """
            return EnteroModuloP(0)

        @classmethod
        def uno(cls):
            """Devuelve el uno.

            Return:
                EnteroModuloP: el uno.
            """
            return EnteroModuloP(1)

        def __new__(cls, entero):
            if EnteroModuloP.p is None:
                raise RuntimeError("Instancie usando la función Zp()")
            return super().__new__(cls, entero % EnteroModuloP.p)

        # descomentar este método solo para la documentación
        # def __init__(self, entero):
        #    pass

        def __eq__(self, m):
            return super().__eq__(EnteroModuloP(m))

        def __ne__(self, m):
            return not self.__eq__(m)

        def __add__(self, m):
            return EnteroModuloP(super().__add__(m))

        __radd__ = __add__

        def __neg__(self):
            return EnteroModuloP(super().__neg__())

        def __sub__(self, m):
            return EnteroModuloP(self + (-m))

        def __rsub__(self, m):
            return -self.__sub__(m)

        def __mul__(self, m):
            return EnteroModuloP(super().__mul__(m))

        __rmul__ = __mul__

        def __pow__(self, m):
            if m < 0:
                inverso = self.inverso()
                return inverso ** (-m)
            else:
                return EnteroModuloP(super().__pow__(m, EnteroModuloP.p))

        def inverso(self):
            """Devuelve el inverso módulo p.

                >>> Z7 = Zp(7)
                >>> Z7(6).inverso()
                6

            Return:
                EnteroModuloP: el inverso.
            """
            if self == 0:
                raise ZeroDivisionError

            x, y, d = alg_euclides(self, EnteroModuloP.p)
            return EnteroModuloP(x)

        def __truediv__(self, m):
            return self * EnteroModuloP(m).inverso()

        def __rtruediv__(self, m):
            return EnteroModuloP(m).__truediv__(self)

        # Necesario para @lru_cache
        def __hash__(self):
            return super().__hash__()

    EnteroModuloP.p = p
    EnteroModuloP.__name__ = "Z{0}".format(p)
    return EnteroModuloP


class PolinomioZp:
    """Representa un polinomio con coeficientes enteros módulo un primo p.

        >>> f = PolinomioZp([0, 0, 1], p=2)
        >>> f
        X^2
        >>> g = PolinomioZp([1, 1], p=2)
        >>> g
        X + 1
        >>> f + g
        X^2 + X + 1
        >>> f * g
        X^3 + X^2
        >>> f ** 3
        X^6

    Soporta los operadores ``+``, ``-``, ``*``, ``/``, ``%`` y ``**`` con su
    significado habitual.

    Los operandos pueden ser ambos de tipo :class:`PolinomioZp` o bien uno de
    tipo :class:`PolinomioZp` y otro de tipo :py:class:`int`. En ambos casos el
    resultado será de tipo :class:`PolinomioZp`.

    Args:
        coeficientes (List[int]): los coeficientes del polinomio ordenados
            de forma ascendente, esto es, el primero el término constante y
            el último el coeficiente líder.
        p (int): el primo p.
    """
    def __init__(self, coeficientes, p):
        # Queremos que el último coeficiente no sea nulo
        # (excepto si es el polinomio cero)
        Z_p = Zp(p)
        ultimo_coef = None
        for indice, coef in enumerate(reversed(coeficientes)):
            if Z_p(coef) != 0:
                ultimo_coef = len(coeficientes) - indice - 1
                break
        else:
            ultimo_coef = 0

        self._coeficientes = [Z_p(c) for c in coeficientes[:ultimo_coef + 1]]

    @property
    def coeficientes(self):
        """List[EnteroModuloP]: los coeficientes del polinomio ordenados de
        forma ascendente, esto es, el primero es el término constante y el
        último el coeficiente líder. Es un atributo de solo lectura.
        """
        return self._coeficientes

    def primo(self):
        """Devuelve el primo p."""
        return self._coeficientes[0].p

    @classmethod
    def monomio(cls, coef, grado, p):
        """Devuelve el monomio con coeficiente *coef* y de grado *grado*.

            >>> PolinomioZp.monomio(-1, 7, 2)
            X^7

        Args:
            coef (int): el coeficiente líder del monomio.
            grado (int): el exponente del monomio.

        Returns:
            PolinomioZp: el monomio con dicho coeficiente y grado.

        """
        return cls([0 for i in range(grado)] + [coef], p)

    def grado(self):
        """Devuelve el grado del polinomio.

            >>> f = PolinomioZp([1, 0, 0, 1], p=2)
            >>> f
            X^3 + 1
            >>> f.grado()
            3

        El grado puede ser:

            * \- :py:data:`math.inf` : si el polinomio es el polinomio cero.
            * ``n`` : si el término lider tiene exponente n.

        Returns:
            int: el grado del polinomio.
        """
        if self == PolinomioZp([0], self.primo()):
            return -math.inf
        else:
            return len(self.coeficientes) - 1

    def coeficiente_lider(self):
        """Devuelve el coeficiente asociado al término de mayor exponente.

            >>> f = PolinomioZp([2, 0, 0, 1], p=3)
            >>> f
            X^3 + 2
            >>> f.coeficiente_lider()
            1

        Returns:
            EnteroModuloP: el coeficiente asociado al mayor exponente.
        """
        return self.coeficientes[-1]

    def es_irreducible(self):
        """Comprueba si el polinomio es irreducible.

            >>> f = PolinomioZp([1, 0, 1, 1], p=2)
            >>> f
            X^3 + X^2 + 1
            >>> f.es_irreducible()
            True

        Returns:
            bool: verdadero o falso.
        """
        p = self.primo()
        f = PolinomioZp(self.coeficientes, p)
        if f.coeficiente_lider() != 1:
            f = f / PolinomioZp([f.coeficiente_lider()], p)  # lo hacemos mónico
        m = f.grado()

        x = PolinomioZp([0, 1], p)
        u = copy.deepcopy(x)
        for i in range(1, m // 2 + 1):
            u = (u ** p) % f
            _, _, d = alg_euclides_polinomios(f, u - x, p)
            if d != PolinomioZp([1], p):
                return False
        return True

    @classmethod
    def genera_irreducible(cls, grado, p):
        """Devuelve un polinomio irreducible de dicho grado con coeficientes
        módulo p.

            >>> f = PolinomioZp.genera_irreducible(3, 2)
            >>> f.es_irreducible()
            True

        Returns:
            PolinomioZp: el polinomio irreducible.
        """
        Z_p = Zp(p)
        while True:
            a_0 = Z_p(1 + random.randrange(p - 1))  # lo queremos != 0
            a_m = Z_p(1)
            coeficientes = [a_0]
            coeficientes += [random.randrange(p) for i in range(1, grado)]
            coeficientes += [a_m]
            f = PolinomioZp(coeficientes, p)
            if f.es_irreducible():
                return f

    def __eq__(self, q):
        if isinstance(q, PolinomioZp):
            return self.coeficientes == q.coeficientes
        else:
            # Si el polinomio es una constante, hacemos
            # la comparación con el coeficiente
            if self.grado() < 1:
                return self.coeficientes[0] == q
            else:
                return False

    def __ne__(self, q):
        return not self.__eq__(q)

    def __add__(self, q):
        if isinstance(q, PolinomioZp):
            return PolinomioZp([a + b for a, b in zip_longest(self.coeficientes,
                                                q.coeficientes,
                                                fillvalue=0)], self.primo())
        else:
            coeficientes = [self.coeficientes[0] + q] + self.coeficientes[1:]
            return PolinomioZp(coeficientes, self.primo())

    __radd__ = __add__

    def __neg__(self):
        return PolinomioZp([-a for a in self.coeficientes], self.primo())

    def __sub__(self, q):
        return self + (-q)

    def __rsub__(self, q):
        return -self.__sub__(q)

    def __mul__(self, q):
        if isinstance(q, PolinomioZp):
            cero = PolinomioZp([0], self.primo())
            if self == cero or q == cero:
                return cero

            maximo_grado = len(self.coeficientes) + len(q.coeficientes)
            multiplicacion = [0 for _ in range(maximo_grado)]
            for i, a in enumerate(self.coeficientes):
                for j, b in enumerate(q.coeficientes):
                    multiplicacion[i + j] += a * b

            return PolinomioZp(multiplicacion, self.primo())
        else:
            return PolinomioZp([a * q for a in self.coeficientes], self.primo())

    __rmul__ = __mul__

    def __pow__(self, n):
        potencia = PolinomioZp([1], self.primo())
        for _ in range(n):
            potencia *= self
        return potencia

    def __divmod__(self, q):
        p = self.primo()
        cero = PolinomioZp([0], self.primo())
        cociente, divisor, resto = cero, copy.deepcopy(q), copy.deepcopy(self)
        while resto != cero and resto.grado() >= divisor.grado():
            monomio_grado = resto.grado() - divisor.grado()
            monomio_cl = resto.coeficiente_lider() / divisor.coeficiente_lider()
            monomio_cociente = PolinomioZp.monomio(monomio_cl, monomio_grado, p)

            cociente += monomio_cociente
            resto -= monomio_cociente * divisor

        return cociente, resto

    def __truediv__(self, q):
        return divmod(self, q)[0]

    def __mod__(self, q):
        if self.grado() < q.grado():
            return self
        else:
            return divmod(self, q)[1]

    def __str__(self):
        if self == PolinomioZp([0], self.primo()):
            return str(0)
        else:
            monomios = []
            # Se imprime los monomios en orden descedente respecto al grado
            for indice, coef in enumerate(reversed(self.coeficientes)):
                if coef != 0:
                    exponente = len(self.coeficientes) - indice - 1
                    # La siguiente casuística es escribir X
                    # en lugar de 1*X^1 y casos similares
                    if exponente == 0:
                        monomios.append(str(coef))
                    elif exponente == 1:
                        if coef == 1:
                            monomios.append("X")
                        else:
                            monomios.append("{0}*X".format(coef))
                    else:
                        if coef == 1:
                            monomios.append("X^{0}".format(exponente))
                        else:
                            monomios.append("{0}*X^{1}".format(coef, exponente))
            return ' + '.join(monomios)

    __repr__ = __str__

    # Necesario para @functools.lru_cache de Fq()
    def __hash__(self):
        return hash(tuple(self.coeficientes))
