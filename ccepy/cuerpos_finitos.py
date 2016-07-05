"""Cuerpos finitos.

Este módulo permite operar con elementos de un cuerpos finito de q elementos,
donde q será la potencia de un primo.

Para utilizar las funciones y las clases de este módulo, debe importarlo
previamente: ::

    # reemplace ... por la función/clase que desea utilizar
    from ccepy.cuerpos_finitos import ...

Para operar con elementos de un cuerpo finito con q elementos, use la función
:func:`Fq` y los operadores aritméticos habituales.

    >>> F16 = Fq(2, 4)
    >>> F16([1, 1, 0, 1]) + F16([1, 0, 0, 1])
    {[0, 1, 0, 0]; 16}
    >>> F16([1, 0, 1, 1]) * F16([1, 0, 0, 1])
    {[1, 0, 0, 0]; 16}
    >>> F16([1, 1, 0, 1]) ** (-1)
    {[0, 1, 0, 1]; 16}
    >>> F16([0, 0, 1, 0]) ** 2
    {[1, 0, 0, 1]; 16}

Se está utilizando la representación polinomial para los elementos de un
cuerpo finito. En particular, para la creación y la representación de
un elemento se utilizan los coeficientes del elemento visto como polinomio.
"""
import copy
import functools

from ccepy.aritmetica_elemental import Zp, PolinomioZp, alg_euclides_polinomios


@functools.lru_cache()
def Fq(p, n=1, pol_irreducible=None):
    """Devuelve el constructor de elementos del cuerpo finito con p**n elementos.

        >>> F16 = Fq(2, 4)
        >>> F16
        <class 'ccepy.cuerpos_finitos.Fq.<locals>.ElementoFq'>
        >>> F16([0, 0, 0, 0, 1])
        {[1, 0, 0, 1]; 16}

    Se puede especificar el polinomio con el cual se hace módulo:

        >>> pol_irreducible = PolinomioZp([1, 1, 1, 0, 1], p=2)
        >>> pol_irreducible
        X^4 + X^2 + X + 1
        >>> F16 = Fq(2, 4, pol_irreducible)
        >>> F16([0, 0, 0, 0, 1])  # X^4 mod X^4 + X^2 + X + 1
        {[1, 1, 1, 0]; 16}


    Args:
        p (int): un número primo.
        n (Optional[int]): un número natural.
        pol_irreducible (Optional[PolinomioZp]): un polinomio de grado
            *n* irreducible.

    Return:
        Si n es uno, devuelve :class:`.EnteroModuloP`.

        Si n es mayor que uno, devuelve :class:`ElementoFq`.
    """
    if n == 1:
        return Zp(p)

    # Copiar la clase fuera de la función para que aparezca en la documentación
    class ElementoFq(PolinomioZp):
        """Representa un elemento del cuerpo finito con q elementos.

            >>> F16 = Fq(2, 4)
            >>> p, q = F16([1, 1, 0, 1]), F16([1, 0, 0, 1])
            >>> type(p)
            <class 'ccepy.cuerpos_finitos.Fq.<locals>.ElementoFq'>
            >>> p
            {[1, 1, 0, 1]; 16}
            >>> q
            {[1, 0, 0, 1]; 16}
            >>> p + q
            {[0, 1, 0, 0]; 16}
            >>> p * q
            {[1, 0, 1, 0]; 16}
            >>> p ** (-1)
            {[0, 1, 0, 1]; 16}

        Soporta los operadores ``+``, ``-``, ``*``, ``/`` y ``**`` con su
        significado habitual.

        Args:
            coeficientes (List[int]): los coeficientes del elemento visto
                como polinomio.

        Attributes:
            Zp (EnteroModuloP): el constructor de enteros módulo un primo p.
                (*atributo de clase*)
            p (int): el primo p. (*atributo de clase*)
            n (int): el grado del irreducible *pol_irreducible*.
                (*atributo de clase*)
            q (int): el número de elementos del cuerpo finito.
                (*atributo de clase*)
        """
        Zp = None
        p = None
        n = None
        pol_irreducible = None

        @classmethod
        def cero(cls):
            """Devuelve el cero del cuerpo finito.

            Return:
                ElementoFq: el cero.
            """
            return ElementoFq([0])

        @classmethod
        def uno(cls):
            """Devuelve el uno del cuerpo finito.

            Return:
                ElementoFq: el uno.
            """
            return ElementoFq([1])

        def __init__(self, coeficientes):
            if isinstance(coeficientes, int):
                pol = PolinomioZp([coeficientes], ElementoFq.p)
            elif isinstance(coeficientes, list):
                pol = PolinomioZp(coeficientes, ElementoFq.p)
            else:
                # PolinomioZp o ElementoFq
                pol = coeficientes

            coeficientes_nuevos = (pol % ElementoFq.pol_irreducible).coeficientes
            super().__init__(coeficientes_nuevos, ElementoFq.p)

        def __eq__(self, alfa):
            return super().__eq__(ElementoFq(alfa))

        def __ne__(self, alfa):
            return not self.__eq__(alfa)

        def __add__(self, alfa):
            return ElementoFq(super().__add__(alfa))

        __radd__ = __add__

        def __neg__(self):
            return ElementoFq(super().__neg__())

        def __sub__(self, alfa):
            return ElementoFq(self + (-alfa))

        def __rsub__(self, alfa):
            return -self.__sub__(alfa)

        def __mul__(self, alfa):
            return ElementoFq(super().__mul__(alfa))

        __rmul__ = __mul__

        @classmethod
        def _exponenciacion_binaria(cls, g, k):
            """Calcula la potencia k-ésima del elemento g eficientemente.

            Args:
                g (ElementoFq): un elemento del cuerpo finito.
                k (int): un exponente natural entre el 0 y q-2 (ambos inclusive).

            Return:
                ElementoFq: la potencia k-ésima del elemento g.
            """
            rep_binaria_k = "".join(reversed(bin(k)[2:]))
            t = len(rep_binaria_k) - 1

            s = ElementoFq([1])
            if k == 0:
                return s
            G = copy.deepcopy(g)
            if rep_binaria_k[0] == "1":
                s = copy.deepcopy(G)
            for i in range(1, t + 1):
                G = G * G
                if rep_binaria_k[i] == "1":
                    s = G * s
            return s

        def __pow__(self, k):
            if self == ElementoFq.cero():
                return self
            if self == ElementoFq.uno() or k == 0:
                return ElementoFq.uno()

            q = ElementoFq.q
            if k < 0:
                inverso = self.inverso()
                return ElementoFq._exponenciacion_binaria(inverso, -k % (q - 1))
            else:
                return ElementoFq._exponenciacion_binaria(self, k % (q - 1))

        def inverso(self):
            """Devuelve el inverso del elemento del cuerpo finito.

                >>> F16 = Fq(2, 4)
                >>> F16([1, 1, 0, 1]).inverso()
                {[0, 1, 0, 1]; 16}

            Returns:
                ElementoFq: el inverso.
            """
            if self == ElementoFq.cero():
                raise ZeroDivisionError

            s, t, d = alg_euclides_polinomios(self, ElementoFq.pol_irreducible, ElementoFq.p)
            return ElementoFq(s)

        def __truediv__(self, alfa):
            return self * ElementoFq(alfa).inverso()

        def __rtruediv__(self, alfa):
            return ElementoFq(alfa).__truediv__(self)

        def __str__(self):
            tope = ElementoFq.n - len(self.coeficientes)
            coeficientes = self.coeficientes + [0 for _ in range(0, tope)]
            return "{{{0}; {1}}}".format(coeficientes, ElementoFq.q)

        __repr__ = __str__

    ElementoFq.Zp = Zp(p)
    ElementoFq.p = p
    ElementoFq.n = n
    ElementoFq.q = p ** n
    if pol_irreducible is None:
        pol_irreducible = PolinomioZp.genera_irreducible(grado=n, p=p)
    ElementoFq.pol_irreducible = pol_irreducible
    ElementoFq.__name__ = "F{0}".format(p**n)
    return ElementoFq


# TODO: remover clase tras realizar documentación
class ElementoFq(PolinomioZp):
    """Representa un elemento del cuerpo finito con q elementos.

        >>> F16 = Fq(2, 4)
        >>> p, q = F16([1, 1, 0, 1]), F16([1, 0, 0, 1])
        >>> type(p)
        <class 'ccepy.cuerpos_finitos.Fq.<locals>.ElementoFq'>
        >>> p
        {[1, 1, 0, 1]; 16}
        >>> q
        {[1, 0, 0, 1]; 16}
        >>> p + q
        {[0, 1, 0, 0]; 16}
        >>> p * q
        {[1, 0, 1, 0]; 16}
        >>> p ** (-1)
        {[0, 1, 0, 1]; 16}

    Soporta los operadores ``+``, ``-``, ``*``, ``/`` y ``**`` con su
    significado habitual.

    Args:
        coeficientes (List[int]): los coeficientes del elemento visto
            como polinomio.

    Attributes:
        Zp (EnteroModuloP): el constructor de enteros módulo un primo p.
            (*atributo de clase*)
        p (int): el primo p. (*atributo de clase*)
        n (int): el grado del irreducible *pol_irreducible*.
            (*atributo de clase*)
        q (int): el número de elementos del cuerpo finito.
            (*atributo de clase*)
    """
    Zp = None
    p = None
    n = None
    pol_irreducible = None

    @classmethod
    def cero(cls):
        """Devuelve el cero del cuerpo finito.

        Return:
            ElementoFq: el cero.
        """
        return ElementoFq([0])

    @classmethod
    def uno(cls):
        """Devuelve el uno del cuerpo finito.

        Return:
            ElementoFq: el uno.
        """
        return ElementoFq([1])

    def __init__(self, coeficientes):
        if isinstance(coeficientes, int):
            pol = PolinomioZp([coeficientes], ElementoFq.p)
        elif isinstance(coeficientes, list):
            pol = PolinomioZp(coeficientes, ElementoFq.p)
        else:
            # PolinomioZp o ElementoFq
            pol = coeficientes

        coeficientes_nuevos = (pol % ElementoFq.pol_irreducible).coeficientes
        super().__init__(coeficientes_nuevos, ElementoFq.p)

    def __eq__(self, alfa):
        return super().__eq__(ElementoFq(alfa))

    def __ne__(self, alfa):
        return not self.__eq__(alfa)

    def __add__(self, alfa):
        return ElementoFq(super().__add__(alfa))

    __radd__ = __add__

    def __neg__(self):
        return ElementoFq(super().__neg__())

    def __sub__(self, alfa):
        return ElementoFq(self + (-alfa))

    def __rsub__(self, alfa):
        return -self.__sub__(alfa)

    def __mul__(self, alfa):
        return ElementoFq(super().__mul__(alfa))

    __rmul__ = __mul__

    @classmethod
    def _exponenciacion_binaria(cls, g, k):
        """Calcula la potencia k-ésima del elemento g eficientemente.

        Args:
            g (ElementoFq): un elemento del cuerpo finito.
            k (int): un exponente natural entre el 0 y q-2 (ambos inclusive).

        Return:
            ElementoFq: la potencia k-ésima del elemento g.
        """
        rep_binaria_k = "".join(reversed(bin(k)[2:]))
        t = len(rep_binaria_k) - 1

        s = ElementoFq([1])
        if k == 0:
            return s
        G = copy.deepcopy(g)
        if rep_binaria_k[0] == "1":
            s = copy.deepcopy(G)
        for i in range(1, t + 1):
            G = G * G
            if rep_binaria_k[i] == "1":
                s = G * s
        return s

    def __pow__(self, k):
        if self == ElementoFq.cero():
            return self
        if self == ElementoFq.uno() or k == 0:
            return ElementoFq.uno()

        q = ElementoFq.q
        if k < 0:
            inverso = self.inverso()
            return ElementoFq._exponenciacion_binaria(inverso, -k % (q - 1))
        else:
            return ElementoFq._exponenciacion_binaria(self, k % (q - 1))

    def inverso(self):
        """Devuelve el inverso del elemento del cuerpo finito.

            >>> F16 = Fq(2, 4)
            >>> F16([1, 1, 0, 1]).inverso()
            {[0, 1, 0, 1]; 16}

        Returns:
            ElementoFq: el inverso.
        """
        if self == ElementoFq.cero():
            raise ZeroDivisionError

        s, t, d = alg_euclides_polinomios(self, ElementoFq.pol_irreducible, ElementoFq.p)
        return ElementoFq(s)

    def __truediv__(self, alfa):
        return self * ElementoFq(alfa).inverso()

    def __rtruediv__(self, alfa):
        return ElementoFq(alfa).__truediv__(self)

    def __str__(self):
        tope = ElementoFq.n - len(self.coeficientes)
        coeficientes = self.coeficientes + [0 for _ in range(0, tope)]
        return "{{{0}; {1}}}".format(coeficientes, ElementoFq.q)

    __repr__ = __str__
