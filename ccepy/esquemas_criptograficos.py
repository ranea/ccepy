"""Esquemas criptográficos con curvas elípticas.

Este módulo permite trabajar con protocolos criptográficos asimétricos,
como el esquema Diffie-Hellman conocido como ECDH o el algoritmo
de firmas digitales conocido como ECDSA.

Para utilizar las funciones y las clases de este módulo, debe importarlo
previamente: ::

    # reemplace ... por la función/clase que desea utilizar
    from ccepy.esquemas_criptograficos import ...

Para trabajar con un protocolo, basta iniciar las entidades involucradas
y llamar a los métodos correspondientes. Veamos un ejemplo con EDCH:

    >>> # definimos los parámetros
    >>> E = curva_eliptica_sobre_Fq(a=324, b=1287, p=3851)  # y^2 = x^3 + a x + b sobre Fp
    >>> generador = E(920, 303)
    >>> orden_generador = 8
    >>> # definimos los participantes
    >>> alicia = ECDH(E, generador, orden_generador)
    >>> bob = ECDH(E, generador, orden_generador)
    >>> # la pareja de llaves se genera automáticamente
    >>> alicia.llave_publica
    (2373,2607)
    >>> alicia.llave_privada
    2
    >>> alicia.calcula_secreto_compartido(bob.llave_publica)
    1136
    >>> bob.calcula_secreto_compartido(alicia.llave_publica)
    1136
"""
import random
import hashlib

from ccepy.curvas_elipticas import curva_eliptica_sobre_Fq
from ccepy.aritmetica_elemental import Zp, alg_euclides


class ECDH(object):
    """Representa un participante del protocolo ECDH.

    El protocolo Diffie-Hellman para curvas elípticas o ECDH es un protocolo
    de establecimiento de llaves, esto es, permite a dos entidades establecer
    un secreto compartido sobre un canal inseguro. Para ello se cálcula
    un punto secreto de una curva elíptica que solo es conocido por
    las entidades involucradas.

    Veamos un ejemplo de la creación de un participante:

        >>> # definimos los parámetros
        >>> E = curva_eliptica_sobre_Fq(a=324, b=1287, p=3851)
        >>> generador = E(920, 303)
        >>> orden_generador = 8
        >>> # definimos el participante
        >>> alicia = ECDH(E, generador, orden_generador)
        >>> # la pareja de llaves se genera automáticamente
        >>> alicia.llave_publica
        (2373,2607)
        >>> alicia.llave_privada
        2

    Una vez se tienen dos participantes, se puede establecer un
    secreto compartido llamando a :meth:`calcula_secreto_compartido`:

        >>> bob = ECDH(E, generador, orden_generador)
        >>> alicia.calcula_secreto_compartido(bob.llave_publica)
        1136
        >>> bob.calcula_secreto_compartido(alicia.llave_publica)
        1136

    Para la inicialización de un participante es necesario pasarle
    como primer argumento una curva elíptica. Para ello puede
    usar las funciones del tipo ``curva_eliptica_sobre_*`` del módulo
    :mod:`curvas_elipticas`.

    Además, en la inicialización se genera un par de llaves pública y privada
    de forma aleatoria. Si quiere utilizar un par concreto, basta
    machacar el valor de :attr:`llave_privada` y :attr:`llave_publica`
    tras la inicialización. La única restricción es que la
    llave pública sea la multiplicación escalar de la llave privada
    por el generador:

        >>> eva = ECDH(E, generador, orden_generador)
        >>> eva.llave_privada = 3
        >>> eva.llave_publica = eva.generador * eva.llave_privada

    Args:
        curva_eliptica: el constructor de puntos de una curva elíptica.
        generador: un punto de la curva elíptica.
        orden(int): el orden del generador.

    Attributes:
        curva_eliptica: el constructor de puntos de una curva elíptica.
        generador: un punto de la curva elíptica.
        orden(int): el orden del generador.
        llave_privada(int): un entero que hace de llave privada del participante.
        llave_publica: un punto de la curva elíptica que hace de llave pública
            del participante.

    """
    def __init__(self, curva_eliptica, generador, orden):
        self.curva_eliptica = curva_eliptica
        self.generador = generador
        self.orden = orden

        # generamos las llaves
        self.llave_privada = random.randrange(1, self.orden)
        self.llave_publica = self.llave_privada * self.generador

    def calcula_secreto_compartido(self, otra_llave_publica):
        """Devuelve el punto de la curva elíptica que hace de
        secreto compartido.

        Args:
            otra_llave_publica: la llave pública del otro participante.

        Returns:
            el secreto compartido.
        """
        return (self.llave_privada * otra_llave_publica).x


class ECDSA(object):
    """Representa un participante del protocolo ECDSA.

    El algoritmo de firmas digitales para curvas elípticas o ECDSA
    es un esquema de firmado análogo al esquema DSA pero utilizando
    la multiplicación escalar de una curva elíptica en vez de la
    exponenciación modular.

    Veamos un ejemplo de la creación del participante que hará de firmante:

        >>> # definimos los parámetros
        >>> E = curva_eliptica_sobre_Fq(a=5, b=7, p=113)
        >>> generador = E(16, 51)
        >>> orden_generador = 127
        >>> # definimos el participante
        >>> alicia = ECDH(E, generador, orden_generador)
        >>> # la pareja de llaves se genera automáticamente
        >>> alicia.llave_publica
        (14,83)
        >>> alicia.llave_privada
        32

    Para firmar, basta llamar al método :meth:`firma` pasándole el
    mensaje que se desea firmar:

        >>> r, s = alicia.firma(mensaje)
        >>> r, s
        13, 73

    Para que otra entidad compruebe la firma, basta llamar a
    :meth:`verifica` pasándole como argumentos el mensaje,
    la firma y la llave pública del firmante:

        >>> bob = ECDSA(E, generador, orden)
        >>> bob.verifica(mensaje, r, s, alicia.llave_publica)
        True

    Para la inicialización de un participante es necesario pasarle
    como primer argumento una curva elíptica sobre un cuerpo finito
    de orden un primo (no una potencia de un primo). Para ello puede
    usar la función :func:`.curva_eliptica_sobre_Fq` (del módulo
    :mod:`curvas_elipticas`) paśandole como argumento ``n=1``

    Además, en la inicialización se genera un par de llaves pública y privada
    de forma aleatoria. Si quiere utilizar un par concreto, basta
    machacar el valor de :attr:`llave_privada` y :attr:`llave_publica`
    tras la inicialización. La única restricción es que la
    llave pública sea la multiplicación escalar de la llave privada
    por el generador:

        >>> eva = ECDH(E, generador, orden_generador)
        >>> eva.llave_privada = 3
        >>> eva.llave_publica = eva.generador * eva.llave_privada

    Args:
        curva_eliptica: el constructor de puntos de una curva elíptica.
        generador: un punto de la curva elíptica.
        orden(int): el orden del generador.

    Attributes:
        curva_eliptica: el constructor de puntos de una curva elíptica.
        generador: un punto de la curva elíptica.
        orden(int): el orden del generador.
        llave_privada(int): un entero que hace de llave privada del participante.
        llave_publica: un punto de la curva elíptica que hace de llave pública
            del participante.

    """
    def __init__(self, curva_eliptica, generador, orden):
        if hasattr(curva_eliptica.Fq, 'n'):
            raise ValueError("El cardinal del cuerpo finito debe ser un número primo.")

        self.curva_eliptica = curva_eliptica
        self.generador = generador
        self.orden = orden  # debe ser primo

        # generamos las llaves
        self.llave_privada = random.randrange(1, self.orden)
        self.llave_publica = self.llave_privada * self.generador

    def firma(self, mensaje):
        """Firma el mensaje utilizando la llave pública del participante.

        No se firma todo el mensaje, sino que previo al firmado
        se le aplica la función hash ``SHA-1`` (del módulo :mod:`hashlib`)
        al mensaje y se firma sobre dicho hash.

        Args:
            mensaje(str): el mensaje que se desea firmar.

        Returns:
            Tuple[int]: el par ``(r, s)`` que forma la firma del mensaje.
        """
        # renombramos las variables
        P = self.generador
        n = self.orden
        d = self.llave_privada
        m = mensaje

        while True:
            k = random.randrange(1, self.orden - 1)
            kP = k * P

            Zn = Zp(n)
            r = Zn(kP.x)
            if r == 0:
                continue

            hash_mensaje = hashlib.sha1(bytes(m, 'utf-8')).digest()
            e = int.from_bytes(hash_mensaje[:n.bit_length()], byteorder='big')

            inverso_k = Zn(k).inverso()
            s = inverso_k * (e + d * r)
            if s == 0:
                continue

            return int(r), int(s)

    def verifica(self, mensaje, r, s, llave_publica_firmante):
        """Comprueba la firma de un mensaje.

        No se firma todo el mensaje, sino que previo al firmado
        se le aplica la función hash ``SHA-1`` (del módulo :mod:`hashlib`)
        al mensaje y se firma sobre dicho hash.

        Los parámetros ``(r, s)`` es el par que devuelve el método
        :meth:`firma`.

        Args:
            mensaje(str): el mensaje que se desea comprobar su firma.
            r(int): la primera componente de la firma.
            s(int): la segunda componente de la firma.
            llave_publica_firmante: la llave pública del firmado.

        Returns:
            bool: verdadero o falso.
        """
        # renombramos las variables
        P = self.generador
        n = self.orden
        Q = llave_publica_firmante
        m = mensaje

        if not (1 <= r <= n - 1 and 1 <= s <= n - 1):
            return False

        hash_mensaje = hashlib.sha1(bytes(m, 'utf-8')).digest()
        e = int.from_bytes(hash_mensaje[:n.bit_length()], byteorder='big')

        Zn = Zp(n)
        w = Zn(s).inverso()
        u1 = int(e * w)
        u2 = int(r * w)

        X = u1 * P + u2 * Q

        if X.es_elemento_neutro():
            return False

        v = Zn(X.x)
        if v == r:
            return True
        else:
            return False
