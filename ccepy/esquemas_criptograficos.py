import random
import hashlib
# from math import gcd

# from ccepy.cuerpos_finitos import Fq
from ccepy.curvas_elipticas import curva_eliptica_sobre_Fq
from ccepy.aritmetica_elemental import Zp, alg_euclides


# añadir ref: sec 1 3.3
class ECDH(object):
    """docstring for ECDH"""
    def __init__(self, curva_eliptica, generador, orden):
        self.curva_eliptica = curva_eliptica
        self.generador = generador
        self.orden = orden

        # generamos las llaves
        self.llave_privada = random.randrange(1, self.orden)
        self.llave_publica = self.llave_privada * self.generador

    def calcula_secreto_compartido(self, otra_llave_publica):
        self.secreto_compartido = (self.llave_privada * otra_llave_publica).x


class ECDSA(object):
    """docstring for ECDH"""
    def __init__(self, curva_eliptica, generador, orden):
        if hasattr(curva_eliptica.Fq, 'n'):
            raise ValueError("El cardinal del cuerpo finito debe ser un número primo.")

        self.curva_eliptica = curva_eliptica
        self.generador = generador
        self.orden = orden  # debe ser primo

        # generamos las llaves
        self.llave_privada = random.randrange(1, self.orden)
        self.llave_publica = self.llave_privada * self.generador

    # TODO: ref 4.29 guide
    def firma(self, mensaje):
        """..."""
        # renombramos las variables para usar la notación de la referencia
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

    # TODO: ref 4.30 guide
    def verifica(self, mensaje, r, s, llave_publica_firmador):
        # renombramos las variables para usar la notación de la referencia
        P = self.generador
        n = self.orden
        Q = llave_publica_firmador
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
