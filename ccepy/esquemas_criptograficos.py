import random
# from math import gcd

# from ccepy.cuerpos_finitos import Fq
from ccepy.curvas_elipticas import curva_eliptica_sobre_Fq


# añadir ref: sec 1 3.3
class ECDH(object):
    """docstring for ECDH"""
    def __init__(self, curva_eliptica, generador, orden):
        self.curva_eliptica = curva_eliptica
        self.generador = generador
        self.orden = orden

        # generamos las llaves
        self.llave_publica = random.randrange(1, self.orden)
        self.llave_privada = self.llave_publica * self.generador

    def calcula_secreto_compartido(self, otra_llave_publica):
        self.secreto_compartido = self.llave_privada * otra_llave_publica


# ParametrosCurvaElipticaSobreFq = namedtuple('CurvaElípticaFq', ['nombre', 'p', 'a', 'b', 'x1', 'y1'])
#
# ejemplo_curva_eliptica = ParametrosCurvaElipticaSobreFq(
#     nombre="ejemplo1",
#     p=1287,
#     a=3851,
#     b=324,
#     x1=920,
#     y1=303,
#     orden=8
# )
#
# def fname(ce):
#     E = curva_eliptica_sobre_Fq(ce.a, ce.b, ce.p)
#     generador = E(ce.x1, ce.y1)
#     orden = ce.orden
#     return E, generador, orden

if __name__ == '__main__':
    pass



    # # Example 2: Secp256k1
    # p = int("0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f", 16)
    # a = 0
    # b = 7
    # order = int("0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
    # x_g = int("0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", 16)
    # y_g = int("0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 16)
    # ec = EllipticCurveOverZp(p, a, b)
    # alice = ECDH(ec, ec(x_g, y_g), order, 8)
    # print(alice)
    #
    # bob = ECDH(ec, ec(x_g, y_g), order, 8)
    # print(bob)
    #
    # alice_public_key = alice.send_public_key()
    # bob_public_key = bob.send_public_key()
    #
    # alice.get_shared_secret(bob_public_key)
    # bob.get_shared_secret(alice_public_key)
    # print()
    # print(alice.shared_secret)
    # print(bob.shared_secret)
    # print(alice.shared_secret == bob.shared_secret)
