import sys
sys.path.append('../ccepy')
import unittest
import random

from hypothesis import given, assume
from hypothesis.strategies import integers, sampled_from, text

# from ccepy import esquemas_criptograficos  # para cargar los docstring
from ccepy.esquemas_criptograficos import ECDH, ECDSA
from ccepy.curvas_elipticas import curva_eliptica_sobre_Fq
from ccepy.listado_curvas_elipticas import curvas_eliptipcas_sobre_Fq_famosas
from ccepy.listado_curvas_elipticas import procesar_parametros_curva_eliptica


class TestECDH(unittest.TestCase):
    """Conjuto de test para ECDH"""
    @classmethod
    def setUpClass(cls):
        random.seed(5040)

    @given(sampled_from(curvas_eliptipcas_sobre_Fq_famosas))
    def test_ECDH(self, ce):
        E, generador, orden = procesar_parametros_curva_eliptica(ce)
        # print(ce.nombre)

        alicia = ECDH(E, generador, orden)
        bob = ECDH(E, generador, orden)
        # print("Llaves de Alicia:\n\tPública: {0}\n\tPrivada: {1}".format(alicia.llave_publica,
        #                                                             alicia.llave_privada))
        # print("Llaves de Bob:\n\tPública: {0}\n\tPrivada: {1}".format(bob.llave_publica,
        #                                                             bob.llave_privada))
        secreto_alicia = alicia.calcula_secreto_compartido(bob.llave_publica)
        secreto_bob = bob.calcula_secreto_compartido(alicia.llave_publica)
        # print("Secreto compartido calculado por:\n\tAlicia: {0}\n\tBob: {1}".format(alicia.secreto_compartido,
        #                                                                             bob.secreto_compartido))
        assert secreto_alicia == secreto_bob

        eva = ECDH(E, generador, orden)  # eva impersona a alicia
        assume(eva.llave_privada != alicia.llave_privada)
        assert secreto_alicia != eva.calcula_secreto_compartido(bob.llave_publica)


class TestECDSA(unittest.TestCase):
    """Conjuto de test para ECDSA"""
    @classmethod
    def setUpClass(cls):
        random.seed(5040)

    @given(sampled_from(curvas_eliptipcas_sobre_Fq_famosas), text(), text())
    def test_ECDSA(self, ce, mensaje, otro_mensaje):
        E, generador, orden = procesar_parametros_curva_eliptica(ce)
        # print(ce.nombre, "|", mensaje, "|")

        alicia = ECDSA(E, generador, orden)
        bob = ECDSA(E, generador, orden)
        # print("\n\nMensaje:", mensaje)
        # print("Llaves de Alicia:\n\tPública: {0}\n\tPrivada: {1}".format(alicia.llave_publica,
        #                                                             alicia.llave_privada))
        r, s = alicia.firma(mensaje)
        # print("r:", r, "\ns:", s)
        assert bob.verifica(mensaje, r, s, alicia.llave_publica)

        assume(mensaje != otro_mensaje)
        # print("Otro mensaje:", otro_mensaje)
        assert not bob.verifica(otro_mensaje, r, s, alicia.llave_publica)

        eva = ECDSA(E, generador, orden)
        assume(eva.llave_privada != alicia.llave_privada)
        assert not bob.verifica(mensaje, r, s, eva.llave_publica)
        rr, ss = eva.firma(mensaje)
        assert not bob.verifica(mensaje, rr, ss, alicia.llave_publica)


if __name__ == '__main__':
    # random.seed(100)
    # E = curva_eliptica_sobre_Fq(a=324, b=1287, p=3851)
    # generador = E(920, 303)
    # # TODO: añadir metodo orden?
    # orden = 8
    #
    # alicia = ECDH(E, generador, orden)
    # bob = ECDH(E, generador, orden)
    # print("Llaves de Alicia:\n\tPública: {0}\n\tPrivada: {1}".format(alicia.llave_publica,
    #                                                             alicia.llave_privada))
    # print("Llaves de Bob:\n\tPública: {0}\n\tPrivada: {1}".format(bob.llave_publica,
    #                                                             bob.llave_privada))
    #
    # alicia.calcula_secreto_compartido(bob.llave_publica)
    # bob.calcula_secreto_compartido(alicia.llave_publica)
    # print("Secreto compartido calculado por:\n\tAlicia: {0}\n\tBob: {1}".format(alicia.secreto_compartido,
    #                                                                             bob.secreto_compartido))
    # print("¿Es el mismo secreto? {0}".format(alicia.secreto_compartido == bob.secreto_compartido))

    E = curva_eliptica_sobre_Fq(a=5, b=7, p=113)
    generador = E(16, 51)

    def orden_punto_curva_eliptica(punto):
        i = 1
        P = 1 * punto  # una copia
        for j in range(1, 128):
            if P.es_elemento_neutro():
                return i
            P += punto
            i += 1

    print(orden_punto_curva_eliptica(generador))
    orden = orden_punto_curva_eliptica(generador)
    mensaje = "It is possible to write endlessly on elliptic curves."
    alicia = ECDSA(E, generador, orden)
    print(alicia.llave_publica, alicia.llave_privada)
    bob = ECDSA(E, generador, orden)

    r, s = alicia.firma(mensaje)
    print(r, s)
    print(bob.verifica(mensaje, r, s, alicia.llave_publica))

    # unittest.main()
