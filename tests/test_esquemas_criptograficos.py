import sys
sys.path.append('../ccepy')
import unittest
import random

from hypothesis import given, assume
from hypothesis.strategies import sampled_from, text

from ccepy.esquemas_criptograficos import ECDH, ECDSA
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

        alicia = ECDH(E, generador, orden)
        bob = ECDH(E, generador, orden)
        secreto_alicia = alicia.calcula_secreto_compartido(bob.llave_publica)
        secreto_bob = bob.calcula_secreto_compartido(alicia.llave_publica)
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

        alicia = ECDSA(E, generador, orden)
        bob = ECDSA(E, generador, orden)
        r, s = alicia.firma(mensaje)
        assert bob.verifica(mensaje, r, s, alicia.llave_publica)

        assume(mensaje != otro_mensaje)
        assert not bob.verifica(otro_mensaje, r, s, alicia.llave_publica)

        eva = ECDSA(E, generador, orden)
        assume(eva.llave_privada != alicia.llave_privada)
        assert not bob.verifica(mensaje, r, s, eva.llave_publica)
        rr, ss = eva.firma(mensaje)
        assert not bob.verifica(mensaje, rr, ss, alicia.llave_publica)


if __name__ == '__main__':
    unittest.main()
