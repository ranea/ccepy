import sys
sys.path.append('../ccepy')
import unittest
import random

from hypothesis import given, assume
from hypothesis.strategies import integers, sampled_from

# from ccepy import esquemas_criptograficos  # para cargar los docstring
from ccepy.esquemas_criptograficos import ECDH
from ccepy.curvas_elipticas import curva_eliptica_sobre_Fq
from ccepy.listado_curvas_elipticas import curvas_eliptipcas_sobre_Fq_famosas
from ccepy.listado_curvas_elipticas import procesar_parametros_curva_eliptica


class TestECDH(unittest.TestCase):
    """Conjuto de test para ECDH"""
    @classmethod
    def setUpClass(cls):
        random.seed(5040)

    @given(sampled_from(curvas_eliptipcas_sobre_Fq_famosas))
    def test(self, ce):
        E, generador, orden = procesar_parametros_curva_eliptica(ce)
        print(ce.nombre)

        alicia = ECDH(E, generador, orden)
        bob = ECDH(E, generador, orden)
        # print("Llaves de Alicia:\n\tPública: {0}\n\tPrivada: {1}".format(alicia.llave_publica,
        #                                                             alicia.llave_privada))
        # print("Llaves de Bob:\n\tPública: {0}\n\tPrivada: {1}".format(bob.llave_publica,
        #                                                             bob.llave_privada))
        alicia.calcula_secreto_compartido(bob.llave_publica)
        bob.calcula_secreto_compartido(alicia.llave_publica)
        # print("Secreto compartido calculado por:\n\tAlicia: {0}\n\tBob: {1}".format(alicia.secreto_compartido,
        #                                                                             bob.secreto_compartido))
        assert alicia.secreto_compartido == bob.secreto_compartido



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
    unittest.main()
