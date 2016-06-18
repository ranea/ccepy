import sys
sys.path.append('../ccepy')
import unittest
from collections import namedtuple
from fractions import Fraction
import doctest

from hypothesis import given, assume  # , settings
from hypothesis.strategies import integers, sampled_from

from ccepy import curvas_elipticas  # para cargar los docstring
from ccepy.curvas_elipticas import curva_eliptica_sobre_Fq, curva_eliptica_sobre_Q, curva_eliptica_sobre_F2m
from ccepy.cuerpos_finitos import Fq, PolinomioZp

ParametrosCurvaElipticaSobreFq = namedtuple('CurvaElípticaFq', ['nombre', 'p', 'a', 'b', 'x1', 'y1'])

# TODO: ref safecurves
curvas_eliptipcas_sobre_Fq_famosas = [
    ParametrosCurvaElipticaSobreFq(
        "Anomalous",
        17676318486848893030961583018778670610489016512983351739677143,
        15347898055371580590890576721314318823207531963035637503096292,
        7444386449934505970367865204569124728350661870959593404279615,
        1619092589586542907492569170434842128165755668543894279235270,
        3436949547626524920645513316569700140535482973634182925459687
    ),
    ParametrosCurvaElipticaSobreFq(
        "NIST P-224",
        26959946667150639794667015087019630673557916260026308143510066298881,
        -3,
        18958286285566608000408668544493926415504680968679321075787234672564,
        19277929113566293071110308034699488026831934219452440156649784352033,
        19926808758034470970197974370888749184205991990603949537637343198772
    ),
    ParametrosCurvaElipticaSobreFq(
        "BN(2,254)",
        16798108731015832284940804142231733909889187121439069848933715426072753864723,
        0,
        2,
        -1,
        1
    ),
    ParametrosCurvaElipticaSobreFq(
        "brainpoolP256t1",
        76884956397045344220809746629001649093037950200943055203735601445031516197751,
        -3,
        46214326585032579593829631435610129746736367449296220983687490401182983727876,
        74138526386500101787937404544159543470173440588427591213843535686338908194292,
        20625154686056605250529482107801269759951443923312408063441227608803066104254
    ),
    ParametrosCurvaElipticaSobreFq(
        "ANSSI FRP256v1",
        109454571331697278617670725030735128145969349647868738157201323556196022393859,
        -3,
        107744541122042688792155207242782455150382764043089114141096634497567301547839,
        82638672503301278923015998535776227331280144783487139112686874194432446389503,
        43992510890276411535679659957604584722077886330284298232193264058442323471611
    ),
    ParametrosCurvaElipticaSobreFq(
        "NIST P-256",
        115792089210356248762697446949407573530086143415290314195533631308867097853951,
        -3,
        41058363725152142129326129780047268409114441015993725554835256314039467401291,
        48439561293906451759052585252797914202762949526041747995844080717082404635286,
        36134250956749795798585127919587881956611106672985015071877198253568414405109
    ),
    ParametrosCurvaElipticaSobreFq(
        "secp256k1",
        115792089237316195423570985008687907853269984665640564039457584007908834671663,
        0,
        7,
        55066263022277343669578718895168534326250603453777594175500187360389116729240,
        32670510020758816978083085130507043184471273380659243275938904335757337482424
    ),
    ParametrosCurvaElipticaSobreFq(
        "brainpoolP384t1",
        21659270770119316173069236842332604979796116387017648600081618503821089934025961822236561982844534088440708417973331,
        -3,
        19596161053329239268181228455226581162286252326261019516900162717091837027531392576647644262320816848087868142547438,
        3827769047710394604076870463731979903132904572714069494181204655675960538951736634566672590576020545838501853661388,
        5797643717699939326787282953388004860198302425468870641753455602553471777319089854136002629714659021021358409132328
    ),
    ParametrosCurvaElipticaSobreFq(
        "NIST P-384",
        39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319,
        -3,
        27580193559959705877849011840389048093056905856361568521428707301988689241309860865136260764883745107765439761230575,
        26247035095799689268623156744566981891852923491109213387815615900925518854738050089022388053975719786650872476732087,
        8325710961489029985546751289520108179287853048861315594709205902480503199884419224438643760392947333078086511627871
    ),
]

# TODO: ver que hacer con esto
# ParametrosCurvaElipticaSobreF2m = namedtuple('CurvaElípticaF2m', ['nombre', 'm', 'pol_irreducible', 'a', 'b', 'x1', 'y1'])
#
# # TODO: ref guide to elliptic curve
# curvas_eliptipcas_sobre_F2m_famosas = [
#     ParametrosCurvaElipticaSobreF2m(
#         "B-163",
#         163,
#         "z^163 + z^7 + z^6 + z^3 + 1",
#         1,
#         0x000000020A601907B8C953CA1481EB10512F78744A3205FD,
#         0x00000003F0EBA16286A2D57EA0991168D4994637E8343E36,
#         0x00000000D51FBC6C71A0094FA2CDD545B11C5C0C797324F1
#     ),
# ]


class TestCurvaElipticaFq(unittest.TestCase):
    """Conjuto de test para PuntosFqRacionales"""
    @given(sampled_from(curvas_eliptipcas_sobre_Fq_famosas), integers(), integers(), integers())
    def test_propiedades_grupo(self, ce, k1, k2, k3):
        assume(k1 >= 0)
        assume(k2 >= 0)
        assume(k3 >= 0)

        E = curva_eliptica_sobre_Fq(ce.a, ce.b, ce.p)
        generador = E(ce.x1, ce.y1)
        # print(ce.nombre)

        P = generador * k1
        Q = generador * k2
        R = generador * k3
        O = E.elemento_neutro()

        # print("\tP: {0}\n\tQ: {1}\n\tR: {2}".format(P, Q, R))
        assert P + O == P == O + P
        assert P + (-P) == O
        assert P + Q == Q + P
        assert P + (Q + R) == (P + Q) + R

    @given(sampled_from(curvas_eliptipcas_sobre_Fq_famosas), integers(), integers(min_value=-100, max_value=100))
    def test_multiplicacion_por_duplicacion(self, ce, k, e):
        assume(k >= 0)
        assume(e)

        E = curva_eliptica_sobre_Fq(ce.a, ce.b, ce.p)
        generador = E(ce.x1, ce.y1)
        # print(ce.nombre)

        P = generador * k

        multiplicacion = E.elemento_neutro()
        if e < 0:
            for i in range(-e):
                multiplicacion += -P
        elif e > 0:
            for i in range(e):
                multiplicacion += P

        # print("\tP: {0}\n\tP': {1}".format(P, multiplicacion))
        assert P * e == multiplicacion


# TODO: decidir si poner o no
# class TestCurvaElipticaF2m(unittest.TestCase):
#     """Conjuto de test para PuntosF2mRacionales"""
#     @given(sampled_from(curvas_eliptipcas_sobre_F2m_famosas))
#     def setUp(self, parametros_curva_eliptica):
#         nombre = parametros_curva_eliptica.nombre
#         m = parametros_curva_eliptica.m
#
#         rep_pol_irreducible = parametros_curva_eliptica.pol_irreducible
#         monomios = rep_pol_irreducible.split(" + ")
#         pol_irreducible = PolinomioZp([0], p=2)
#         for monomio in monomios:
#             if monomio == "1":
#                 pol_irreducible += PolinomioZp([1], p=2)
#                 continue
#
#             grado = int(monomio.split("^")[1])
#             pol_irreducible += PolinomioZp.monomio(1, grado, 2)
#
#         a = parametros_curva_eliptica.a
#
#         F2m = Fq(2, m, pol_irreducible)
#
#         b = parametros_curva_eliptica.b
#         b = F2m([int(i) for i in reversed(bin(b)[2:])])
#         x1 = parametros_curva_eliptica.x1
#         x1 = F2m([int(i) for i in reversed(bin(x1)[2:])])
#         y1 = parametros_curva_eliptica.y1
#         y1 = F2m([int(i) for i in reversed(bin(y1)[2:])])
#
#         self.nombre, self.m, self.pol_irreducible = nombre, m, pol_irreducible
#         self.a, self.b, self.x1, self.y1 = a, b, x1, y1
#
#     @given(integers(min_value=0, max_value=3), integers(min_value=0, max_value=3), integers(min_value=0, max_value=3))
#     def test_propiedades_grupo(self, k1, k2, k3):
#         assume(k1 >= 0)
#         assume(k2 >= 0)
#         assume(k3 >= 0)
#
#         nombre, m, pol_irreducible = self.nombre, self.m, self.pol_irreducible
#         a, b, x1, y1 = self.a, self.b, self.x1, self.y1
#
#         E = curva_eliptica_sobre_F2m(a, b, m, pol_irreducible)
#         generador = E(x1, y1)
#         print(nombre)
#
#         P = generador * k1
#         Q = generador * k2
#         R = generador * k3
#         O = E.elemento_neutro()
#
#         print("\tP: {0}\n\tQ: {1}\n\tR: {2}".format(P, Q, R))
#         assert P + O == P == O + P
#         assert P + (-P) == O
#         assert P + Q == Q + P
#         assert P + (Q + R) == (P + Q) + R
#
#     @given(integers(min_value=-10, max_value=10), integers(min_value=-10, max_value=10))
#     def test_multiplicacion_por_duplicacion(self, k, e):
#         assume(k >= 0)
#         assume(e)
#
#         nombre, m, pol_irreducible = self.nombre, self.m, self.pol_irreducible
#         a, b, x1, y1 = self.a, self.b, self.x1, self.y1
#
#         E = curva_eliptica_sobre_F2m(a, b, m, pol_irreducible)
#         generador = E(x1, y1)
#         print(nombre)
#
#         P = generador * k
#
#         multiplicacion = E.elemento_neutro()
#         if e < 0:
#             for i in range(-e):
#                 multiplicacion += -P
#         elif e > 0:
#             for i in range(e):
#                 multiplicacion += P
#
#         print("\tP: {0}\n\tP': {1}".format(P, multiplicacion))
#         assert P * e == multiplicacion


# Añade los ejemplos de los docstring
def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(curvas_elipticas))
    return tests


if __name__ == '__main__':
    unittest.main()
