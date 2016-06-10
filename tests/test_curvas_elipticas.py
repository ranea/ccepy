import sys
sys.path.append('../ccepy')
import unittest
from collections import namedtuple
# import doctest

from hypothesis import given, assume  # , settings
from hypothesis.strategies import integers, sampled_from

# from ccepy import cuerpos_finitos
from ccepy.curvas_elipticas import curva_eliptica_sobre_Fq
from ccepy.cuerpos_finitos import Fq

ParametrosCurvaEliptica = namedtuple('CurvaElíptica', ['nombre', 'p', 'a', 'b', 'x1', 'y1'])

# TODO: ref safecurves
curvas_eliptipcas_famosas = [
    ParametrosCurvaEliptica(
        "Anomalous",
        17676318486848893030961583018778670610489016512983351739677143,
        15347898055371580590890576721314318823207531963035637503096292,
        7444386449934505970367865204569124728350661870959593404279615,
        1619092589586542907492569170434842128165755668543894279235270,
        3436949547626524920645513316569700140535482973634182925459687
    ),
    ParametrosCurvaEliptica(
        "NIST P-224",
        26959946667150639794667015087019630673557916260026308143510066298881,
        -3,
        18958286285566608000408668544493926415504680968679321075787234672564,
        19277929113566293071110308034699488026831934219452440156649784352033,
        19926808758034470970197974370888749184205991990603949537637343198772
    ),
    ParametrosCurvaEliptica(
        "BN(2,254)",
        16798108731015832284940804142231733909889187121439069848933715426072753864723,
        0,
        2,
        -1,
        1
    ),
    ParametrosCurvaEliptica(
        "brainpoolP256t1",
        76884956397045344220809746629001649093037950200943055203735601445031516197751,
        -3,
        46214326585032579593829631435610129746736367449296220983687490401182983727876,
        74138526386500101787937404544159543470173440588427591213843535686338908194292,
        20625154686056605250529482107801269759951443923312408063441227608803066104254
    ),
    ParametrosCurvaEliptica(
        "ANSSI FRP256v1",
        109454571331697278617670725030735128145969349647868738157201323556196022393859,
        -3,
        107744541122042688792155207242782455150382764043089114141096634497567301547839,
        82638672503301278923015998535776227331280144783487139112686874194432446389503,
        43992510890276411535679659957604584722077886330284298232193264058442323471611
    ),
    ParametrosCurvaEliptica(
        "NIST P-256",
        115792089210356248762697446949407573530086143415290314195533631308867097853951,
        -3,
        41058363725152142129326129780047268409114441015993725554835256314039467401291,
        48439561293906451759052585252797914202762949526041747995844080717082404635286,
        36134250956749795798585127919587881956611106672985015071877198253568414405109
    ),
    ParametrosCurvaEliptica(
        "secp256k1",
        115792089237316195423570985008687907853269984665640564039457584007908834671663,
        0,
        7,
        55066263022277343669578718895168534326250603453777594175500187360389116729240,
        32670510020758816978083085130507043184471273380659243275938904335757337482424
    ),
    ParametrosCurvaEliptica(
        "brainpoolP384t1",
        21659270770119316173069236842332604979796116387017648600081618503821089934025961822236561982844534088440708417973331,
        -3,
        19596161053329239268181228455226581162286252326261019516900162717091837027531392576647644262320816848087868142547438,
        3827769047710394604076870463731979903132904572714069494181204655675960538951736634566672590576020545838501853661388,
        5797643717699939326787282953388004860198302425468870641753455602553471777319089854136002629714659021021358409132328
    ),
    ParametrosCurvaEliptica(
        "NIST P-384",
        39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319,
        -3,
        27580193559959705877849011840389048093056905856361568521428707301988689241309860865136260764883745107765439761230575,
        26247035095799689268623156744566981891852923491109213387815615900925518854738050089022388053975719786650872476732087,
        8325710961489029985546751289520108179287853048861315594709205902480503199884419224438643760392947333078086511627871
    ),
]


# TODO: rellenar los docstring?
class TestCurvaElipticaFq(unittest.TestCase):
    """Conjuto de test para PuntosFqRacionales"""
    @given(sampled_from(curvas_eliptipcas_famosas), integers(), integers(), integers())
    def test_propiedades_grupo(self, ce, k1, k2, k3):
        E = curva_eliptica_sobre_Fq(ce.a, ce.b, ce.p)
        generador = E(ce.x1, ce.y1)
        print(ce.nombre)

        assume(k1 >= 0)
        assume(k2 >= 0)
        assume(k3 >= 0)

        P = generador * k1
        Q = generador * k2
        R = generador * k3
        O = E.elemento_neutro()

        print("\tP: {0}\n\tQ: {1}\n\tR: {2}".format(P, Q, R))
        assert P + O == P == O + P
        assert P + (-P) == O
        assert P + Q == Q + P
        assert P + (Q + R) == (P + Q) + R

        # TODO: borrar debug
        # x1, y1 = Fpn(x1), Fpn(y1)
        # print("(x1, y1): ({0},{1})".format(x1, y1))
        # assume(E.contiene(x1, y1))
        # print("contiene?:", E.contiene(x1, y1))
        # P = E(x1, y1)
        # print("P:", P)

if __name__ == '__main__':
    # TODO: añadirlo al doctest (añadir tb una curva famosa)
    # E = curva_eliptica_sobre_Q(0, 4)
    # print(E.coeficientes)
    # print("Discriminante:", E.discriminante)
    # P = E(0, -2)
    # print("P =", P)
    # print("-P =", -P)
    # for i in range(3):
    #     print(i, "P =", i * P)
    #
    # E = curva_eliptica_sobre_Q(-18, -72)
    # print(E.coeficientes)
    # print("Discriminante:", E.discriminante)
    # P = E(6, 6)
    # print("P =", P)
    # print("-P =", -P)
    # for i in range(5):
    #     print(i, "P =", i * P)

    # E = curva_eliptica_sobre_Fq(2, 3, 97, 1)
    # print(E.coeficientes)
    # print("Discriminante:", E.discriminante)
    # F97 = Fq(97, 1)
    # P = E(0, 10)
    # print("P =", P)
    # print("-P =", -P)
    # for i in range(5):
    #     print(i, "P =", i * P)
    # Q = E(3, 6)
    # print("Q =", Q)
    # print("P + Q =", P + Q)

    # E = curva_eliptica_sobre_Fq(1, 1, 5, 2)
    # print(E.coeficientes)
    # print("Discriminante:", E.discriminante)
    # F25 = Fq(5, 2)
    # P = E(F25.cero(), F25.uno())
    # print("P =", P)
    # print("-P =", -P)
    # for i in range(5):
    #     print(i, "P =", i * P)
    unittest.main()
