"""Listado de curvas elípticas.

Este módulo ofrece una serie de curvas elípticas listas para utilizar.
Estas curvas se han sacado de varios estándares.

Para utilizar las curvas aquí listadas, debe importar la función
:func:`parametros_dominio` ::

    from ccepy.listado_curvas_elipticas import parametros_dominio

Para trabajar con una curva, basta llamar a esta función pasándole
como parámetro el nombre de la curva a utilizar:

    >>> E, P, n = parametros_dominio("secp256k1")
    >>> E  # la curva elíptica
    ccepy.curvas_elipticas.curva_eliptica_sobre_Fq.<locals>.PuntoFqRacional
    >>> E.Fq.p  # el numero de elementos del cuerpo finito
    115792089237316195423570985008687907853269984665640564039457584007908834671663
    >>> E.coeficientes
    Coeficientes(a=0, b=7)
    >>> P  # un punto de E
    (55066263022277343669578718895168534326250603453777594175500187360389116729240,32670510020758816978083085130507043184471273380659243275938904335757337482424)
    >>> n  # el orden del punto P
    115792089237316195423570985008687907852837564279074904382605163141518161494337

Los nombres de las curvas disponibles son:

- Anomalous
- NIST P-224
- BN(2,254)
- brainpoolP256t1
- ANSSI FRP256v1
- NIST P-256
- secp256k1
- brainpoolP384t1
- NIST P-384
"""
from collections import namedtuple


from ccepy.curvas_elipticas import curva_eliptica_sobre_Fq


ParametrosCurvaElipticaSobreFq = namedtuple('CurvaElípticaSobreFq', ['nombre', 'p', 'a', 'b', 'x1', 'y1', 'orden'])


# Listado de curvas
curvas_eliptipcas_sobre_Fq_famosas = [
    ParametrosCurvaElipticaSobreFq(
        "Anomalous",  # nombre
        17676318486848893030961583018778670610489016512983351739677143,  # p
        15347898055371580590890576721314318823207531963035637503096292,  # a
        7444386449934505970367865204569124728350661870959593404279615,  # b
        1619092589586542907492569170434842128165755668543894279235270,  # x1
        3436949547626524920645513316569700140535482973634182925459687,  # y1
        17676318486848893030961583018778670610489016512983351739677143,  # orden
    ),
    ParametrosCurvaElipticaSobreFq(
        "NIST P-224",
        26959946667150639794667015087019630673557916260026308143510066298881,
        -3,
        18958286285566608000408668544493926415504680968679321075787234672564,
        19277929113566293071110308034699488026831934219452440156649784352033,
        19926808758034470970197974370888749184205991990603949537637343198772,
        26959946667150639794667015087019625940457807714424391721682722368061
    ),
    ParametrosCurvaElipticaSobreFq(
        "BN(2,254)",
        16798108731015832284940804142231733909889187121439069848933715426072753864723,
        0,
        2,
        -1,
        1,
        16798108731015832284940804142231733909759579603404752749028378864165570215949,
    ),
    ParametrosCurvaElipticaSobreFq(
        "brainpoolP256t1",
        76884956397045344220809746629001649093037950200943055203735601445031516197751,
        -3,
        46214326585032579593829631435610129746736367449296220983687490401182983727876,
        74138526386500101787937404544159543470173440588427591213843535686338908194292,
        20625154686056605250529482107801269759951443923312408063441227608803066104254,
        76884956397045344220809746629001649092737531784414529538755519063063536359079,
    ),
    ParametrosCurvaElipticaSobreFq(
        "ANSSI FRP256v1",
        109454571331697278617670725030735128145969349647868738157201323556196022393859,
        -3,
        107744541122042688792155207242782455150382764043089114141096634497567301547839,
        82638672503301278923015998535776227331280144783487139112686874194432446389503,
        43992510890276411535679659957604584722077886330284298232193264058442323471611,
        109454571331697278617670725030735128146004546811402412653072203207726079563233,
    ),
    ParametrosCurvaElipticaSobreFq(
        "NIST P-256",
        115792089210356248762697446949407573530086143415290314195533631308867097853951,
        -3,
        41058363725152142129326129780047268409114441015993725554835256314039467401291,
        48439561293906451759052585252797914202762949526041747995844080717082404635286,
        36134250956749795798585127919587881956611106672985015071877198253568414405109,
        115792089210356248762697446949407573529996955224135760342422259061068512044369,
    ),
    ParametrosCurvaElipticaSobreFq(
        "secp256k1",
        115792089237316195423570985008687907853269984665640564039457584007908834671663,
        0,
        7,
        55066263022277343669578718895168534326250603453777594175500187360389116729240,
        32670510020758816978083085130507043184471273380659243275938904335757337482424,
        115792089237316195423570985008687907852837564279074904382605163141518161494337,
    ),
    ParametrosCurvaElipticaSobreFq(
        "brainpoolP384t1",
        21659270770119316173069236842332604979796116387017648600081618503821089934025961822236561982844534088440708417973331,
        -3,
        19596161053329239268181228455226581162286252326261019516900162717091837027531392576647644262320816848087868142547438,
        3827769047710394604076870463731979903132904572714069494181204655675960538951736634566672590576020545838501853661388,
        5797643717699939326787282953388004860198302425468870641753455602553471777319089854136002629714659021021358409132328,
        21659270770119316173069236842332604979796116387017648600075645274821611501358515537962695117368903252229601718723941
    ),
    ParametrosCurvaElipticaSobreFq(
        "NIST P-384",
        39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319,
        -3,
        27580193559959705877849011840389048093056905856361568521428707301988689241309860865136260764883745107765439761230575,
        26247035095799689268623156744566981891852923491109213387815615900925518854738050089022388053975719786650872476732087,
        8325710961489029985546751289520108179287853048861315594709205902480503199884419224438643760392947333078086511627871,
        39402006196394479212279040100143613805079739270465446667946905279627659399113263569398956308152294913554433653942643,
    ),
]


def procesar_parametros_curva_eliptica(parametros_curva_eliptica):
    ce = parametros_curva_eliptica
    E = curva_eliptica_sobre_Fq(ce.a, ce.b, ce.p)
    generador = E(ce.x1, ce.y1)
    orden = ce.orden
    return E, generador, orden


def parametros_dominio(nombre):
    """Devuelve los parámetros de dominio de una curva en el listado.:

        >>> E, P, n = parametros_dominio("secp256k1")
        >>> E.Fq.p
        115792089237316195423570985008687907853269984665640564039457584007908834671663
        >>> E.coeficientes
        Coeficientes(a=0, b=7)
        >>> P
        (55066263022277343669578718895168534326250603453777594175500187360389116729240,32670510020758816978083085130507043184471273380659243275938904335757337482424)
        >>> n
        115792089237316195423570985008687907852837564279074904382605163141518161494337

    Los nombres de las curvas en el listado son:

    - Anomalous
    - NIST P-224
    - BN(2,254)
    - brainpoolP256t1
    - ANSSI FRP256v1
    - NIST P-256
    - secp256k1
    - brainpoolP384t1
    - NIST P-384

    Args:
        nombre (str): uno de los nombres de la lista de arriba.

    Returns:
        Tuple: una tupla (E, P, n) donde E es una curva elítipca, P un punto
        y n el orden de P. Si no existe ninguna curva con dicho nombre,
        se devuelve ``None``.
    """
    for param_dominio in curvas_eliptipcas_sobre_Fq_famosas:
        if param_dominio.nombre == nombre:
            return procesar_parametros_curva_eliptica(param_dominio)
    else:
        return None
