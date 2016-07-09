Primeros pasos
**************

===========
Instalación
===========

Para instalar la última versión de *ccepy*, simplemente ejecute: ::

    pip install ccepy

*pip* se encuentra ya instalado por defecto en las últimas distribuciones
de python. Si no es su caso, `aquí
<https://packaging.python.org/installing/>`_
se explica como instalar pip.


===
Uso
===

Una vez instalado *ccepy*, para utilizarlo basta importarlo
y usar las distintas funciones o clases presentes.

Por ejemplo, dada la curva elíptica definida por
la ecuación :math:`y^2 = x^3 + 2 x + 3` sobre el cuerpo
finito de 97 elementos :math:`\mathbb{F}_{97}`, el siguiente trozo
de código calcula la suma :math:`(0, 10) + (3, 6)`::

    >>> from ccepy.curvas_elipticas import curva_eliptica_sobre_Fq
    >>> E = curva_eliptica_sobre_Fq(a=2, b=3, p=97)
    >>> E(0, 10) + E(3, 6)
    (85,71)

O por otro lado, supongamos que quiere simular el protocolo
Diffie-Hellman de intercambio de llaves utilizando como
parámetros de dominio la curva elíptica  :math:`y^2 = x^3 + 324 x + 1287`
definida sobre :math:`\mathbb{F}_{3851}` y el subgrupo generado
por el punto :math:`P = (0, 10)` de orden :math:`8`. ::

    >>> from ccepy.curvas_elipticas import curva_eliptica_sobre_Fq
    >>> from ccepy.esquemas_criptograficos import ECDH
    >>> # definimos los parámetros de dominio
    >>> E = curva_eliptica_sobre_Fq(a=324, b=1287, p=3851)
    >>> P = E(920, 303)
    >>> orden_P = 8
    >>> # definimos los participantes
    >>> alicia = ECDH(E, P, orden_P)
    >>> bob = ECDH(E, P, orden_P)
    >>> # la pareja de llaves se genera automáticamente y aleatoriamente
    >>> alicia.llave_publica
    (2373,2607)
    >>> alicia.llave_privada
    2
    >>> alicia.calcula_secreto_compartido(bob.llave_publica)
    1136
    >>> bob.calcula_secreto_compartido(alicia.llave_publica)
    1136

Puede encontrar la funcionalidad completa de esta programa
estructurada en los siguientes módulos:

.. autosummary::
   aritmetica_elemental
   cuerpos_finitos
   curvas_elipticas
   esquemas_criptograficos
   listado_curvas_elipticas

Le recomendamos que lea la documentación de cada módulo en este
orden para aprender totalmente a usar *ccepy*.
