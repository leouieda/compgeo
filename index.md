*^a*Professores:** [Leonardo Uieda](http://www.leouieda.com/) e Eder Molina

## Conteúdo da parte em Python

Nesta parte da matéria aprenderemos técnicas de engenharia de software e programação em Python em um nível mais avançado. O problema geofísico que utilizaremos de motivação é o cálculo do campo geomagnético segundo o modelo [IGRF](https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field).

O IGRF é um modelo de harmônicos esféricos. Com esses modules, podemos calcular as 3 componentes do campo geomagnético em qualquer local da Terra através de uma tabela de números chamados **coeficientes de Gauss**. Para um determinado conjunto de coeficientes de Gauss $g_n^m$ e $h_n^m$, as 3 componentes do campo geomagnético são:

$$
B_n(r, \theta, \lambda) = -\dfrac{R\sin\theta}{r} \sum\limits_{n=1}^{N}\sum\limits_{m=0}^{n} \left(\dfrac{R}{r}\right)^{n+1} [ g_n^m \cos m\lambda + h_n^m \sin m\lambda ] \dfrac{\partial P_n^m(\cos\theta)}{\partial \cos\theta}
$$

$$
B_e(r, \theta, \lambda) = -\dfrac{R}{r\sin\theta} \sum\limits_{n=1}^{N}\sum\limits_{m=0}^{n} \left(\dfrac{R}{r}\right)^{n+1} [ -m g_n^m \sin m\lambda + m h_n^m \cos m\lambda ] P_n^m(\cos\theta)
$$

$$
B_r(r, \theta, \lambda) = \sum\limits_{n=1}^{N}\sum\limits_{m=0}^{n} (n + 1)\left(\dfrac{R}{r}\right)^{n+2} [ g_n^m \cos m\lambda + h_n^m \sin m\lambda ] P_n^m(\cos\theta)
$$

em que $r$ é a direção radial, $\theta$ é a colatitude geocêntrica, $\lambda$ é a longitude, $n$ é o grau, $m$ é a ordem, $R$ é o raio médio da Terra e $P_n^m(x)$ são [funções associadas de Legendre](https://en.wikipedia.org/wiki/Associated_Legendre_polynomials).

A tabela com os coeficientes de Gauss para cada 5 anos pode ser baixada do [site do NOAA](https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf13coeffs.txt). Com os coeficientes, podemos utilizar as fórmulas acima para calcular o campo geomagnético para qualquer data. Com a tabela de coeficientes também podemos calcular o momento magnético de um dipolo geocêntrico:

$$
m_X = \dfrac{4\pi}{\mu_0} R^3 g_1^1
$$

$$
m_Y = \dfrac{4\pi}{\mu_0} R^3 h_1^1
$$

$$
m_Z = \dfrac{4\pi}{\mu_0} R^3 g_1^0
$$

O momento pode ser convertido para coordenadas esféricas e providenciar a longitude e latitude do polo geomagnético.

As [funções associadas de Legendre](https://en.wikipedia.org/wiki/Associated_Legendre_polynomials) podem ser calculadas com as fórmulas recursivas:

$$
P_n^m (x) = \dfrac{(2n - 1) x P_{n-1}^m - (n + m - 1)P_{n-2}^m}{n - m}\quad \text{para}\ m < n - 1
$$

$$
P_n^m (x) = \dfrac{-\sqrt{1 - x^2}P_n^{m+1} - \sqrt{1-x^2}(n + m)(n - m + 1)P_n^{m-1}}{2mx}\quad \text{para}\ m = n - 1
$$

$$
P_n^m (x) = -(2n -1)\sqrt{1-x^2}P_{n-1}^{m-1} \quad \text{para}\ m = n
$$

Sabendo os valores iniciais: $P_0^0 = 1$, $P_1^0 = x$, $P_1^1=-\sqrt{1 - x^2}$

As derivadas $\partial P_n^m / \partial x$ também são calculadas assim:

$$
\dfrac{\partial P_n^m}{\partial x} (x) = \dfrac{n x P_n^m - (n + m) P_{n-1}^m}{x^2 - 1} \quad \text{para}\ m \neq n
$$

$$
\dfrac{\partial P_n^m}{\partial x} (x) = \dfrac{-(n + m)(n - m + 1)\sqrt{1 - x^2} P_n^{m-1} - m x P_n^m}{x^2 - 1} \quad \text{para}\ m = n
$$

com os valores iniciais: $\partial P_0^0 / \partial x=0$, $\partial P_1^0 / \partial x=1$, $\partial P_1^1 / \partial x=\frac{x}{\sqrt{1 - x^2}}$.

Nas aplicações em geomagnetismo, as funções de Legendre e suas derivadas são normalizadas com a normalização de Schmidt multiplicando-os pelo fator:

$$
S_n^m = \sqrt{(2 - \delta_{m,0)\dfrac{(n - m)!}{(n + m)!}}
$$

na qual $\delta_(m,0)$ é o [delta de Kronecker](https://en.wikipedia.org/wiki/Kronecker_delta).

Esse problema pode ser decomposto em diversas etapas computacionais. Melhor ainda, cada etapa pode ser encapsulada em 1 ou mais funções que podem ser testadas independentemente. **Ao final dessa disciplina, teremos funções e programas de linha de comando que serão capazes de calcular o campo geomagnético em qualquer lugar da Terra e em qualquer data entre 1900 e o presente.**


## Cronograma

> O cronograma provavelmente sofrerá alterações ao longo do semestre.

| Semana | Tema                                 | Produto |
|:------:|:-------------------------------------|:--------|
| 6      | Apresentação do problema / Revisão de Python / Lendo dados de arquivos | Código que lê o arquivo de coeficientes de Gauss |
| 7      | Programação defensiva / Refatoração de código em funções | Testes para parâmetros de entrada e saída + função com o código da aula anterior |
| 8      | Cálculo do polo magnético / Mapas com [PyGMT](https://www.pygmt.org) | Função que calcula o momento de dipolo + latitude e longitude do polo magnético dadas coeficientes de Gauss |
| 9      | Lindando com datas e interpolação linear | Função que produz coeficientes de Gauss para qualquer data |
| 10     | As tais das funções de Legendre | Função que retorna as funções associadas de Legendre |
| 11     | Calculando o campo magnético | Função que calculam Be, Bn, Bu para uma data |
| 12     | Malhas regulares em Python / xarray + netCDF | Função que produz malhas regulares do campo magnético |
| 13     | Optimização e perfilagem de código / acelerando as funções de Legendre com numba | Função que calcula as funções de Legendre de forma mais rápida |
| 14     | Programas de linha de comando em Python | Transformação das funções em um programa de linha de comando |
| 15     | Livre |
