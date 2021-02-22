# TermoQuimica.jl 🧪

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
|[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://EmilioAlvizo.github.io/TermoQuimica.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://EmilioAlvizo.github.io/TermoQuimica.jl/dev)| [![Build Status](https://github.com/EmilioAlvizo/TermoQuimica.jl/workflows/CI/badge.svg)](https://github.com/EmilioAlvizo/TermoQuimica.jl/actions) [![Coverage](https://codecov.io/gh/EmilioAlvizo/TermoQuimica.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/EmilioAlvizo/TermoQuimica.jl)|


Programa para el equilibrio liquido/vapor y calculo de propiedades usando ecuaciones cubicas de estado.


## Caracteristicas
TermoQuimica puede calcular el equilibrio liquido y vapor de n componentes, haciendo uso de distintos métodos y sin hacer molestas conversiones de unidades.

Métodos utilizados:
* Ideal
* Van Laar
* Margules
* Wilson
* Wilson mod.
* NRTL
* UNIQUAC
* UNIFAC

Ecuaciones de estado:
* Van Der Waals
* Redlich Kwong
* Soave Redlich Kwong
* Peng Robinson

## Ejemplo de uso
Para el equilibrio **P-X-Y** la mezcla **etanol(1)-agua(2)** a 273.15 K, usando el modelo de actividad de **Van Laar**, las coeficientes de Van Laar son: `A12 = 1.6798` y `A21 = 0.9227`, las constantes de Antoine para el etanol son: `A = 8.12875`, `B = 1660.8713`, `C = 238.131` y para el agua `A = 8.05573`, `B = 1723.6425`, `C = 233.08`

````
vanlaar.pxy(n,T,CA,Λ,xx;uni=u"Torr")
````
Donde: 
* **n** Es el numero de componentes
* **T** Es la temperatura del sistema
* **CA** Constantes de Antoine para los `n` componentes
* **Λ** Parámetros de interacción entre los componentes en grados kelvin
* **xx** Son los puntos en el liquido donde se buscara el equilibrio con el vapor
* **uni** Son las unidades en las que se desea el resultado, por defecto esta en Torr

Entonces se deberá hacer lo siguiente
````
n = 2
T = 273.15u"K"
CA = [8.12875, 1660.8713, 238.131; 8.05573, 1723.6425, 233.08]
Λ = [0 , 1.6798; 0.9227, 0]
xx = 0:0.1:1
unidad = "Pa"

vanlaar.pxy(n,T,CA,Λ,xx;uni=u"Torr")
````

## Autor ✒️

* **Emilio Alvizo Velázquez** - (emilio_alvizo@yahoo.com.mx)

## Expresiones de Gratitud 🎁

* Comenta a otros sobre este proyecto 📢
* Invita una cerveza 🍺 a alguien del equipo. 
* Da las gracias públicamente 🤓.
* etc.