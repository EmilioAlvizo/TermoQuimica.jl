# TermoQuimica.jl
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

Ecuaciones de estado
: Van Der Waals
: Redlich Kwong
: Soave Redlich Kwong
: Peng Robinson

## Ejemplo de uso
Evolucion diferencial pose una funcion llamada 
````
ed(fnc,d,L,h,np,gen,n)
````
Donde: 
* **fnc** es la funcion a optimizar
* **d** es la dimencion de la funcion
* **L** es un vector que contiene los valores inferiores que puede tener las variables de la funcion
* **h** es un vector que contiene los valores superiores que puede tener las variables de la funcion
* **np** es el numero de poblacion
* **gen** es el numero de generaciones que habra
* **n** es el numero de veces que corre el programa

## Autor ✒️

* **Emilio Alvizo Velázquez** - (emilio_alvizo@yahoo.com.mx)

## Expresiones de Gratitud 🎁

* Comenta a otros sobre este proyecto 📢
* Invita una cerveza 🍺 a alguien del equipo. 
* Da las gracias públicamente 🤓.
* etc.