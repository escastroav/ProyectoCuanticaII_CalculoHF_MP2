# Proyecto Cuantica II: Calculo Hartree-Fock para una molecula diatomica aplicando teoria de perturbacion a segundo orden de Moller-Plesset

## *Resumen*

En este proyecto se ha creado un programa que calcula la energía de correlación electrónica, para un sistema cuántico compuesto por dos núcleos y dos electrones, o bien una molécula diatómica. Se tomaron como ejemplo las moléculas He H+ y H2. Este cálculo es bastante complicado de realizar de forma analítica, ya que tenemos potenciales de Coulomb electrón-electrón y electrón-núcleo los cuales complican la ecuación de Schrödinguer. Para este tipo de sistemas se recurre a un método variacional conocido como método Hartree-Fock, el cual consiste en tomar el potencial de Coulomb electrón-electrón como la interacción entre un sólo electrón y una distribución de carga provocada por los demás electrones en el sistema. Posterior a ello se elige unas funciones de onda iniciales e hidrogenoides, concebidas como determinantes de Slater para mantener sus propiedades antisimétricas, las cuales harán parte del cálculo de campo autoconsistente (SCF). Aunque con este método es posible obtener los orbitales, se considerará únicamente la energía de correlación electrónica para aplicar una perturbación de segundo orden de Moller-Plesset, la cual realizará una corrección de la energía dependiendo del sistema.

## *Introduccion*

La ecuación de Schrödinguer para un sistema diatómico de varios electrones es en general de la forma

(-nabla2_i/2 - nabla2_A/2 + Za*Zb/RAB + Za/rAi + 1/rij)*phi(r) = E_i*phi(r)

Acudiendo a la aproximación Born-Oppenheimer, la cual permite la separabilidad del Hamiltoniano entre coordenadas de los electrones y de los núcleos, asumiendo el hecho de que el núcleo se encuentra aproximadamente quieto, permite despreciar el término cinético de los núcleos y el potencial núcleo-núcleo, obteniendo así 

(-nabla2_i/2 + Za/rAi + 1/rij)*phi(r) = E_i*phi(r)

Dado que los electrones son fermiones que deben satisfacer el principio de exclusión de Pauli, la función de onda del sistema no basta con definirse a partir de una productoria de las funciones de cada electrón, ya que no cumpliría con la antisimetría de intercambio de partículas en un estado. Lo cual requiere hacer uso de la determinante de Slater para definir la función de onda completa

PHI(r1,2,..) = 1/sqrt(N) Det(Slater..) = phi_i(ri)

### Ecuación de Hartree-Fock
Se recurre a la primera aproximación que aborda el método Hartree-Fock, el cual consiste en considerar el término de potencial Coulombiano electrón-electrón como una interacción de un sólo electrón con una distribución de carga promedio. Por cada estado phi_i(r_i) hay dos electrones con diferente spin, éstos interactúan con una distribución de carga realizando el valor promedio del potencial entre un electrón con posición instantánea r_i y una distribución ubicada en r_j la cual ocupa un diferencial de volumen dr y se pesa con la densidad de probabilidad |phi_j(r)|^2. Luego el operador de Coulomb *J* se define como

*J_i* = < phi_j | (r_ij)^-1 | phi_j > = integral...

Ahora bien, debido al intercambio de electrones se requiere considerar un operador de intercambio debido a la antisimetría de los orbitales electrónicos. Este operador tiene la forma

*H_i*phi_i(r_i) = < phi_i | (r_ij)^-1 | phi_j >phi_i(r_i) = integral...

Con esto se escribe el potencial electrón-electrón V_hf, teniendo en cuenta que hay 2 electrones por estado, como

V_hf = SUM_i(2*J_i*-*H_i*)

Por último definiendo el Hamiltoniano base H_core = -nabla2_i/2 + Za/rAi se define el operador de Fock como

F = H_core + V_hf 

El cual se utiliza para escribir la ecuación de Hartree-Fock

F phi_i(r_i) = E_i phi_i(r_i)

Es posible obtener esta ecuación a partir de la extremal de una funcional, sin embargo se expone una manera más sencilla y basada en argumentos físicos para obtenerla[1].

### Escogencia de orbitales

Las funciones de onda a utilizar se pueden escribir como una combinación lineal de funciones base hidrogenoides

psi(r) = SUM_i c_i phi(ri)

En este caso se utilizarían funciones tipo Slater asociadas al orbital 1S, la cual tiene la foma

PHI SLA ...

Sin embargo esta función base requiere un coste computacional alto a la hora de calcular las integrales, por lo que es más conveniente escribir la función base como una combinación lineal de Gaussianas que dependen de un parámetro genérico alpha y unos coeficientes d_n, obtenidos realizando un método variacional que converga al orbital hidrogenide deseado, denominada por ejemplo STO-G3 para la superposición de 3 gaussianas. Esta función suele expresarse como

PHI STO-3G

y serán las utilizadas como nuevas bases para realizar el cálculo de las integrales, ahorrando tiempo de cómputo. 

### Cálculo de matrices para la ecuación de Hartree-Fock

Dada la definición de la funcion de onda, se escribe la ecuación de Hartree-Fock de forma matricial aplicando el bra-ket <phi_nu | phi_mu>, tal que

HF integrales....

En el miembro derecho es posible definir la matriz de superposición Sij como

intergral Sij

Así como la matriz diagonal donde se almacenan los valores de la energía. En cuanto al miembro derecho se pueden definir dos matrices. En primer lugar se encuentra la matriz Hcore calculada como

integral Hcore...

La cual hace referencia a los términos cinéticos y potenciales electrón-núcleo. Por otro lado la matriz para V_hf se redefine como matriz G y esta se calcula de la siguiente forma

integral G...

Defininiendo la matriz P de densidad de carga de la forma

definido P

Aprovechando los coeficientes que se obtienen al utilizar las funciones de onda base. Con estas matrices la ecuación queda escrita matricialmente como

FC = SCE

Donde F = H + G es la matriz de Fock. La matriz emergente C se calculan iterativamente diagonalizando la matriz F, y posterior a ello se recalcula matriz de densidad electronica P con el fin de que en cada iteración la matriz P converja a unos determinados valores y la matriz E se minimice a los valores más cercanos de energía. Cabe resaltar, que para mayor facilidad en el cálculo computacional se puede realizar las siguientes sustituciones

F' = S-1/2 F S-1/2

C' = S1/2 C

Para que la ecuación se escriba de forma equivalente como una ecuación de autovalores de la forma

F'C' = C' E'

Estas iteraciones se conocen como el método de campo autoconsistente y se explicará a continuación

### Método de campo autoconsistente (SCF)

Este es el algoritmo principal para hallar la energía del sistema, la cual se puede resumir en el siguiente diagrama de flujo.

Diagrama...

Asignando un P inicial arbitrario (que en este caso se asignó una matriz nula) Se realiza el cálculo de la matriz de Fock, para posteriormente diagonalizar F' y obtener los autovalores de energía E y los autovectores C'. Realizando las sustituciones adecuadas se recalcula la matriz P con los coeficientes obtenidos en C. Finalmente se establece un criterio de convergencia, que evalúa entre un P_old de la iteración anterior y el P actual la siguiente condición

P_old - P < epsilon 
>

Donde epsilon es un valor pequeño y se define en el computador, de acuerdo a la precisión del mismo.  Una vez se cumpla ese criterio de convergencia, se habrán obtenido los autovalores reales de energía de la matriz E.

### Perturbación Moller-Plesset

Luego de realizar el cálculo Hatree-Fock, es convieniente utilizar una corrección de esta energía, usando la teoría de perturbaciones independiente del tiempo Moller-Plesset. Esta consiste en tomar el potencial electrón-electrón V_hf como una perturbación, tal que se calcula a segundo orden

E MP2 = 1/4 SUM....

Cabe resaltar que este cálculo requiere reorganizar los valores encontrados y no es en sí un método iterativo, lo cual lo hace relativamente sencillo de implementar y mejora sustancialmente el cálculo de la energía.
## Estructura del código

Para llevar a cabo este calculo se implementó el método en una script escrita en C++ (MP2.cpp), la cual utiliza la librería Eigen 3, esencial para realizar cálculos matriciales [2]. Observando las funciones y métodos escritos en la script, es posible identificar las siguientes etapas.

En primer lugar se utilizan las funciones S_int y T_int para calcular las integrales necesarias en cada elemento matricial del Hamiltoniano nuclear Hcore y la matriz de superposición S, necesarias para calcular la matriz G y F. La función Intgrl realiza la integración y guarda los valores a sus matrices correspondientes por medio de la función Colect, así como los elementos de G que deben ser almacenados en un arreglo de 4 dimensiones, ya que según su ecuación requiere de varias integrales. Esta función a su vez permite realizar el cálculo de S1/2, la cual se denomina matriz X, y será necesaria para el método SCF.Este método utiliza los valores calculados en H y S para ejecutar el algoritmo de campo autoconsistente, tal como se explica en la sección anterior. Finalmente se procede a realizar el cálculo Moller-Plesset a las energías obtenidas por Hartree-Fock, utilizando un arreglo de 4 dimensiones para reorganizar los valores de energía obtenidos. 

Finalmente se varía la distancia R electrón-nucleo con el fin de observar el comportamiento de la energía en funcion de esta distancia, utilizando los métodos previamente mencionados. En las figuras 1 y 2 se puede observar la energía para H2 y HeH+ respectivamente.







 




 





