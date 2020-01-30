# Proyecto Cuantica II: Calculo Hartree-Fock para una molecula diatomica aplicando teoria de perturbacion a segundo orden de Moller-Plesset

*Resumen*

En este proyecto se ha creado un programa que calcula la energía de correlación electrónica, para un sistema cuántico compuesto por dos núcleos y dos electrones, o bien una molécula diatómica. Este cálculo es bastante complicado de realizar de forma analítica, ya que tenemos potenciales de Coulomb electrón-electrón y electrón-núcleo los cuales complican la ecuación de Schrödinguer. Para este tipo de sistemas se recurre a un método variacional conocido como método Hartree-Fock, el cual consiste en tomar el potencial de Coulomb electrón-electrón como la interacción entre un sólo electrón y una distribución de carga provocada por los demás electrones en el sistema. Posterior a ello se elige unas funciones de onda iniciales e hidrogenoides, concebidas como determinantes de Slater para mantener sus propiedades antisimétricas, las cuales harán parte del cálculo de campo autoconsistente (SCF). Aunque con este método es posible obtener los orbitales, se considerará únicamente la energía de correlación electrónica para aplicar una perturbación de segundo orden de Moller-Plesset, la cual realizará una corrección de la energía dependiendo del sistema.


