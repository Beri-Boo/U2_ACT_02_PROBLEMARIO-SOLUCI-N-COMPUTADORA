# U2_ACT_02_PROBLEMARIO-SOLUCION-COMPUTADORA
# Calculadora de Criterios de Falla por Fatiga ⚙️

## Descripción
Este repositorio contiene el código en Python para automatizar los cálculos de diseño por fatiga bajo cargas combinadas (flexión y torsión). El script utiliza interpolación bidimensional para obtener de manera precisa los factores de concentración de esfuerzo (Kt, Kts) y las sensibilidades a la muesca (q, qc) a partir de gráficas estandarizadas. 

Posteriormente, calcula los esfuerzos máximos de Von Mises y determina el factor de diseño (n) utilizando los siguientes criterios:
* ED-Goodman modificado
* ED-Soderberg
* ED-Gerber
* ED-ASME elíptica

## Requisitos y Dependencias
Para ejecutar el código principal, es necesario contar con un entorno de Python 3.x y las siguientes librerías:
* `math` (Nativa)
* `numpy`
* `scipy`
* `matplotlib`

Puedes instalar las dependencias externas ejecutando el siguiente comando en tu terminal:
```bash
pip install numpy scipy matplotlib
