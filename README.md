# Modelado de COVID19

Modelado de la transmisión del COVID19 empleando el modelo epidemiológico SEIQRD, SEIQRDP y SEIRD. Luego de cierta experimentación con los tres modelos, se determinó que aquel que refleja de mejor manera las condiciones en las que se propaga el COVID19 era el modelo SEIQRD. Debido a esto, para el mismo (encontrado en la carpeta `Matlab/Modelo SEIQRD/`) se presentan dos tipos de simulaciones: Una estocástica (dependiente de factores aleatorios) y una determinística (sin aleatoriedad). 

## Simulación Determinista

Simulación basada en el modelado de la enfermedad como un sistema de ecuaciones diferenciales. Carece de elementos aleatorios, ya que la evolución del modelo depende únicamente de las condiciones iniciales que se le provean al mismo. La visualización del mismo es en la forma de un conjunto de "plots" que despliegan la evolución de los valores de susceptibles (S), expuestos (E), infectados (I), en cuarentenados (Q), recuperados (R) y muertos (D). A este modelo fue al que se le aplicaron estrategias de control, lo que permitió demostrar que al reducir el número de personas con el que un infectado tiene contacto, la curva de la enfermedad se reduce drásticamente.

## Simulación Estocástica

Simulación basada en el modelado de la enferemedad como un sistema dinámico de partículas (que representan personas) que se mueven con direcciones y velocidades aleatorias en una región cerrada. Esta simulación cuenta con muchos más parámetros, con la capacidad de cambiar factores como el número de partículas, la frecuencia de infección, el factor de riesgo asociado a diferentes grupos de edad, el rango de infección, la probabilidad de muerte, entre otros. La visualización de este modelo es por medio de una animación que despliega el movimiento en tiempo real de las partículas. Para facilitar el proceso de prueba y error, el script primero genera una versión simplificada de las animaciones que corren mucho más rápido. Si al usuario le gustaron los resultados obtenidos, puede optar por repetir la simulación empleando la misma "seed" para así generar una animación que será grabada en la forma de un video de forma mucho más presentable.

## Resultados

Se editó un video junto con Gabriela Iriarte presentando todos los resultados que se obtuvieron durante la realización de esta investigación. Hacer click en la imagen para ver el video. 

[![Controlando la Pandemia](https://img.youtube.com/vi/vihgAMi1RC0/0.jpg)](https://www.youtube.com/watch?v=vihgAMi1RC0 "Controlando la Pandemia")
