/*

Tenemos que hacer varias funciones para no perder precisión.

Reescalamiento: distancia, tiempo y masa

- Reescalamiento 1: 1 UA = 1
- Reescalamiento 2: 24 horas = 1 día
- Reescalamiento 3: 1 Masa Solar = 1


(por ejemplo)

Importante que luego tenemos que deshacer el reescalamiento para obtener los valores reales.


Hay que reescalar también las derivadas y las segundas derivadas

si r -> 1/c * r entonces dr -> 1/c * dr

*/


/*
Creamos un vector de N posiciones de planetas
Otro vector del tiempo, R
Otro vector de masas, N.

La energía podemos calcularla en cada paso de la simulación. La energía TIENE QUE CONSERVARSE. Si crece o decrece algo malo hemos hecho.

La energía potencial también tenemos que reescalarla. Como U=G*M*m/r, si reescalamos r -> 1/c * r, entonces U -> c * U, lo mismo con las masas.

*/



