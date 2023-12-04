# ThermalParallelization
Proyecto de paralelización de código de difusión de calor. 
Trabajo final de la asignatura Sistemas de Cómputo Paralelo, ingeniería de computadores, facultad de informática, EHU/UPV.
Código base: Olatz Viñáspre (olatz.perezdevinaspre@ehu.eus) 
Asignatura SCP: Olatz Viñáspre (olatz.perezdevinaspre@ehu.eus) y Javier Navaridas (javier.navaridas@ehu.eus)

El problema tiene escala variable. La placa tiene un tamaño base de (100,200) multiplicado por la escala (max 12).  
En un fichero de entrada se determinan las posiciones, los tamaños y las temepraturas de los chips, además de otros parámetros.  
El código trabaja un conjunto de configuraciones ya dado, y realiza un cálculo iterativo de la temperatura media. En cada iteración, los chips inyectan calor según su temperatura, en sus posiciones. Por otro lado, la placa cuenta con refrigeración en el centro. El algoritmo converge una vez alcanzadas un número de iteraciones o cuando el cambio de temperatura media no sea significativo.

Proyecto llevado a cabo por Alejandro Perez (https://github.com/Aleeja137) y Alejandro Martin (https://github.com/almaqu)
# Tareas a realizar
En la primera fase, el objetivo era paralelizar el cálculo de la temperatura media, usando N procesos mediante MPI. Para ello, se divide la placa en conjuntos de filas, y cada procesador calcula un conjunto. Los procesadores con filas consecutivas deben comunicarse entre sí al final de cada 'barrida' (la actualización de la temperatura en cada celda), ya que la temperatura de una celda depende de las celdas adjacentes. Se pedían son enfoques para esto, usar Ssend en la primera versión, Isend después.

Para la segunda fase, se pedía crear grupos de trabajo, donde un solo proceso se encargaba de distribuír configuraciones a cada grupo, donde un proceso lo recibía y lo distribuía entre su grupo. Es un modelos de manager/worker. Después este mismo proceso devolvía el resultado y se quedaba a la espera de otra configuraciṕn que trabajar

# Código
En este repositorio está únicamente el código de la fase 2. Es posible encontrar todo el código [aquí](https://drive.google.com/drive/folders/1NHJrw9s4qAJeTJD-5bwrLU8Au7l_3QyH?usp=sharing), así como el informe final.

# Entorno
Para ejecutar el código paralelo se ha utilizado un cluster de 53 nodos (uno utilizado para el desarrollo
del código y los otros 52 para realizar su ejecución. Cada nodo tiene 2 procesadores Intel Xeon E5462
de 4 cores, 2,80 GHz, con 16GB RAM. Por lo tanto, en total disponemos de 416 núcleos.
Realmente durante el proceso solamente se han utilizado menos debido a que no todos estaban
encendidos.

# Ejecución
Para compilar el código se ha utilizado la herramienta mpicc, el compilador de C para programas escritos en MPI, el cual se encarga de compilar y vincular programas con las bibliotecas MPI necesarias.
Un ejemplo de uso de mpicc para compilar el proyecto:  
mpicc -o heat_p heat_p.c diffusion_p.c faux_p.c  

Para ejecutar el archivo binario del programa se ha utilizado un script denominado exepar, el cual obtiene como parámetros cuyo contenido es el siguiente:  
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n $1 $2 $3

# Informe de resultados
En general, con Ssend se ha alcanzado un speed-up de 15x usando 16 procesadores y con Isend del 25x usando 32 procesadores.
Usando el modelo manager/worker, se ha alcanzado un speed-up máximo de 63x con 64 procesadores con grupos de trabajo de 8 procesadores cada. Es muy importante determinar el tamaño de los grupos de trabajo, para más información consultar el [informe](https://drive.google.com/file/d/1yeiIJ7J-iFaMlf4B1YmluMrh9pPuLLYw/view?usp=sharing)
