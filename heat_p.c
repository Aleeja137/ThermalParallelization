/* heat_s.c 

     Difusion del calor en 2 dimensiones      Version en serie

     Se analizan las posiciones de los chips en una tarjeta, para conseguir la temperatura minima
     de la tarjeta. Se utiliza el metodo de tipo Poisson, y la tarjeta se discretiza en una rejilla 
     de puntos 2D.
     
     Entrada: card > la definicion de la tarjeta y las configuraciones a simular
     Salida: la mejor configuracion y la temperatura media
	      card_s.chips: situacion termica inicial
        card_s.res: la situacion termica final
 
     defines.h: definiciones de ciertas variables y estructuras de datos

     Compilar con estos dos ficheros: 
       diffusion.c: insertar calor, difundir y calcular la temperatura media hasta que se estabilice
       faux.c: ciertas funciones auxiliares

************************************************************************************************/

#include <stdio.h>
#include <values.h>
#include <time.h>
#include <mpi.h>

#include "defines.h"
#include "faux_p.h"
#include "diffusion_p.h"


/************************************************************************************/
void init_grid_chips (int conf, struct info_param param, struct info_chips *chips, int **chip_coord, float *grid_chips)
{
  int i, j, n;

  for (i=0; i<NROW; i++)
  for (j=0; j<NCOL; j++)  
    grid_chips[i*NCOL+j] = param.t_ext;

  
  for (n=0; n<param.nchip; n++)
  for (i = chip_coord[conf][2*n]   * param.scale; i < (chip_coord[conf][2*n] + chips[n].h) * param.scale; i++)
  for (j = chip_coord[conf][2*n+1] * param.scale; j < (chip_coord[conf][2*n+1]+chips[n].w) * param.scale; j++) 
    grid_chips[(i+1)*NCOL+(j+1)] = chips[n].tchip;
  
}

/************************************************************************************/
void init_grids (float t_ext, float *grid, float *grid_aux, int NROW_glob, int NCOL_glob)
{
  int i, j;

  for (i=0; i<NROW_glob; i++)
  for (j=0; j<NCOL_glob; j++) 
    grid[i*NCOL_glob+j] = grid_aux[i*NCOL_glob+j] = t_ext;
}

/************************************************************************************/
/************************************************************************************/
int main (int argc, char *argv[])
{
  struct info_param param;  // P0
  struct info_chips *chips; // P0
  int	 **chip_coord;        // P0
  struct info_results BT;   // P0
  double t0, t1, tej;       // P0

  float *grid, *grid_chips;
  int    conf, i;
  double Tmean = -1;


  // Variables nuevas
  int pid = -1, npr;
  int *tam, *tam_completo, *dis, *dis_completo;
  int nconf, max_iter, scale, pos, pack_read_data_size;
  float t_delta, t_ext;
  char *pack_read_data;
  float *grid_loc, *grid_chips_loc, *grid_aux_loc;  
  int NROW_glob, NCOL_glob;

  // Variables fase 2
  int tam_grupos_trabajo, *all_workers_pids, all_workers_pid = -1, npr_workers;
  int pid_worker = -1, proc_acabados, pack_results_size;
  int origen, tag, destino, tarea, n, acabado, nchip;

  char *pack_results;

  MPI_Group WORLD_GRP, ALL_WORKERS_GRP;
  MPI_Comm ALL_WORKERS_COMM, WORKER_COMM;
  MPI_Status info;


  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&pid);
  MPI_Comm_size(MPI_COMM_WORLD,&npr);

 // Comprobar que se ha dado un archivo con los parámetros y el tamaño de los grupos de trabajo
  if (argc != 3) {
    if (pid == 0) printf ("\n\nERROR: needs a card description file and size of workers group\n\n");
    MPI_Finalize();
    return (-1);
  } 

  // Comprobar número de procesadores
  tam_grupos_trabajo = atoi(argv[2]);
  if (npr % tam_grupos_trabajo != 1){
    if (pid == 0) printf ("\n\nERROR: Size of workers group must be multiple of npr plus one. Example: Size of workers group = 8; npr = 17\n\n");
    MPI_Finalize();
    return (-1);
  } 

  // Crear grupos de trabajo

    // Se obtiene el grupo original
    MPI_Comm_group(MPI_COMM_WORLD, &WORLD_GRP);

    // Se crea el comunicador de los procesadores trabajadores
    all_workers_pids = (int*) malloc ((npr-1)*sizeof(int));
    for ( i = 0; i < npr; i++) all_workers_pids[i] = i;
    MPI_Group_excl(WORLD_GRP,1,all_workers_pids,&ALL_WORKERS_GRP);
    MPI_Comm_create(MPI_COMM_WORLD,ALL_WORKERS_GRP,&ALL_WORKERS_COMM);

    // Se crean los grupos de trabajo (en total (npr-1)/tam_grupos_trabajo)
    if (pid != 0)
    {
      MPI_Comm_rank(ALL_WORKERS_COMM,&all_workers_pid);
      MPI_Comm_split(ALL_WORKERS_COMM,all_workers_pid%((npr-1)/tam_grupos_trabajo),all_workers_pid,&WORKER_COMM);
      MPI_Comm_size(WORKER_COMM,&npr_workers);
      MPI_Comm_rank(WORKER_COMM,&pid_worker);
    } 

  // P0 de cada trabajador de grupo lee los datos de card y manda datos necesarios a todos los trabajadores (P0 también los lee para tener NCOL y NROW para reserva de memoria)
  MPI_Pack_size(6,MPI_INT,MPI_COMM_WORLD,&pack_read_data_size);
  pack_read_data = (char*) malloc (pack_read_data_size*sizeof(char));
  pos = 0;

  if (pid_worker == 0) 
  {
    read_data (argv[1], &param, &chips, &chip_coord);
    MPI_Pack(&param.nconf,    1,MPI_INT,pack_read_data,pack_read_data_size,&pos,WORKER_COMM);
    MPI_Pack(&param.max_iter, 1,MPI_INT,pack_read_data,pack_read_data_size,&pos,WORKER_COMM);
    MPI_Pack(&param.t_ext,    1,MPI_INT,pack_read_data,pack_read_data_size,&pos,WORKER_COMM);
    MPI_Pack(&param.t_delta,  1,MPI_INT,pack_read_data,pack_read_data_size,&pos,WORKER_COMM);
    MPI_Pack(&param.scale,    1,MPI_INT,pack_read_data,pack_read_data_size,&pos,WORKER_COMM);
    MPI_Pack(&param.nchip,    1,MPI_INT,pack_read_data,pack_read_data_size,&pos,WORKER_COMM);
    nconf = param.nconf; max_iter = param.max_iter; t_ext = param.t_ext; t_delta = param.t_delta; scale = param.scale; nchip = param.nchip;
  } else if (pid == 0) read_data (argv[1], &param, &chips, &chip_coord); 
  
  if (pid != 0) MPI_Bcast(pack_read_data,pack_read_data_size,MPI_PACKED,0,WORKER_COMM);

  if (pid != 0 && pid_worker > 0)
  {
    MPI_Unpack(pack_read_data,pack_read_data_size,&pos,&nconf,1,MPI_INT,WORKER_COMM);
    MPI_Unpack(pack_read_data,pack_read_data_size,&pos,&max_iter,1,MPI_INT,WORKER_COMM);
    MPI_Unpack(pack_read_data,pack_read_data_size,&pos,&t_ext,1,MPI_INT,WORKER_COMM);
    MPI_Unpack(pack_read_data,pack_read_data_size,&pos,&t_delta,1,MPI_INT,WORKER_COMM);
    MPI_Unpack(pack_read_data,pack_read_data_size,&pos,&scale,1,MPI_INT,WORKER_COMM);
    MPI_Unpack(pack_read_data,pack_read_data_size,&pos,&nchip,1,MPI_INT,WORKER_COMM);
  }

  // No es código elegante pero a estas alturas resuelve el problema
  if (pid == 0) 
  {
    NROW_glob = (200*param.scale) + 2;
    NCOL_glob = (100*param.scale) + 2;
  }
  else
  {
    NROW_glob = (200*scale)+2; 
    NCOL_glob = (100*scale)+2;
  }

  // Reservar memoria para enviar y recibir los resultados de cada grupo de workers hacia P0
  // Conf (1 INT), Tmean (1 DOUBLE), Grid (NROW*NCOL*FLOAT), Grid_chips (NROW*NCOL*FLOAT)
  if (pid == 0 || pid_worker == 0)
  {
    pack_results_size = (4+8+(2*NROW_glob*NCOL_glob*4));
    pack_results = (char*) malloc (pack_results_size*sizeof(char));
  }

  if (pid != 0) // Los procesos workers calculan el tamaño de sus grid y grid_chips locales
  {
    // Cada proceso calcula cuántas líneas necesita en tam y dis
    int cociente = (NROW_glob-2) / npr_workers;
    int resto    = (NROW_glob-2) % npr_workers;

    tam          = (int*) malloc (npr_workers*sizeof(int));
    tam_completo = (int*) malloc (npr_workers*sizeof(int));
    dis          = (int*) malloc (npr_workers*sizeof(int));
    dis_completo = (int*) malloc (npr_workers*sizeof(int));

    for (i=0;i<npr_workers;i++)
    {
      // Cuántas filas trabaja cada procesador
      // tam_completo tiene las líneas de los vecinos añadidas, de utiliza para el scatter y gather
      tam[i] = cociente;
      if (i<resto) tam[i]++;
      tam_completo[i] = tam[i] + 2;

      // Se multiplica por los elementos que trabaja
      tam[i]          *= NCOL_glob;
      tam_completo[i] *= NCOL_glob;

      // A qué distancia están desde la fila de trabajo inicial
      if (i==0) dis[i] = NCOL_glob;
      else dis[i] = tam[i-1] + dis[i-1];  

      // A qué distancia están desde la fila inicial
      dis_completo[i] = dis[i] - NCOL_glob;
    }
  }
  
  if (pid == 0) // Solo P0 maneja los resultados 
  {
    BT.bgrid = malloc(NROW_glob * NCOL_glob * sizeof(float));
    BT.cgrid = malloc(NROW_glob * NCOL_glob * sizeof(float));
    BT.Tmean = MAXDOUBLE;

    grid       = malloc(NROW_glob * NCOL_glob * sizeof(float));
    grid_chips = malloc(NROW_glob * NCOL_glob * sizeof(float));
  } 
  else 
  {
    if (pid_worker == 0) // Grids completos para recoger los valores 
    {
      grid       = malloc(NROW_glob * NCOL_glob * sizeof(float));
      grid_chips = malloc(NROW_glob * NCOL_glob * sizeof(float));
    } 
    
    // Cada proceso inicializa grid_loc, grid_chips_loc y grid_aux_loc
    grid_loc       = (float*) malloc(tam_completo[pid_worker] * sizeof(float));
    grid_chips_loc  = (float*) malloc(tam_completo[pid_worker] * sizeof(float));
    grid_aux_loc   = (float*) malloc(tam_completo[pid_worker] * sizeof(float));
  }   

  if (pid == 0)
  {
    printf ("\n  ===================================================================");
    printf ("\n    Thermal diffusion - SERIAL version ");
    printf ("\n    %d x %d points, %d chips", RSIZE*param.scale, CSIZE*param.scale, param.nchip);
    printf ("\n    T_ext = %1.1f, Tmax_chip = %1.1f, T_delta: %1.3f, Max_iter: %d", param.t_ext, param.tmax_chip, param.t_delta, param.max_iter);
    printf ("\n  ===================================================================\n\n");
    t0 = MPI_Wtime(); // Solo P0, cambiar por Wtime()
  } 
  
  // Loop to process chip configurations
  if (pid == 0)
  {
    proc_acabados = 0; conf = 0;
    while(proc_acabados < npr-1)
    {
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG,MPI_COMM_WORLD,&info);
      origen = info.MPI_SOURCE;
      tag = info.MPI_TAG;

      // Petición inicial, recibir mensaje vacío
      if (tag == 3) MPI_Recv(&tarea,0,MPI_INT,origen,tag,MPI_COMM_WORLD,&info); 
      else if (tag == 4) // Respuesta de resultados, recibir y procesar los resultados
      {
        // Recibir resultados
        MPI_Recv(pack_results,pack_results_size,MPI_PACKED,origen,tag,MPI_COMM_WORLD,&info);
        // Desempaquetar resultado
        pos = 0;
        int conf_recv;
        MPI_Unpack(pack_results,pack_results_size,&pos,&conf_recv,1,MPI_INT,MPI_COMM_WORLD);
        MPI_Unpack(pack_results,pack_results_size,&pos,&Tmean,1,MPI_DOUBLE,MPI_COMM_WORLD);
        MPI_Unpack(pack_results,pack_results_size,&pos,grid,NROW_glob*NCOL_glob,MPI_FLOAT,MPI_COMM_WORLD);
        MPI_Unpack(pack_results,pack_results_size,&pos,grid_chips,NROW_glob*NCOL_glob,MPI_FLOAT,MPI_COMM_WORLD);
        // Processing configuration results 
        results_conf (conf_recv, Tmean, param, grid, grid_chips, &BT);
      }

      if (conf < param.nconf) 
      {
        tag = 1; // Hay tarea
        destino = origen;
        MPI_Send(&conf,1,MPI_INT,origen,tag,MPI_COMM_WORLD);
        conf++;
      } 
      else
      {
        tag = 2; // No hay tarea
        destino = origen;
        MPI_Send(&conf,0,MPI_INT,origen,tag,MPI_COMM_WORLD);
        proc_acabados += (npr-1) / ((npr-1)/tam_grupos_trabajo);
      }
    }
  } 
  else
  {
    conf = -1; destino = 0;
    if (pid_worker == 0 && conf == -1) // P0 de cada grupo de trabajo manda la primera petición
    {
      tag = 3;
      MPI_Send(&tarea,0,MPI_INT,destino,tag,MPI_COMM_WORLD);
      conf = 0;
    }

    while(1)
    {
      if (pid_worker == 0)
      {
        MPI_Probe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&info);
        origen = info.MPI_SOURCE;
        tag = info.MPI_TAG;
        if (tag == 1) // Hay tarea
        {
          MPI_Recv(&conf,1,MPI_INT,origen,tag,MPI_COMM_WORLD,&info);
        }
        else if (tag == 2) // No hay tarea
        { 
          MPI_Recv(&conf,0,MPI_INT,origen,tag,MPI_COMM_WORLD,&info); 
          conf = -1;
        }
      }

      // Comprobar que hay tarea, terminar sino
      MPI_Bcast(&conf,1,MPI_INT,0,WORKER_COMM);
      if (conf == -1) break;

      // Trabajar la configuración
        // P0 lee grid_chips entero y lo distribuye
        if (pid_worker == 0) init_grid_chips (conf, param, chips, chip_coord, grid_chips);
        MPI_Scatterv(&grid_chips[NCOL_glob],tam,dis_completo,MPI_FLOAT,&grid_chips_loc[NCOL_glob],tam[pid_worker],MPI_FLOAT,0,WORKER_COMM);
        
        // Inintial values for grids
        init_grids (t_ext, grid_loc, grid_aux_loc, tam_completo[pid_worker]/NCOL_glob, NCOL_glob);

        // Main loop: thermal injection/disipation until convergence (t_delta or max_iter)
        Tmean = calculate_Tmean (grid_loc, grid_chips_loc, grid_aux_loc, t_delta, max_iter, t_ext, tam[pid_worker]/NCOL_glob, NROW_glob, NCOL_glob, pid_worker, npr_workers, WORKER_COMM);
        
        MPI_Gatherv(&grid_loc[NCOL_glob],tam[pid_worker],MPI_FLOAT,grid,tam,dis,MPI_FLOAT,0,WORKER_COMM);
        
      // Preparar y mandar resultado
      if (pid_worker == 0)
      {
        tag = 4; pos = 0;

        MPI_Pack(&conf,1,MPI_INT,pack_results,pack_results_size,&pos,MPI_COMM_WORLD);
        MPI_Pack(&Tmean,1,MPI_DOUBLE,pack_results,pack_results_size,&pos,MPI_COMM_WORLD);
        MPI_Pack(grid,NROW_glob*NCOL_glob,MPI_FLOAT,pack_results,pack_results_size,&pos,MPI_COMM_WORLD);
        MPI_Pack(grid_chips,NROW_glob*NCOL_glob,MPI_FLOAT,pack_results,pack_results_size,&pos,MPI_COMM_WORLD);

        printf ("  Config: %2d    Tmean: %1.2f\n", conf + 1, Tmean);
        
        MPI_Send(pack_results,pack_results_size,MPI_PACKED,origen,tag,MPI_COMM_WORLD);
      }
        
    }
  }
    

  
  
  if (pid == 0) 
  {
    t1 = MPI_Wtime(); // Cambiar por Wtime()
    tej = t1-t0;
    printf ("\n\n >>> Best configuration: %2d    Tmean: %1.2f\n", BT.conf + 1, BT.Tmean);
    printf ("   > Time (paralel): %1.3f s \n\n", tej);
    // Writing best configuration results
    results (param, &BT, argv[1]);
  }

  if (pid == 0 || pid_worker == 0)
  {
    free (chips); free(grid); free (grid_chips);
    for (i=0; i<param.nconf; i++) free (chip_coord[i]);
    free (chip_coord); free (pack_results);
    
  }
  if (pid == 0)
  {
    free (BT.bgrid); free (BT.cgrid);
  } 
  else
  {
    free(tam); free(dis); free(tam_completo); free(dis_completo); free(pack_read_data);
    free(grid_loc); free(grid_chips_loc); free(grid_aux_loc);
    MPI_Comm_free(&WORKER_COMM); MPI_Comm_free(&ALL_WORKERS_COMM);
  }

  free(all_workers_pids);
  MPI_Group_free(&WORLD_GRP); MPI_Group_free(&ALL_WORKERS_GRP);
   

  MPI_Finalize();
  return (0);
}

