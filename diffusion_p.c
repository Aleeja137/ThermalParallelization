/* File: diffusion_s.c */ 

#include "defines.h"
#include <mpi.h>

/************************************************************************************/
void thermal_update (float *grid, float *grid_chips, float t_ext, int NROW_loc, int NCOL_glob, int pid, int npr, MPI_Request *req_send_up, MPI_Request *req_send_down, MPI_Comm com)
{
  int i, j, a, b;
  
  // Heat injection at chip positions
  for (i=1; i<=NROW_loc; i++)
  for (j=1; j<NCOL_glob-1; j++) 
    if (grid_chips[i*NCOL_glob+j] > grid[i*NCOL_glob+j])
      grid[i*NCOL_glob+j] += 0.05 * (grid_chips[i*NCOL_glob+j] - grid[i*NCOL_glob+j]);
      
  // Air cooling at the middle of the card
  a = 0.44*(NCOL_glob-2) + 1;
  b = 0.56*(NCOL_glob-2) + 1;
  
  for (i=1; i<=NROW_loc; i++)
    for (j=a; j<b; j++) grid[i*NCOL_glob+j] -= 0.01 * (grid[i*NCOL_glob+j] - t_ext);

  if (pid !=0)        MPI_Isend(&grid[NCOL_glob],         NCOL_glob,MPI_FLOAT,pid-1,0,com,req_send_up);
  if (pid != (npr-1)) MPI_Isend(&grid[NROW_loc*NCOL_glob],NCOL_glob,MPI_FLOAT,pid+1,0,com,req_send_down);
}

/************************************************************************************/
double thermal_diffusion (float *grid, float *grid_aux, int NROW_loc, int NCOL_glob, int pid, int npr, MPI_Request *req_recv_up, MPI_Request *req_recv_down, MPI_Request *req_send_up, MPI_Request *req_send_down)
{
  int    i, j;
  double  T;
  double Tfull = 0.0;

  // Se procesa el bloque 'interno' del grid
  for (i=2; i<=NROW_loc-1; i++)
    for (j=1; j<NCOL_glob-1; j++)
    {
      T = grid[i*NCOL_glob+j] + 
          0.10 * (grid[(i+1)*NCOL_glob+j]   + grid[(i-1)*NCOL_glob+j]   + grid[i*NCOL_glob+(j+1)]     + grid[i*NCOL_glob+(j-1)] + 
                  grid[(i+1)*NCOL_glob+j+1] + grid[(i-1)*NCOL_glob+j+1] + grid[(i+1)*NCOL_glob+(j-1)] + grid[(i-1)*NCOL_glob+(j-1)] 
                  - 8*grid[i*NCOL_glob+j]);

      grid_aux[i*NCOL_glob+j] = T;
      Tfull += T;
    }
    
  // Se procesa la primera línea del grid
  MPI_Wait(req_recv_up  , MPI_STATUS_IGNORE);
  
  for (j=1; j<NCOL_glob-1; j++)
  {
    T = grid[1*NCOL_glob+j] + 
        0.10 * (grid[(1+1)*NCOL_glob+j]   + grid[(1-1)*NCOL_glob+j]   + grid[1*NCOL_glob+(j+1)]     + grid[1*NCOL_glob+(j-1)] + 
                grid[(1+1)*NCOL_glob+j+1] + grid[(1-1)*NCOL_glob+j+1] + grid[(1+1)*NCOL_glob+(j-1)] + grid[(1-1)*NCOL_glob+(j-1)] 
                - 8*grid[1*NCOL_glob+j]);

    grid_aux[1*NCOL_glob+j] = T;
    Tfull += T;
  }
  
  // Se procesa la última línea del grid
  MPI_Wait(req_recv_down, MPI_STATUS_IGNORE);
  
  for (j=1; j<NCOL_glob-1; j++)
  {
    T = grid[NROW_loc*NCOL_glob+j] + 
        0.10 * (grid[(NROW_loc+1)*NCOL_glob+j]   + grid[(NROW_loc-1)*NCOL_glob+j]   + grid[NROW_loc*NCOL_glob+(j+1)]     + grid[NROW_loc*NCOL_glob+(j-1)] + 
                grid[(NROW_loc+1)*NCOL_glob+j+1] + grid[(NROW_loc-1)*NCOL_glob+j+1] + grid[(NROW_loc+1)*NCOL_glob+(j-1)] + grid[(NROW_loc-1)*NCOL_glob+(j-1)] 
                - 8*grid[NROW_loc*NCOL_glob+j]);

    grid_aux[NROW_loc*NCOL_glob+j] = T;
    Tfull += T;
  }

  // Esperar a que se hayan mandado las filas del grid arriba y abajo antes de modificarlas
  MPI_Wait(req_send_down,MPI_STATUS_IGNORE);
  MPI_Wait(req_send_up,  MPI_STATUS_IGNORE);

  // // New values for the grid
  // for (i=1; i<=NROW_loc; i++)
  // for (j=1; j<NCOL_glob-1; j++)
  //   grid[i*NCOL_glob+j] = grid_aux[i*NCOL_glob+j]; 

  return (Tfull);
}

/************************************************************************************/
double calculate_Tmean (float *grid, float *grid_chips, float *grid_aux, float t_delta, int max_iter, float t_ext, int NROW_loc, int NROW_glob, int NCOL_glob, int pid, int npr, MPI_Comm com)
{
  int    i, j, end, niter;
  double  Tfull, Tfull_loc;
  double Tmean, Tmean0 = t_ext;
  niter = 0; end = 0;

  MPI_Request req_send_up = MPI_REQUEST_NULL, req_send_down = MPI_REQUEST_NULL, req_recv_up = MPI_REQUEST_NULL, req_recv_down = MPI_REQUEST_NULL;

  while (end == 0)
  {
    niter++;
    Tmean = 0.0;
    
    // Heat injection and air cooling 
    if (niter % 2 == 0) thermal_update (grid    , grid_chips, t_ext, NROW_loc, NCOL_glob, pid, npr, &req_send_up, &req_send_down, com);
    else                thermal_update (grid_aux, grid_chips, t_ext, NROW_loc, NCOL_glob, pid, npr, &req_send_up, &req_send_down, com);
    
    if (pid != 0)
    {
      if (niter % 2 == 0) MPI_Irecv(&grid[0],    NCOL_glob,MPI_FLOAT,pid-1,0,com,&req_recv_up);
      else                MPI_Irecv(&grid_aux[0],NCOL_glob,MPI_FLOAT,pid-1,0,com,&req_recv_up);
    }       
    if (pid != (npr-1))
    {
      if (niter % 2 == 0) MPI_Irecv(&grid[(NROW_loc + 1) * NCOL_glob],    NCOL_glob,MPI_FLOAT,pid+1,0,com,&req_recv_down);
      else                MPI_Irecv(&grid_aux[(NROW_loc + 1) * NCOL_glob],NCOL_glob,MPI_FLOAT,pid+1,0,com,&req_recv_down);
    } 
    
    // Thermal diffusion
    if (niter % 2 == 0) Tfull_loc = thermal_diffusion(grid, grid_aux, NROW_loc, NCOL_glob, pid, npr, &req_recv_up, &req_recv_down, &req_send_up, &req_send_down);
    else                Tfull_loc = thermal_diffusion(grid_aux, grid, NROW_loc, NCOL_glob, pid, npr, &req_recv_up, &req_recv_down, &req_send_up, &req_send_down);
     
    
    // Convergence every 10 iterations
    if (niter % 10 == 0)
    {
      for (i=1; i<=NROW_loc; i++)
      for (j=1; j<NCOL_glob-1; j++)
        grid[i*NCOL_glob+j] = grid_aux[i*NCOL_glob+j]; 

      // MPI_Allreduce de Tfull_loc en Tfull, cuidado con los tipos (float vs double)
      MPI_Allreduce(&Tfull_loc,&Tfull,1,MPI_DOUBLE,MPI_SUM,com);

      Tmean = Tfull / ((NCOL_glob-2)*(NROW_glob-2));
      if ((fabs(Tmean - Tmean0) < t_delta) || (niter > max_iter)) end = 1;
      else Tmean0 = Tmean;
    }
  } // End while 
  if (pid == 0) printf ("Iter (par): %d\t", niter);
  
  if (req_recv_up != MPI_REQUEST_NULL)   MPI_Request_free(&req_recv_up);
  if (req_send_up != MPI_REQUEST_NULL)   MPI_Request_free(&req_send_up);
  if (req_recv_down != MPI_REQUEST_NULL) MPI_Request_free(&req_recv_down);
  if (req_send_down != MPI_REQUEST_NULL) MPI_Request_free(&req_send_down);

  return (Tmean);
}

