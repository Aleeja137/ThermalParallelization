/* File: diffusion.h */ 

extern 
double calculate_Tmean (float *grid, float *grid_chips, float *grid_aux, float t_delta, int max_iter, float t_ext, int NROW_loc, int NROW_glob, int NCOL_glob, int pid, int npr, MPI_Comm com);

