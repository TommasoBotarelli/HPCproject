#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#if defined(_OPENMP)
#include <omp.h>
#define GET_TIME() (omp_get_wtime()) // wall time
#else
#define GET_TIME() ((double)clock() / CLOCKS_PER_SEC) // cpu time
#endif

struct parameters {
  double dx, dy, dt, max_t;
  double g, gamma;
  int source_type;
  int sampling_rate;
  char input_h_filename[256];
  char output_eta_filename[256];
  char output_u_filename[256];
  char output_v_filename[256];
};

struct data {
  int nx, ny;
  double dx, dy;
  double *values;
};

struct data_message {
  int message_size;
  double *values;
};

#define GET(data, i, j) ((data)->values[(data)->nx * (j) + (i)])
#define SET(data, i, j, val) ((data)->values[(data)->nx * (j) + (i)] = (val))

#define GET_X_COORD(data, i) ((data)->dx * (i))
#define GET_Y_COORD(data, i) ((data)->dy * (i))

#define ERROR_VALUE -1
#define MAX_ITERATION 10

struct time_records {
  int counter;
  double times[MAX_ITERATION];
};

int read_parameters(struct parameters *param, const char *filename)
{
  FILE *fp = fopen(filename, "r");
  if(!fp) {
    printf("Error: Could not open parameter file '%s'\n", filename);
    return 1;
  }
  int ok = 1;
  if(ok) ok = (fscanf(fp, "%lf", &param->dx) == 1);
  if(ok) ok = (fscanf(fp, "%lf", &param->dy) == 1);
  if(ok) ok = (fscanf(fp, "%lf", &param->dt) == 1);
  if(ok) ok = (fscanf(fp, "%lf", &param->max_t) == 1);
  if(ok) ok = (fscanf(fp, "%lf", &param->g) == 1);
  if(ok) ok = (fscanf(fp, "%lf", &param->gamma) == 1);
  if(ok) ok = (fscanf(fp, "%d", &param->source_type) == 1);
  if(ok) ok = (fscanf(fp, "%d", &param->sampling_rate) == 1);
  if(ok) ok = (fscanf(fp, "%256s", param->input_h_filename) == 1);
  if(ok) ok = (fscanf(fp, "%256s", param->output_eta_filename) == 1);
  if(ok) ok = (fscanf(fp, "%256s", param->output_u_filename) == 1);
  if(ok) ok = (fscanf(fp, "%256s", param->output_v_filename) == 1);
  fclose(fp);
  if(!ok) {
    printf("Error: Could not read one or more parameters in '%s'\n", filename);
    return 1;
  }
  return 0;
}

void print_parameters(const struct parameters *param)
{
  printf("Parameters:\n");
  printf(" - grid spacing (dx, dy): %g m, %g m\n", param->dx, param->dy);
  printf(" - time step (dt): %g s\n", param->dt);
  printf(" - maximum time (max_t): %g s\n", param->max_t);
  printf(" - gravitational acceleration (g): %g m/s^2\n", param->g);
  printf(" - dissipation coefficient (gamma): %g 1/s\n", param->gamma);
  printf(" - source type: %d\n", param->source_type);
  printf(" - sampling rate: %d\n", param->sampling_rate);
  printf(" - input bathymetry (h) file: '%s'\n", param->input_h_filename);
  printf(" - output elevation (eta) file: '%s'\n", param->output_eta_filename);
  printf(" - output velocity (u, v) files: '%s', '%s'\n",
         param->output_u_filename, param->output_v_filename);
}

int read_data(struct data *data, const char *filename)
{
  FILE *fp = fopen(filename, "rb");
  if(!fp) {
    printf("Error: Could not open input data file '%s'\n", filename);
    return 1;
  }
  int ok = 1;
  if(ok) ok = (fread(&data->nx, sizeof(int), 1, fp) == 1);
  if(ok) ok = (fread(&data->ny, sizeof(int), 1, fp) == 1);
  if(ok) ok = (fread(&data->dx, sizeof(double), 1, fp) == 1);
  if(ok) ok = (fread(&data->dy, sizeof(double), 1, fp) == 1);
  if(ok) {
    int N = data->nx * data->ny;
    if(N <= 0) {
      printf("Error: Invalid number of data points %d\n", N);
      ok = 0;
    }
    else {
      data->values = (double*)malloc(N * sizeof(double));
      if(!data->values) {
        printf("Error: Could not allocate data (%d doubles)\n", N);
        ok = 0;
      }
      else {
        ok = (fread(data->values, sizeof(double), N, fp) == N);
      }
    }
  }
  fclose(fp);
  if(!ok) {
    printf("Error reading input data file '%s'\n", filename);
    return 1;
  }
  return 0;
}

int write_data(const struct data *data, const char *filename, int step)
{
  char out[512];
  if(step < 0)
    sprintf(out, "%s.dat", filename);
  else
    sprintf(out, "%s_%d.dat", filename, step);
  FILE *fp = fopen(out, "wb");
  if(!fp) {
    printf("Error: Could not open output data file '%s'\n", out);
    return 1;
  }
  int ok = 1;
  if(ok) ok = (fwrite(&data->nx, sizeof(int), 1, fp) == 1);
  if(ok) ok = (fwrite(&data->ny, sizeof(int), 1, fp) == 1);
  if(ok) ok = (fwrite(&data->dx, sizeof(double), 1, fp) == 1);
  if(ok) ok = (fwrite(&data->dy, sizeof(double), 1, fp) == 1);
  int N = data->nx * data->ny;
  if(ok) ok = (fwrite(data->values, sizeof(double), N, fp) == N);
  fclose(fp);
  if(!ok) {
    printf("Error writing data file '%s'\n", out);
    return 1;
  }
  return 0;
}

int write_data_vtk(const struct data *data, const char *name,
                   const char *filename, int step)
{
  char out[512];
  if(step < 0)
    sprintf(out, "%s.vti", filename);
  else
    sprintf(out, "%s_%d.vti", filename, step);

  FILE *fp = fopen(out, "wb");
  if(!fp) {
    printf("Error: Could not open output VTK file '%s'\n", out);
    return 1;
  }

  unsigned long num_points = data->nx * data->ny;
  unsigned long num_bytes = num_points * sizeof(double);

  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" "
          "byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp, "  <ImageData WholeExtent=\"0 %d 0 %d 0 0\" "
          "Spacing=\"%lf %lf 0.0\">\n",
          data->nx - 1, data->ny - 1, data->dx, data->dy);
  fprintf(fp, "    <Piece Extent=\"0 %d 0 %d 0 0\">\n",
          data->nx - 1, data->ny - 1);

  fprintf(fp, "      <PointData Scalars=\"scalar_data\">\n");
  fprintf(fp, "        <DataArray type=\"Float64\" Name=\"%s\" "
          "format=\"appended\" offset=\"0\">\n", name);
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </PointData>\n");

  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </ImageData>\n");

  fprintf(fp, "  <AppendedData encoding=\"raw\">\n_");

  fwrite(&num_bytes, sizeof(unsigned long), 1, fp);
  fwrite(data->values, sizeof(double), num_points, fp);

  fprintf(fp, "  </AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");

  fclose(fp);
  return 0;
}

int write_manifest_vtk(const char *name, const char *filename,
                       double dt, int nt, int sampling_rate)
{
  char out[512];
  sprintf(out, "%s.pvd", filename);

  FILE *fp = fopen(out, "wb");
  if(!fp) {
    printf("Error: Could not open output VTK manifest file '%s'\n", out);
    return 1;
  }

  fprintf(fp, "<VTKFile type=\"Collection\" version=\"0.1\" "
          "byte_order=\"LittleEndian\">\n");
  fprintf(fp, "  <Collection>\n");
  for(int n = 0; n < nt; n++) {
    if(sampling_rate && !(n % sampling_rate)) {
      double t = n * dt;
      fprintf(fp, "    <DataSet timestep=\"%g\" file='%s_%d.vti'/>\n", t,
              filename, n);
    }
  }
  fprintf(fp, "  </Collection>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
  return 0;
}

int init_data(struct data *data, int nx, int ny, double dx, double dy,
              double val)
{
  data->nx = nx;
  data->ny = ny;
  data->dx = dx;
  data->dy = dy;
  data->values = (double*)malloc(nx * ny * sizeof(double));
  if(!data->values){
    printf("Error: Could not allocate data\n");
    return 1;
  }
  for(int i = 0; i < nx * ny; i++) data->values[i] = val;
  return 0;
}

void free_data(struct data *data)
{
  free(data->values);
}

void print_to_file(const struct time_records *times, const char *filename){
  FILE *file = fopen(filename, "w");

  if (file == NULL){
      printf("Error opening file!\n");
      exit(1);
  }

  double total_time = 0.0;

  for (int i = 0; i < times->counter; i++){
    fprintf(file, "%f\n", times->times[i]);
    total_time += times->times[i];
  }

  double mean_time = total_time/(double)times->counter;

  fprintf(file, "Total time: %f\n", total_time);
  fprintf(file, "Mean time for execution: %f\n", mean_time);
}

double interpolate_data(const struct data *data, double x, double y)
{
  // TODO: this returns the nearest neighbor, should implement actual
  // interpolation instead
  double real_i = x / data->dx;
  double real_j = y / data->dy;

  int i = floor(real_i);
  int j = floor(real_j);

  if(i < 0) i = 0;
  else if(i > data->nx - 1) i = data->nx - 1;
  if(j < 0) j = 0;
  else if(j > data->ny - 1) j = data->ny - 1;

  double val;

  // TODO: rivedere i controlli se sono necessari o meno
  
  if (i >= data->nx-1 && j >= data->ny-1)
  {
    val = GET(data, i, j);
  }

  else if ((0 <= i && i < data->nx-1) && (0 <= j && j < data->ny-1)){
    val = (GET(data, i, j) * (GET_X_COORD(data, i+1) - x) * (GET_Y_COORD(data, j+1) - y)
                + GET(data, i+1, j) * (x - GET_X_COORD(data, i)) * (GET_Y_COORD(data, j+1) - y)
                + GET(data, i, j+1) * (GET_X_COORD(data, i+1) - x) * (y - GET_Y_COORD(data, j))
                + GET(data, i+1, j+1) * (x - GET_X_COORD(data, i)) * (y - GET_Y_COORD(data, j))) / (data->dx * data->dy);
  }

  else if (i == data->nx-1){
    val = (GET(data, i, j) * (GET_Y_COORD(data, j+1) - y) +
          GET(data, i, j+1) * (y - GET_Y_COORD(data, j))) / data->dy;
  }

  else if (j == data->ny-1){
    val = (GET(data, i, j) * (GET_X_COORD(data, i+1) - x) +
          GET(data, i+1, j) * (x - GET_X_COORD(data, i))) / data->dx;
  }

  else {
    val = ERROR_VALUE;
  }

  return val;
}

void save_coordinate(struct data *data, const char *filename){
    FILE* f = fopen(filename, "w");

    fprintf(f, "nx: %d\n", data->nx);
    fprintf(f, "ny: %d\n", data->ny);
    fprintf(f, "dx: %f\n", data->dx);
    fprintf(f, "dx: %f\n", data->dy);

    for(int j = 0; j < data->ny; j++) {
        for(int i = 0; i < data->nx; i++) {
            fprintf(f, "%f\n", GET(data, i, j));
        }
    }

    fclose(f);
}

void add_values_by_line(struct data* line_data, 
                        struct data* values_to_add, 
                        int size_i, 
                        int size_j, 
                        int dims[],
                        int coord){
  for (int j = 0; j < size_j; j++){
    for (int i = 0; i < size_i; i++){
      line_data->values[j * dims[1] + i] = GET(values_to_add, i, j);
    }
  }
}

int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  
  int world_size;
  int rank, cart_rank;

  int dims[2]    = {0, 0};
  int periods[2] = {0, 0};

  int reorder = 0;

  int coords[2];

  MPI_Comm cart_comm;

  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Dims_create(world_size, 2, dims);

  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cart_comm);
  MPI_Comm_rank(cart_comm, &cart_rank);

  MPI_Cart_coords(cart_comm, cart_rank, 2, coords);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  printf("Hello world from node %s, I'm rank %d out of %d ranks\n",
         processor_name, rank, world_size);

  struct data eta, u, v;
  // interpolate bathymetry
  struct data h_interp_u;
  struct data h_interp_v;

  if (rank == 0){
    if(argc != 2) {
      printf("Usage: %s parameter_file\n", argv[0]);
      return 1;
    }
  }

  struct parameters param;
  if(read_parameters(&param, argv[1])) return 1;
  if (rank == 0){
    print_parameters(&param);
  }

  struct data h;
  if(read_data(&h, param.input_h_filename)) return 1;

  // infer size of domain from input elevation data
  double hx = h.nx * h.dx;
  double hy = h.ny * h.dy;
  int nx = floor(hx / param.dx);
  int ny = floor(hy / param.dy);
  if(nx <= 0) nx = 1;
  if(ny <= 0) ny = 1;
  int nt = floor(param.max_t / param.dt);

  if (rank == 0){
    printf(" - grid size: %g m x %g m (%d x %d = %d grid points)\n",
          hx, hy, nx, ny, nx * ny);
    printf(" - number of time steps: %d\n", nt);
  }

  unsigned int start_i = (nx * coords[1] / dims[1]) + 1;
  unsigned int   end_i = nx * (coords[1]+1) / dims[1];

  unsigned int start_j = (ny * coords[0] / dims[0]) + 1;
  unsigned int   end_j = ny * (coords[0]+1) / dims[0];

  unsigned int mysize_i = end_i - start_i;
  unsigned int mysize_j = end_j - start_j;

  init_data(&eta, mysize_i, mysize_j, param.dx, param.dy, 0.);
  init_data(&u, mysize_i+1, mysize_j, param.dx, param.dy, 0.);
  init_data(&v, mysize_i, mysize_j+1, param.dx, param.dy, 0.);

  init_data(&h_interp_u, mysize_i + 1, mysize_j, param.dx, param.dy, 0.);
  for(int j = start_j; j < end_j; j++) {
    for(int i = start_i; i < end_i + 1; i++) {
      double x = i * param.dx;
      double y = ((double)j + 0.5) * param.dy;
      double val = interpolate_data(&h, x, y);
      
      SET(&h_interp_u, i - start_i, j - start_j, val);
    }
  }

  init_data(&h_interp_v, mysize_i, mysize_j + 1, param.dx, param.dy, 0.);
  for(int j = start_j; j < end_j + 1; j++) {
    for(int i = start_i; i < end_i; i++) {
      double x = ((double)i + 0.5) * param.dx;
      double y = j * param.dy;
      double val = interpolate_data(&h, x, y);
      
      SET(&h_interp_v, i - start_i, j - start_j, val);
    }
  }

  double start = GET_TIME();

  for(int n = 0; n < nt; n++) {

    if(n && (n % (nt / 10)) == 0 && rank == 0) {
      double time_sofar = GET_TIME() - start;
      double eta = (nt - n) * time_sofar / n;
      printf("Computing step %d/%d (ETA: %g seconds)     \r", n, nt, eta);
      fflush(stdout);
    }

    // output solution
    if(param.sampling_rate && !(n % param.sampling_rate)) {
      /* TODO: vedere se togliere questa parte
      if (coords[1] == 0){
        
        struct data line_data;
        init_data(&line_data, nx, mysize_j, param.dx, param.dy, 0.);
        //line_data.values = malloc(sizeof(double) * nx * mysize_j);

        add_values_by_line(&line_data, &eta, mysize_i, mysize_j, dims);
        
        double* sender_values = malloc(sizeof(double) * mysize_i * mysize_j);

        int rank_null;
        int rank_sender;

        for (int j_rank_index = 1; j_rank_index < dims[1]; j_rank_index++){
          MPI_Cart_shift(cart_comm, 1, j_rank_index, &rank_null, &rank_sender);

          MPI_Recv(sender_values, mysize_i * mysize_j, MPI_DOUBLE, rank_sender);

          add_values_by_line(line_data, )
        }

      }

      else{
        MPI_Send()
      }
      */
      
      char str[10];
      sprintf(str, "0%d_0%d_", coords[0], coords[1]);

      char* filename_eta = strcat(str, param.output_eta_filename);

      // TODO: togliere ghost cell prima di stampare
      write_data_vtk(&eta, "water elevation", filename_eta, n);
      //write_data_vtk(&u, "x velocity", param.output_u_filename, n);
      //write_data_vtk(&v, "y velocity", param.output_v_filename, n);
    }

    // mysize_i ultima riga 
    // mysize_j ultima colonna

    // impose boundary conditions
    double t = n * param.dt;
    if(param.source_type == 1) {
      // sinusoidal velocity on top boundary
      double A = 5;
      double f = 1. / 20.;
      //if (coords[1] == 0){
      if (1){
        for(int j = 0; j < mysize_j; j++) {
          // a sinistra
          SET(&u, 0, j, 0.);
        }
      }
      //if (coords[1] == dims[1] - 1)
      if (1){
        for(int j = 0; j < mysize_j; j++) {
          // a destra
          SET(&u, mysize_i, j, 0.);
        }
      }
      //if (coords[0] == 0)
      if (1){
        for(int i = 0; i < mysize_i; i++) {
          // sopra
          SET(&v, i, 0, 0.);
        }
      }
      //if (coords[0] == dims[0]-1)
      if (1){
        for(int i = 0; i < mysize_i; i++) {
          // sotto
          SET(&v, i, mysize_j, A * sin(2 * M_PI * f * t));
        }
      }
      printf("Rank %d: Ho messo le boundary condition", rank);
    }
    else if(param.source_type == 2) {
      // sinusoidal elevation in the middle of the domain
      double A = 5;
      double f = 1. / 20.;
      SET(&eta, nx / 2, ny / 2, A * sin(2 * M_PI * f * t));
    }
    else {
      // TODO: add other sources
      printf("Error: Unknown source type %d\n", param.source_type);
      exit(0);
    }

    // update eta
    for(int j = 0; j < mysize_j ; j++) {
      for(int i = 0; i < mysize_i; i++) {
        // TODO: this does not evaluate h at the correct locations
        double h_x_u_i1_j = GET(&h_interp_u, i+1, j);
        double h_x_u_i_j = GET(&h_interp_u, i, j);

        double h_x_v_i_j1 = GET(&h_interp_v, i, j+1);
        double h_x_v_i_j = GET(&h_interp_v, i, j);

        if (h_x_u_i1_j != 20.0 || h_x_u_i_j != 20.0 || h_x_v_i_j1 != 20.0 || h_x_v_i_j != 20){
          printf("i: %d, j: %d, value: %f | %f | %f | %f\n", i, j, h_x_u_i1_j, h_x_u_i_j, h_x_v_i_j1, h_x_v_i_j);
          return 1;
        }

        double eta_ij = GET(&eta, i, j) 
                - param.dt / param.dx * (h_x_u_i1_j * GET(&u, i+1, j) - h_x_u_i_j * GET(&u, i, j))
                - param.dt / param.dy * (h_x_v_i_j1 * GET(&v, i, j+1) - h_x_v_i_j * GET(&v, i, j));
        
        if (eta_ij > 10.0 || eta_ij < -10.0 || isnan(eta_ij)){
          printf("i: %d, j: %d, value: %f", i, j, eta_ij);
          return 1;
        }
        
        SET(&eta, i, j, eta_ij);
      }
    }

    // update u and v
    for(int j = 0; j < mysize_j; j++) {
      for(int i = 0; i < mysize_i; i++) {
        double c1 = param.dt * param.g;
        double c2 = param.dt * param.gamma;
        double eta_ij = GET(&eta, i, j);

        double eta_imj = GET(&eta, (i == 0) ? 0 : i - 1, j);
        double eta_ijm = GET(&eta, i, (j == 0) ? 0 : j - 1);
        
        double u_ij = (1. - c2) * GET(&u, i, j)
          - c1 / param.dx * (eta_ij - eta_imj);

        double v_ij = (1. - c2) * GET(&v, i, j)
          - c1 / param.dy * (eta_ij - eta_ijm);

        SET(&u, i, j, u_ij);
        SET(&v, i, j, v_ij);
      }
    }

  }

  MPI_Finalize();

  //write_manifest_vtk("water elevation", param.output_eta_filename,
  //                  param.dt, nt, param.sampling_rate);
  //write_manifest_vtk("x velocity", param.output_u_filename,
  //                   param.dt, nt, param.sampling_rate);
  //write_manifest_vtk("y velocity", param.output_v_filename,
  //                   param.dt, nt, param.sampling_rate);

  double time = GET_TIME() - start;
  printf("\nDone: %g seconds (%g MUpdates/s)\n", time,
        1e-6 * (double)eta.nx * (double)eta.ny * (double)nt / time);

  free_data(&h_interp_u);
  free_data(&h_interp_v);
  free_data(&eta);
  free_data(&u);
  free_data(&v);

  return 0;
}
