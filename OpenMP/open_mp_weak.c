#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

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

#define GET(data, i, j) ((data)->values[(data)->nx * (j) + (i)])
#define SET(data, i, j, val) ((data)->values[(data)->nx * (j) + (i)] = (val))

#define GET_X_COORD(data, i) ((data)->dx * (i))
#define GET_Y_COORD(data, i) ((data)->dy * (i))

#define ERROR_VALUE -1

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

int update_parameters_for_weak_experiment(struct parameters* param,
                                           const char *ranks_number)
{
  int ranks = strtol(ranks_number, NULL, 10);

  if (ranks == 2){
    param->dx = param->dx / 2;
  }
  if (ranks == 4){
    param->dx = param->dx / 2;
    param->dy = param->dy / 2;
  }
  if (ranks == 8){
    param->dx = param->dx / 4;
    param->dy = param->dy / 2;
  }
  if (ranks == 16){
    param->dx = param->dx / 4;
    param->dy = param->dy / 4;
  }
  if (ranks == 32){
    param->dx = param->dx / 8;
    param->dy = param->dy / 4;
  }
  if (ranks == 64){
    param->dx = param->dx / 8;
    param->dy = param->dy / 8;
  }
  if (ranks > 64){
    return 1;
  }
  return 0;
}

int main(int argc, char **argv)
{
  
  if(argc != 3) {
    printf("Usage: %s parameter_file, number of ranks\n", argv[0]);
    return 1;
  }

  struct parameters param;
  if(read_parameters(&param, argv[1])) return 1;

  int result = update_parameters_for_weak_experiment(&param, argv[2]);
  
  if (result == 1){
    printf("Number of ranks greater (%s) than expected", argv[2]);
    return 1;
  }

  print_parameters(&param);

  struct data h;
  if(read_data(&h, param.input_h_filename)) return 1;

  printf("h->dx = %f\n", h.dx);
  printf("h->dy = %f\n", h.dy);

  // infer size of domain from input elevation data
  double hx = h.nx * h.dx;
  double hy = h.ny * h.dy;
  int nx = floor(hx / param.dx);
  int ny = floor(hy / param.dy);
  if(nx <= 0) nx = 1;
  if(ny <= 0) ny = 1;
  int nt = floor(param.max_t / param.dt);

  printf(" - grid size: %g m x %g m (%d x %d = %d grid points)\n",
         hx, hy, nx, ny, nx * ny);
  printf(" - number of time steps: %d\n", nt);

  struct data eta, u, v;
  init_data(&eta, nx, ny, param.dx, param.dx, 0.);
  init_data(&u, nx + 1, ny, param.dx, param.dy, 0.);
  init_data(&v, nx, ny + 1, param.dx, param.dy, 0.);

  // interpolate bathymetry
  struct data h_interp_u;
  struct data h_interp_v;

  init_data(&h_interp_u, u.nx, u.ny, param.dx, param.dy, 0.);

  #pragma omp parallel for collapse(2)
  for(int j = 0; j < u.ny ; ++j) {
    for(int i = 0; i < u.nx; ++i) {
      double x = i * param.dx;
      double y = ((double)j + 0.5) * param.dy;
      double val = interpolate_data(&h, x, y);
      
      SET(&h_interp_u, i, j, val);
    }
  }

  init_data(&h_interp_v, v.nx, v.ny, param.dx, param.dy, 0.);

  #pragma omp parallel for collapse(2)
  for(int j = 0; j < v.ny; ++j) {
    for(int i = 0; i < v.nx ; ++i) {
      double x = ((double)i + 0.5) * param.dx;
      double y = j * param.dy;
      double val = interpolate_data(&h, x, y);
      
      SET(&h_interp_v, i, j, val);
    }
  }

  double start = GET_TIME();

  for(int n = 0; n < nt; n++) {

    if(n && (n % (nt / 10)) == 0) {
      double time_sofar = GET_TIME() - start;
      double eta = (nt - n) * time_sofar / n;
      printf("Computing step %d/%d (ETA: %g seconds)     \r", n, nt, eta);
      fflush(stdout);
    }

    // output solution
    if(param.sampling_rate && !(n % param.sampling_rate)) {
      //write_data_vtk(&eta, "water elevation", param.output_eta_filename, n);
      //write_data_vtk(&u, "x velocity", param.output_u_filename, n);
      //write_data_vtk(&v, "y velocity", param.output_v_filename, n);
    }

    
    // impose boundary conditions

    double t = n * param.dt;
    if(param.source_type == 1) {
      double A = 5;
      double f = 1. / 20.;
      // sinusoidal velocity on top boundary

      for(int j = 0; j < ny; j++) {
        SET(&u, 0, j, 0.);
        SET(&u, nx, j, 0.);
      }
      for(int i = 0; i < nx; i++) {
        SET(&v, i, 0, 0.);
        SET(&v, i, ny, A * sin(2 * M_PI * f * t));
      }
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
    #pragma omp parallel for collapse(2)
    for(int j = 0; j < ny ; ++j) {
      for(int i = 0; i < nx; ++i) {

        double h_x_u_i1_j = GET(&h_interp_u, i+1, j);
        double h_x_u_i_j = GET(&h_interp_u, i, j);

        double h_x_v_i_j1 = GET(&h_interp_v, i, j+1);
        double h_x_v_i_j = GET(&h_interp_v, i, j);

        //Fare la divisione fuori dai for

        double eta_ij = GET(&eta, i, j) 
                - param.dt / param.dx * (h_x_u_i1_j * GET(&u, i+1, j) - h_x_u_i_j * GET(&u, i, j))
                - param.dt / param.dy * (h_x_v_i_j1 * GET(&v, i, j+1) - h_x_v_i_j * GET(&v, i, j));
        
        
        SET(&eta, i, j, eta_ij);
      }
    }

    // Spostare fuori da for n

    double c1_dx = param.dt * param.g / param.dx;
    double c1_dy = param.dt * param.g / param.dy;
    double one_minus_c2 = 1.0 - param.dt * param.gamma;

    // update u and v
    #pragma omp parallel for collapse(2)
    for(int j = 0; j < ny; ++j) {
      for(int i = 0; i < nx; ++i) {
        double eta_ij = GET(&eta, i, j);
        double eta_imj = GET(&eta, (i == 0) ? 0 : i - 1, j);
        double eta_ijm = GET(&eta, i, (j == 0) ? 0 : j - 1);
        double u_ij = one_minus_c2 * GET(&u, i, j)
          -c1_dx * (eta_ij - eta_imj);
        double v_ij = one_minus_c2 * GET(&v, i, j)
          - c1_dy * (eta_ij - eta_ijm);
        SET(&u, i, j, u_ij);
        SET(&v, i, j, v_ij);
      }
    }

  }

  //write_manifest_vtk("water elevation", param.output_eta_filename,
  //                   param.dt, nt, param.sampling_rate);
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
