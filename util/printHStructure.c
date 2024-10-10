#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>

struct data
{
    int nx, ny;
    double dx, dy;
    double *values;
};

int read_data(struct data *data, const char *filename)
{
    FILE *fp = fopen(filename, "rb");
    if (!fp)
    {
        printf("Error: Could not open input data file '%s'\n", filename);
        return 1;
    }
    int ok = 1;
    if (ok)
        ok = (fread(&data->nx, sizeof(int), 1, fp) == 1);
    if (ok)
        ok = (fread(&data->ny, sizeof(int), 1, fp) == 1);
    if (ok)
        ok = (fread(&data->dx, sizeof(double), 1, fp) == 1);
    if (ok)
        ok = (fread(&data->dy, sizeof(double), 1, fp) == 1);
    if (ok)
    {
        int N = data->nx * data->ny;
        if (N <= 0)
        {
            printf("Error: Invalid number of data points %d\n", N);
            ok = 0;
        }
        else
        {
            data->values = (double *)malloc(N * sizeof(double));
            if (!data->values)
            {
                printf("Error: Could not allocate data (%d doubles)\n", N);
                ok = 0;
            }
            else
            {
                ok = (fread(data->values, sizeof(double), N, fp) == N);
            }
        }
    }
    fclose(fp);
    if (!ok)
    {
        printf("Error reading input data file '%s'\n", filename);
        return 1;
    }
    return 0;
}

void free_data(struct data *data)
{
    free(data->values);
}

void print_data(const struct data *data)
{
    printf("nx: %d\n", data->nx);
    printf("ny: %d\n", data->ny);
    printf("dx: %f\n", data->dx);
    printf("dy: %f\n", data->dy);
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

struct data print_data_to_console(){
    struct data h;
    char input_h_filename[20] = "h_simple.dat";

    read_data(&h, input_h_filename);

    print_data(&h);

    free_data(&h);

    return h;
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

void save_coordinate(int nx, int ny, double dx, double dy, const char *filename){
    FILE* f = fopen(filename, "w");

    fprintf(f, "x;y\n");
    for(int j = 0; j < ny ; j++) {
        for(int i = 0; i < nx; i++) {
            double x = ((double)i + 0.5) * dx;
            double y = j * dy;

            fprintf(f, "%f;%f\n", x, y);
        }
    }

    fclose(f);
}

int main(int argc, char **argv)
{
    struct data h = print_data_to_console();

    double hx = h.nx * h.dx;
    double hy = h.ny * h.dy;
    int nx = floor(hx / 5.0);
    int ny = floor(hy / 5.0);

    printf("hx: %f, hy: %f, nx: %d, ny: %d\n", hx, hy, nx, ny);

    save_coordinate(nx, ny, 5.0, 5.0, "v_coordinates.csv");
}