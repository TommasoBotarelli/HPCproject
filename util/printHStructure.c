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

int main(int argc, char **argv)
{
    struct data h;
    char input_h_filename[20] = "h_simple.dat";

    if (read_data(&h, input_h_filename))
        return 1;

    print_data(&h);

    free_data(&h);
}