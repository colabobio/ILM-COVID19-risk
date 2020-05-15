#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static int *indices = NULL;
static double *contacts = NULL;
static int max_t = 0;
static int num_v = 0;

int load_data() {
  if (indices == NULL) {
    printf("Reading data...\n");

    FILE *file;

    file = fopen("indices", "r");

    int idx;
    while (fscanf(file, "%d", &idx) > 0) max_t++;
    rewind(file);

    indices = (int *)malloc(sizeof(int)*max_t);
    int i = 0;
    while (fscanf(file, "%d", &idx) > 0) {
      indices[i] = idx;
      i++;
    }
    fclose(file);

    file = fopen("contacts", "r");
    float val;
    while (fscanf(file, "%f", &val) > 0) num_v++;
    rewind(file);

    contacts = (double *)malloc(sizeof(double)*num_v);
    i = 0;
    while (fscanf(file, "%f", &val) > 0) {
      contacts[i] = val;
      i++;
    }
    fclose(file);

    printf("Done.\n");
  }

  return max_t;
}

void print_data() {
  for (int t = 0; t < max_t; t++) {
    printf("=======> time %d\n", t);
    int idx = indices[t];
    while (-1 < contacts[idx]) {
      int ncont = (int) contacts[idx++];
      double y = contacts[idx++];
      printf("%.1f [", y);
      for (int i = 0; i < ncont; i++) {
        double x = contacts[idx++];
        printf("%.1f", x);
        if (i < ncont - 1) printf(", ");
      }
      printf("]\n");
    }
  }
}

void print_beta() {
  double a0 = 0.1; // should be low 
  double a1 = 10.0; // should be high according to the importance of covariate.
  double b0 = 0.1; // only considering 2 covariates as of now for both infectivity and succeptibility
  double b1 = 10.0;

  for (int t = 0; t < max_t; t++) {

    printf("=======> time %d\n", t);
    int idx = indices[t];
    int ninf = 0;
    double beta = 0;
    while (-1 < contacts[idx]) {
      int ncont = (int) contacts[idx++];
      double y = contacts[idx++];
      for (int i = 0; i < ncont; i++) {
        double x = contacts[idx++];
        double p = (a0 + a1 * x) * (b0 + b1 * y);
        beta += 1 - exp(-p);
      }
      ninf++;
    }
    if (0 < ninf) {
      beta /= ninf;
    }
    printf("%.3f\n", beta);
  }
}

int main() {
   int res = load_data();
   printf("Got %d\n", res);
   print_data();
   print_beta();
   return 0;
}