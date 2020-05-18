#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static int *indices = NULL;
static double *contacts = NULL;
static int max_t = 0;
static int num_v = 0;

int load_data(const char *dir) {
  if (indices == NULL) {
    printf("Reading data...\n");

    char *indices_fn = (char *) malloc(1 + strlen(dir)+ strlen("/indices"));
    strcpy(indices_fn, dir);
    strcat(indices_fn, "/indices");

    char *contacts_fn = (char *) malloc(1 + strlen(dir)+ strlen("/contacts"));
    strcpy(contacts_fn, dir);
    strcat(contacts_fn, "/contacts");

    FILE *file;

    file = fopen(indices_fn, "r");

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

    file = fopen(contacts_fn, "r");
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

void print_beta(double a0, double a1, double b0, double b1) {
  printf("\nBETA ESTIMATES\n");

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

int main(int argc, char *argv[]) {
  if (argc < 2) { // argc should be at least 2 for correct execution
    printf("usage: %s directory\n", argv[0]);
  } else {
    int res = load_data(argv[1]);
    printf("Got %d\n", res);
    print_data();
    if (argc == 6) {
      double a0 = strtod(argv[2], NULL);
      double a1 = strtod(argv[3], NULL);
      double b0 = strtod(argv[4], NULL);
      double b1 = strtod(argv[5], NULL);
      print_beta(a0, a1, b0, b1);
    }    
    return 0;
  }
}