#include <stdio.h>
#include <stdlib.h>

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

int main() {
   int res = load_data();
   printf("Got %d\n", res);
   print_data();
   return 0;
}