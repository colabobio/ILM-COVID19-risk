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
    double val;
    while (fscanf(file, "%lf", &val) > 0) num_v++;
    rewind(file);

    contacts = (double *)malloc(sizeof(double)*num_v);
    i = 0;
    while (fscanf(file, "%lf", &val) > 0) {
      contacts[i] = val;
      i++;
    }
    fclose(file);

    printf("Done.\n");
  }

  return max_t;
}

void print_data() {
    // for t in range(0, len(indices)):
    //    idx = indices[t]
    //    if idx < 0: continue
    //    print("=======> time", t)
    //    while -1 < values[idx]:
    //        ncont = int(values[idx])
    //        idx += 1
    //        inf = values[idx]
    //        idx += 1       
    //        print(inf, end =" ")
    //        susc = values[idx:idx+ncont]
    //        idx += ncont
    //        print(susc)

  for (int t = 0; t < max_t; t++) {
    printf("=======> time %d\n", t);
    int idx = indices[t];
    while (-1 < contacts[idx]) {
      int ncont = (int) contacts[++idx];
      double y = contacts[++idx];
      printf("%f [", y);
      for (int i = 0; i < ncont; i++) {
        double x = contacts[++idx];
        printf("%f,", x);
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