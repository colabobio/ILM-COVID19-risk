/* pomp C snippet file: snippets */
/* Time: 2020-05-16 15:03:38.185 -0400 */
/* Salt: FBF9AC5CDDF8EC6EEDA175C3 */

#include <pomp.h>
#include <R_ext/Rdynload.h>


double calc_beta(double td, double a0, double a1, double b0, double b1) {
  static int *indices = NULL;
  static double *contacts = NULL;
  static int max_t = 0;
  static int num_v = 0;

  if (indices == NULL) {
    FILE *file;

    file = fopen("./gama/indices", "r");

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

    file = fopen("./gama/contacts", "r");
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

    //Rprintf("%d %d\n", max_t, num_v);
  }

  double beta = 0;

  int t = (int) td;
  if (max_t <= t) t = max_t - 1;
  int idx = indices[t];
  int ninf = 0;
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

  //Rprintf("%lg = %lg\n", td, beta);

  return beta;
}

double tol() {
  static double val = 1e-6;
  return val;
} 
 


/* C snippet: 'rinit' */
#define a0		(__p[__parindex[0]])
#define a1		(__p[__parindex[1]])
#define b0		(__p[__parindex[2]])
#define b1		(__p[__parindex[3]])
#define sigma		(__p[__parindex[4]])
#define gamma		(__p[__parindex[5]])
#define pop		(__p[__parindex[6]])
#define S_0		(__p[__parindex[7]])
#define E_0		(__p[__parindex[8]])
#define I_0		(__p[__parindex[9]])
#define R_0		(__p[__parindex[10]])
#define rho		(__p[__parindex[11]])
#define S		(__x[__stateindex[0]])
#define E		(__x[__stateindex[1]])
#define I		(__x[__stateindex[2]])
#define R		(__x[__stateindex[3]])
#define C		(__x[__stateindex[4]])

void __pomp_rinit (double *__x, const double *__p, double t, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars)
{
 
  double m = pop/(S_0 + E_0 + I_0 + R_0);

  S = nearbyint(m*S_0);
  E = nearbyint(m*E_0);
  I = nearbyint(m*I_0);
  R = nearbyint(m*R_0);

  C = 0;
 
}

#undef a0
#undef a1
#undef b0
#undef b1
#undef sigma
#undef gamma
#undef pop
#undef S_0
#undef E_0
#undef I_0
#undef R_0
#undef rho
#undef S
#undef E
#undef I
#undef R
#undef C

/* C snippet: 'step.fn' */
#define a0		(__p[__parindex[0]])
#define a1		(__p[__parindex[1]])
#define b0		(__p[__parindex[2]])
#define b1		(__p[__parindex[3]])
#define sigma		(__p[__parindex[4]])
#define gamma		(__p[__parindex[5]])
#define pop		(__p[__parindex[6]])
#define S_0		(__p[__parindex[7]])
#define E_0		(__p[__parindex[8]])
#define I_0		(__p[__parindex[9]])
#define R_0		(__p[__parindex[10]])
#define rho		(__p[__parindex[11]])
#define S		(__x[__stateindex[0]])
#define E		(__x[__stateindex[1]])
#define I		(__x[__stateindex[2]])
#define R		(__x[__stateindex[3]])
#define C		(__x[__stateindex[4]])

void __pomp_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t, double dt)
{
 
  double beta; 
  double foi;
  double rate[3], trans[3];

  beta = calc_beta(t, a0, a1, b0, b1);

  // expected force of infection
  foi = beta * I/pop;

  rate[0] = foi;      // stochastic force of infection
  rate[1] = sigma;    // rate of ending of latent stage
  rate[2] = gamma;    // recovery

  // transitions between classes
  reulermultinom(1, S, &rate[0], dt, &trans[0]);
  reulermultinom(1, E, &rate[1], dt, &trans[1]);
  reulermultinom(1, I, &rate[2], dt, &trans[2]);

  S += -trans[0];
  E += trans[0] - trans[1];
  I += trans[1] - trans[2];
  R = pop - S - E - I;

  // Assigning the right number to the accumulation variable that's used
  // in the observation model is absolutely critical!!!!
  C += trans[2];
 
}

#undef a0
#undef a1
#undef b0
#undef b1
#undef sigma
#undef gamma
#undef pop
#undef S_0
#undef E_0
#undef I_0
#undef R_0
#undef rho
#undef S
#undef E
#undef I
#undef R
#undef C

/* C snippet: 'rmeasure' */
#define a0		(__p[__parindex[0]])
#define a1		(__p[__parindex[1]])
#define b0		(__p[__parindex[2]])
#define b1		(__p[__parindex[3]])
#define sigma		(__p[__parindex[4]])
#define gamma		(__p[__parindex[5]])
#define pop		(__p[__parindex[6]])
#define S_0		(__p[__parindex[7]])
#define E_0		(__p[__parindex[8]])
#define I_0		(__p[__parindex[9]])
#define R_0		(__p[__parindex[10]])
#define rho		(__p[__parindex[11]])
#define S		(__x[__stateindex[0]])
#define E		(__x[__stateindex[1]])
#define I		(__x[__stateindex[2]])
#define R		(__x[__stateindex[3]])
#define C		(__x[__stateindex[4]])
#define cases		(__y[__obsindex[0]])

void __pomp_rmeasure (double *__y, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t)
{
 
  cases = rbinom(I, rho);
 
}

#undef a0
#undef a1
#undef b0
#undef b1
#undef sigma
#undef gamma
#undef pop
#undef S_0
#undef E_0
#undef I_0
#undef R_0
#undef rho
#undef S
#undef E
#undef I
#undef R
#undef C
#undef cases

/* C snippet: 'dmeasure' */
#define a0		(__p[__parindex[0]])
#define a1		(__p[__parindex[1]])
#define b0		(__p[__parindex[2]])
#define b1		(__p[__parindex[3]])
#define sigma		(__p[__parindex[4]])
#define gamma		(__p[__parindex[5]])
#define pop		(__p[__parindex[6]])
#define S_0		(__p[__parindex[7]])
#define E_0		(__p[__parindex[8]])
#define I_0		(__p[__parindex[9]])
#define R_0		(__p[__parindex[10]])
#define rho		(__p[__parindex[11]])
#define S		(__x[__stateindex[0]])
#define E		(__x[__stateindex[1]])
#define I		(__x[__stateindex[2]])
#define R		(__x[__stateindex[3]])
#define C		(__x[__stateindex[4]])
#define cases		(__y[__obsindex[0]])
#define lik		(__lik[0])

void __pomp_dmeasure (double *__lik, const double *__y, const double *__x, const double *__p, int give_log, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t)
{
 
  lik = dbinom(cases, I, rho, give_log);
 
}

#undef a0
#undef a1
#undef b0
#undef b1
#undef sigma
#undef gamma
#undef pop
#undef S_0
#undef E_0
#undef I_0
#undef R_0
#undef rho
#undef S
#undef E
#undef I
#undef R
#undef C
#undef cases
#undef lik

/* C snippet: 'toEst' */
#define a0		(__p[__parindex[0]])
#define a1		(__p[__parindex[1]])
#define b0		(__p[__parindex[2]])
#define b1		(__p[__parindex[3]])
#define sigma		(__p[__parindex[4]])
#define gamma		(__p[__parindex[5]])
#define pop		(__p[__parindex[6]])
#define S_0		(__p[__parindex[7]])
#define E_0		(__p[__parindex[8]])
#define I_0		(__p[__parindex[9]])
#define R_0		(__p[__parindex[10]])
#define rho		(__p[__parindex[11]])
#define T_a0		(__pt[__parindex[0]])
#define T_a1		(__pt[__parindex[1]])
#define T_b0		(__pt[__parindex[2]])
#define T_b1		(__pt[__parindex[3]])
#define T_sigma		(__pt[__parindex[4]])
#define T_gamma		(__pt[__parindex[5]])
#define T_pop		(__pt[__parindex[6]])
#define T_S_0		(__pt[__parindex[7]])
#define T_E_0		(__pt[__parindex[8]])
#define T_I_0		(__pt[__parindex[9]])
#define T_R_0		(__pt[__parindex[10]])
#define T_rho		(__pt[__parindex[11]])

void __pomp_to_trans (double *__pt, const double *__p, const int *__parindex)
{
 	T_a0 = log(a0);
	T_a1 = log(a1);
	T_b0 = log(b0);
	T_b1 = log(b1);
	T_sigma = log(sigma);
	T_gamma = log(gamma); 
}

#undef a0
#undef a1
#undef b0
#undef b1
#undef sigma
#undef gamma
#undef pop
#undef S_0
#undef E_0
#undef I_0
#undef R_0
#undef rho
#undef T_a0
#undef T_a1
#undef T_b0
#undef T_b1
#undef T_sigma
#undef T_gamma
#undef T_pop
#undef T_S_0
#undef T_E_0
#undef T_I_0
#undef T_R_0
#undef T_rho

/* C snippet: 'fromEst' */
#define a0		(__p[__parindex[0]])
#define a1		(__p[__parindex[1]])
#define b0		(__p[__parindex[2]])
#define b1		(__p[__parindex[3]])
#define sigma		(__p[__parindex[4]])
#define gamma		(__p[__parindex[5]])
#define pop		(__p[__parindex[6]])
#define S_0		(__p[__parindex[7]])
#define E_0		(__p[__parindex[8]])
#define I_0		(__p[__parindex[9]])
#define R_0		(__p[__parindex[10]])
#define rho		(__p[__parindex[11]])
#define T_a0		(__pt[__parindex[0]])
#define T_a1		(__pt[__parindex[1]])
#define T_b0		(__pt[__parindex[2]])
#define T_b1		(__pt[__parindex[3]])
#define T_sigma		(__pt[__parindex[4]])
#define T_gamma		(__pt[__parindex[5]])
#define T_pop		(__pt[__parindex[6]])
#define T_S_0		(__pt[__parindex[7]])
#define T_E_0		(__pt[__parindex[8]])
#define T_I_0		(__pt[__parindex[9]])
#define T_R_0		(__pt[__parindex[10]])
#define T_rho		(__pt[__parindex[11]])

void __pomp_from_trans (double *__p, const double *__pt, const int *__parindex)
{
 	a0 = exp(T_a0);
	a1 = exp(T_a1);
	b0 = exp(T_b0);
	b1 = exp(T_b1);
	sigma = exp(T_sigma);
	gamma = exp(T_gamma); 
}

#undef a0
#undef a1
#undef b0
#undef b1
#undef sigma
#undef gamma
#undef pop
#undef S_0
#undef E_0
#undef I_0
#undef R_0
#undef rho
#undef T_a0
#undef T_a1
#undef T_b0
#undef T_b1
#undef T_sigma
#undef T_gamma
#undef T_pop
#undef T_S_0
#undef T_E_0
#undef T_I_0
#undef T_R_0
#undef T_rho

static int __pomp_load_stack = 0;

void __pomp_load_stack_incr (void) {++__pomp_load_stack;}

void __pomp_load_stack_decr (int *val) {*val = --__pomp_load_stack;}

void R_init_snippets (DllInfo *info)
{
R_RegisterCCallable("snippets", "__pomp_load_stack_incr", (DL_FUNC) __pomp_load_stack_incr);
R_RegisterCCallable("snippets", "__pomp_load_stack_decr", (DL_FUNC) __pomp_load_stack_decr);
R_RegisterCCallable("snippets", "__pomp_rinit", (DL_FUNC) __pomp_rinit);
R_RegisterCCallable("snippets", "__pomp_stepfn", (DL_FUNC) __pomp_stepfn);
R_RegisterCCallable("snippets", "__pomp_rmeasure", (DL_FUNC) __pomp_rmeasure);
R_RegisterCCallable("snippets", "__pomp_dmeasure", (DL_FUNC) __pomp_dmeasure);
R_RegisterCCallable("snippets", "__pomp_to_trans", (DL_FUNC) __pomp_to_trans);
R_RegisterCCallable("snippets", "__pomp_from_trans", (DL_FUNC) __pomp_from_trans);
}
