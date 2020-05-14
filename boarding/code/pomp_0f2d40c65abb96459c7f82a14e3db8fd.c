/* pomp C snippet file: pomp_0f2d40c65abb96459c7f82a14e3db8fd */
/* Time: 2020-05-14 09:38:15.696 -0400 */
/* Salt: 4E7222952DD2C1558AF5CC23 */

#include <pomp.h>
#include <R_ext/Rdynload.h>


double tol() { 
  static double val = 1e-6; 
  return val;
} 
 


/* C snippet: 'rinit' */
#define Beta		(__p[__parindex[0]])
#define Gamma		(__p[__parindex[1]])
#define Rho		(__p[__parindex[2]])
#define N		(__p[__parindex[3]])
#define S0		(__p[__parindex[4]])
#define I0		(__p[__parindex[5]])
#define R0		(__p[__parindex[6]])
#define Sigma		(__p[__parindex[7]])
#define S		(__x[__stateindex[0]])
#define I		(__x[__stateindex[1]])
#define R		(__x[__stateindex[2]])

void __pomp_rinit (double *__x, const double *__p, double t, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars)
{
 
  S = nearbyint(S0);
  I = nearbyint(I0);
  R = nearbyint(R0);
 
}

#undef Beta
#undef Gamma
#undef Rho
#undef N
#undef S0
#undef I0
#undef R0
#undef Sigma
#undef S
#undef I
#undef R

/* C snippet: 'step.fn' */
#define Beta		(__p[__parindex[0]])
#define Gamma		(__p[__parindex[1]])
#define Rho		(__p[__parindex[2]])
#define N		(__p[__parindex[3]])
#define S0		(__p[__parindex[4]])
#define I0		(__p[__parindex[5]])
#define R0		(__p[__parindex[6]])
#define Sigma		(__p[__parindex[7]])
#define S		(__x[__stateindex[0]])
#define I		(__x[__stateindex[1]])
#define R		(__x[__stateindex[2]])

void __pomp_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t, double dt)
{
 
  double dSI = rbinom(S, 1 - exp(-Beta*I/N*dt));
  double dIR = rbinom(I, 1 - exp(-Gamma*dt));
  double dRR2 = rbinom(R, 1 - exp(-Sigma*dt));

  S -= dSI;
  I += dSI - dIR;
  R += dIR - dRR2;
 
}

#undef Beta
#undef Gamma
#undef Rho
#undef N
#undef S0
#undef I0
#undef R0
#undef Sigma
#undef S
#undef I
#undef R

/* C snippet: 'rmeasure' */
#define Beta		(__p[__parindex[0]])
#define Gamma		(__p[__parindex[1]])
#define Rho		(__p[__parindex[2]])
#define N		(__p[__parindex[3]])
#define S0		(__p[__parindex[4]])
#define I0		(__p[__parindex[5]])
#define R0		(__p[__parindex[6]])
#define Sigma		(__p[__parindex[7]])
#define S		(__x[__stateindex[0]])
#define I		(__x[__stateindex[1]])
#define R		(__x[__stateindex[2]])
#define B		(__y[__obsindex[0]])
#define C		(__y[__obsindex[1]])

void __pomp_rmeasure (double *__y, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t)
{
 
  B = rpois(Rho * R + tol());
 
}

#undef Beta
#undef Gamma
#undef Rho
#undef N
#undef S0
#undef I0
#undef R0
#undef Sigma
#undef S
#undef I
#undef R
#undef B
#undef C

/* C snippet: 'dmeasure' */
#define Beta		(__p[__parindex[0]])
#define Gamma		(__p[__parindex[1]])
#define Rho		(__p[__parindex[2]])
#define N		(__p[__parindex[3]])
#define S0		(__p[__parindex[4]])
#define I0		(__p[__parindex[5]])
#define R0		(__p[__parindex[6]])
#define Sigma		(__p[__parindex[7]])
#define S		(__x[__stateindex[0]])
#define I		(__x[__stateindex[1]])
#define R		(__x[__stateindex[2]])
#define B		(__y[__obsindex[0]])
#define C		(__y[__obsindex[1]])
#define lik		(__lik[0])

void __pomp_dmeasure (double *__lik, const double *__y, const double *__x, const double *__p, int give_log, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t)
{
 
  lik = dpois(B, Rho * R + tol(), give_log);
 
}

#undef Beta
#undef Gamma
#undef Rho
#undef N
#undef S0
#undef I0
#undef R0
#undef Sigma
#undef S
#undef I
#undef R
#undef B
#undef C
#undef lik

/* C snippet: 'toEst' */
#define Beta		(__p[__parindex[0]])
#define Gamma		(__p[__parindex[1]])
#define Rho		(__p[__parindex[2]])
#define N		(__p[__parindex[3]])
#define S0		(__p[__parindex[4]])
#define I0		(__p[__parindex[5]])
#define R0		(__p[__parindex[6]])
#define Sigma		(__p[__parindex[7]])
#define T_Beta		(__pt[__parindex[0]])
#define T_Gamma		(__pt[__parindex[1]])
#define T_Rho		(__pt[__parindex[2]])
#define T_N		(__pt[__parindex[3]])
#define T_S0		(__pt[__parindex[4]])
#define T_I0		(__pt[__parindex[5]])
#define T_R0		(__pt[__parindex[6]])
#define T_Sigma		(__pt[__parindex[7]])

void __pomp_to_trans (double *__pt, const double *__p, const int *__parindex)
{
 	T_Beta = log(Beta);
	T_Gamma = log(Gamma);
	T_Sigma = log(Sigma);
	T_Rho = logit(Rho); 
}

#undef Beta
#undef Gamma
#undef Rho
#undef N
#undef S0
#undef I0
#undef R0
#undef Sigma
#undef T_Beta
#undef T_Gamma
#undef T_Rho
#undef T_N
#undef T_S0
#undef T_I0
#undef T_R0
#undef T_Sigma

/* C snippet: 'fromEst' */
#define Beta		(__p[__parindex[0]])
#define Gamma		(__p[__parindex[1]])
#define Rho		(__p[__parindex[2]])
#define N		(__p[__parindex[3]])
#define S0		(__p[__parindex[4]])
#define I0		(__p[__parindex[5]])
#define R0		(__p[__parindex[6]])
#define Sigma		(__p[__parindex[7]])
#define T_Beta		(__pt[__parindex[0]])
#define T_Gamma		(__pt[__parindex[1]])
#define T_Rho		(__pt[__parindex[2]])
#define T_N		(__pt[__parindex[3]])
#define T_S0		(__pt[__parindex[4]])
#define T_I0		(__pt[__parindex[5]])
#define T_R0		(__pt[__parindex[6]])
#define T_Sigma		(__pt[__parindex[7]])

void __pomp_from_trans (double *__p, const double *__pt, const int *__parindex)
{
 	Beta = exp(T_Beta);
	Gamma = exp(T_Gamma);
	Sigma = exp(T_Sigma);
	Rho = expit(T_Rho); 
}

#undef Beta
#undef Gamma
#undef Rho
#undef N
#undef S0
#undef I0
#undef R0
#undef Sigma
#undef T_Beta
#undef T_Gamma
#undef T_Rho
#undef T_N
#undef T_S0
#undef T_I0
#undef T_R0
#undef T_Sigma

static int __pomp_load_stack = 0;

void __pomp_load_stack_incr (void) {++__pomp_load_stack;}

void __pomp_load_stack_decr (int *val) {*val = --__pomp_load_stack;}

void R_init_pomp_0f2d40c65abb96459c7f82a14e3db8fd (DllInfo *info)
{
R_RegisterCCallable("pomp_0f2d40c65abb96459c7f82a14e3db8fd", "__pomp_load_stack_incr", (DL_FUNC) __pomp_load_stack_incr);
R_RegisterCCallable("pomp_0f2d40c65abb96459c7f82a14e3db8fd", "__pomp_load_stack_decr", (DL_FUNC) __pomp_load_stack_decr);
R_RegisterCCallable("pomp_0f2d40c65abb96459c7f82a14e3db8fd", "__pomp_rinit", (DL_FUNC) __pomp_rinit);
R_RegisterCCallable("pomp_0f2d40c65abb96459c7f82a14e3db8fd", "__pomp_stepfn", (DL_FUNC) __pomp_stepfn);
R_RegisterCCallable("pomp_0f2d40c65abb96459c7f82a14e3db8fd", "__pomp_rmeasure", (DL_FUNC) __pomp_rmeasure);
R_RegisterCCallable("pomp_0f2d40c65abb96459c7f82a14e3db8fd", "__pomp_dmeasure", (DL_FUNC) __pomp_dmeasure);
R_RegisterCCallable("pomp_0f2d40c65abb96459c7f82a14e3db8fd", "__pomp_to_trans", (DL_FUNC) __pomp_to_trans);
R_RegisterCCallable("pomp_0f2d40c65abb96459c7f82a14e3db8fd", "__pomp_from_trans", (DL_FUNC) __pomp_from_trans);
}
