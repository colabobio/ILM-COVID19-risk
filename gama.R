args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  main_folder <- "."
} else {
  main_folder <- args[1]  
}

# =============================================================================
# Libraries

library(doRNG)
library(foreach)
library(doParallel)

# tidyr has to be imported before magrittr so extract() from the latter is used
library(tidyr)

library(dplyr)
library(plyr)
library(reshape2)
library(magrittr)

library(ggplot2)
theme_set(theme_bw())

library(pomp)
stopifnot(packageVersion("pomp")>="2")

# =============================================================================
# Folders

output_folder <- file.path(main_folder, "output") 
if (!dir.exists(output_folder)) dir.create(output_folder)
cooking_folder <- file.path(main_folder, "output", "bake")
if (!dir.exists(cooking_folder)) dir.create(cooking_folder)
plotting_folder <- file.path(main_folder, "output", "plots")
if (!dir.exists(plotting_folder)) dir.create(plotting_folder)
code_folder <- file.path(main_folder, "code")
if (!dir.exists(code_folder)) dir.create(code_folder)
file_name <- "snippets"

# =============================================================================
# Observed data
csv_table <- read.csv(file.path(main_folder, "case_counts.csv"))
case_data <- data.frame(day = csv_table$Time, cases = csv_table$Count)
head(case_data)

max_day <- max(case_data$day)
ggplot(data=subset(case_data, day <= max_day), aes(x=day, y=cases, group=1)) + geom_line()

ggsave(file.path(plotting_folder, "0-observed_data.pdf"))

# =============================================================================
# Parameters

pop_size <- 9205
exp0 <- 10

time_step <- 1/4

num_sims <- 100

free_param_names <- c("a0", "a1", "b0", "b1", "sigma", "gamma")
free_param_box <- rbind(
  a0 = c(0, 1),
  a1 = c(0, 20),
  b0 = c(0, 1),
  b1 = c(0, 20),
  sigma = c(0, 2),
  gamma = c(0, 2)
)

#free_param_names <- c("beta", "gamma")
#free_param_box <- rbind(
#  beta = c(0, 4),
#  gamma = c(0, 2)
#)

fixed_param_names <- c("pop", "S_0", "E_0", "I_0", "R_0", "rho")
fixed_param_values <- c(pop=pop_size, S_0=1-exp0/pop_size, E_0=exp0/pop_size, I_0=0, R_0=0, rho=1.0)

all_param_names <- c(free_param_names, fixed_param_names)

# Random seeds, keep unchanged to ensure reproducibilty of results
test_mle_seed <- 0998468235L
full_mle_seed <- 0998468235L
global_search_seed <- 290860873
test_sim_seed <- 157999
full_sim_seed <- 157999
mcap_plik_seed <- 290860873

# =============================================================================
# Csnippets defining the SEIR model

# pomp C API:
# https://kingaa.github.io/pomp/vignettes/C_API.html

sir_step <- Csnippet("
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

  C += trans[2];
")

sir_init <- Csnippet("
  double m = pop/(S_0 + E_0 + I_0 + R_0);

  S = nearbyint(m*S_0);
  E = nearbyint(m*E_0);
  I = nearbyint(m*I_0);
  R = nearbyint(m*R_0);

  C = 0;
")

rmeas <- Csnippet("
  cases = rbinom(I, rho);
")

dmeas <- Csnippet("
  lik = dbinom(cases, I, rho, give_log);
")

extra <- Csnippet("
double calc_beta(double td, double a0, double a1, double b0, double b1) {
  static int *indices = NULL;
  static double *contacts = NULL;
  static int max_t = 0;
  static int num_v = 0;

  if (indices == NULL) {
    FILE *file;

    file = fopen(\"./gama/indices\", \"r\");

    int idx;
    while (fscanf(file, \"%d\", &idx) > 0) max_t++;
    rewind(file);

    indices = (int *)malloc(sizeof(int)*max_t);
    int i = 0;
    while (fscanf(file, \"%d\", &idx) > 0) {
      indices[i] = idx;
      i++;
    }
    fclose(file);

    file = fopen(\"./gama/contacts\", \"r\");
    float val;
    while (fscanf(file, \"%f\", &val) > 0) num_v++;
    rewind(file);

    contacts = (double *)malloc(sizeof(double)*num_v);
    i = 0;
    while (fscanf(file, \"%f\", &val) > 0) {
      contacts[i] = val;
      i++;
    }
    fclose(file);

    //Rprintf(\"%d %d\\n\", max_t, num_v);
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

  //Rprintf(\"%lg = %lg\\n\", td, beta);

  return beta;
}

double tol() {
  static double val = 1e-6;
  return val;
} 
")

# =============================================================================
# The documentation of the pomp object
# https://kingaa.github.io/pomp/manual/pomp.html

case_data %>% 
  pomp(t0 = case_data$day[1],
       time = "day",
       rprocess = euler(sir_step, delta.t=time_step),
       rinit = sir_init,
       rmeasure = rmeas,
       dmeasure = dmeas,
       globals = extra,
       cdir = code_folder,
       cfile = file_name,
       accumvars=c("C"),
       statenames=c("S", "E", "I", "R", "C"),
       partrans=parameter_trans(
         log=c("a0", "a1", "b0", "b1", "sigma", "gamma")),
         #log=c("beta", "gamma")),
       paramnames = c(free_param_names, fixed_param_names),
       #compile=FALSE,
       verbose = TRUE
  ) -> model

plot(model, main="")

# =============================================================================
# IF parameters, see more details in the manual and tutorial:
# https://kingaa.github.io/pomp/manual/mif2.html
# https://kingaa.github.io/sbied/mif/mif.html#choosing-the-algorithmic-settings-for-if2

num_test_runs <- 10

num_guesses <- 100       # Number of starting points for the parameter guesses
num_filter_iter <- 100    # Number of filtering iterations to perform
num_particles <- 5000    # Number of particles to use in filtering.
num_replicates <- 50     # Number of replicated particle filters at each point estimate

perturb_sizes <- list(a0=0.02, a1=0.02, b0=0.02, b1=0.02, sigma=0.02, gamma=0.02)
#perturb_sizes <- list(beta=0.02, gamma=0.02)
cool_frac <- 0.5
cool_type <- "geometric"

# Variables to use in the scatterplot matrix showing the result of the IF search
pair_vars <- ~loglik+a0+a1+b0+b1+sigma+gamma
#pair_vars <- ~loglik+beta+gamma

# =============================================================================
# Test run from single starting point in parameter space and no replicates

registerDoParallel()
#set.seed(test_mle_seed, kind="L'Ecuyer")

guess <- apply(free_param_box, 1, function(x) runif(1, x[1], x[2]))

bake(file=file.path(cooking_folder, "box_search_local.rds"), {
  foreach(i=1:num_test_runs,
          .packages='pomp', .combine=c, .options.multicore=list(set.seed=TRUE)
  ) %dopar%  
  {
    mif2(
      model,
      params=c(guess, fixed_param_values),
      Np=num_particles,
      Nmif=num_filter_iter,
      cooling.type=cool_type,
      cooling.fraction.50=cool_frac,
      rw.sd=do.call(rw.sd, perturb_sizes)
    )
  }
}) -> mifs_test

ggplot(data=melt(traces(mifs_test)),
       aes(x=iteration, y=value, group=L1, color=factor(L1))) +
  geom_line() +
  guides(color=FALSE) +
  facet_wrap(~variable, scales="free_y") +
  theme_bw()

ggsave(file.path(plotting_folder, "1-mle_local_search.pdf"))


# =============================================================================
# Full MLE with multiple starting points for the free parameters

registerDoParallel()
set.seed(full_mle_seed, kind="L'Ecuyer")

stew(file=file.path(cooking_folder, "box_search_global.rda"), {
  param_guesses <- as.data.frame(apply(free_param_box, 1, function(x) runif(num_guesses, x[1], x[2])))
  
  workers <- getDoParWorkers()
  systime <- system.time({
    res_global <- foreach(guess=iter(param_guesses,"row"), 
                         .packages='pomp', .combine=rbind, .options.multicore=list(set.seed=TRUE)
    ) %dopar% {
      mifs_local <- mif2(
        model,
        params=c(unlist(guess), fixed_param_values),
        Np=num_particles,
        Nmif=num_filter_iter,
        cooling.type=cool_type,
        cooling.fraction.50=cool_frac,
        rw.sd=do.call(rw.sd, perturb_sizes)
        )
      ll <- logmeanexp(replicate(num_replicates, logLik(pfilter(mifs_local, Np=num_particles))), se = TRUE)
      data.frame(as.list(coef(mifs_local)), loglik=ll[1], loglik.se=ll[2])
    }
  })
}, seed=global_search_seed, kind="L'Ecuyer")

res_global <- as.data.frame(res_global)

all <- ldply(list(guess=param_guesses, result=subset(res_global, loglik > max(loglik)-50)), .id="type")

pairs(pair_vars, data=all, col=ifelse(all$type=="guess", grey(0.5), "red"), pch=16)
dev.copy(pdf, file.path(plotting_folder, "2-mle_global_search.pdf"))
dev.off()

# =============================================================================
# Getting best parameters

mle_params <- arrange(rbind(res_global), -loglik)
write.csv(mle_params, file=file.path(output_folder, "param_mle_global_search.csv"), row.names=FALSE, na="")

log_idx <- length(mle_params) - 1
mle_global <- mle_params[which.max( mle_params[,log_idx] ), ]
mle_global %>% extract(all_param_names) %>% unlist() -> theta
write.csv(theta, file=file.path(output_folder, "param_point_estimates.csv"), row.names=TRUE, na="")
theta

# =============================================================================
# Running simulations using the MLE parameters

#set.seed(test_sim_seed)

model  %>%
  simulate(params=theta, nsim=9, format="data.frame", include.data=TRUE) %>%
  ggplot(aes(x=day, y=cases, group=.id, color=(.id=="data"))) +
  guides(color=FALSE) +
  geom_line() + facet_wrap(~.id, ncol=2)


ggsave(file.path(plotting_folder, "3-simulations.pdf"))  

# =============================================================================
# Computing a large number of simulations

set.seed(full_sim_seed)

# Getting the median cumulative curve for the simulations
model %>% 
  simulate(params=theta, nsim=num_sims, format = "data.frame", include.data=TRUE) -> sim_data

# =============================================================================
# Compare the observed data with the simulated data between the 5th and 95th percentiles

sim_data %>%
  select(day, .id, cases) %>%
  mutate(data=.id=="data") %>%
  ddply(~day+data, plyr::summarize,
    p=c(0.05, 0.5, 0.95), q=quantile(cases, prob=p, names=FALSE)) %>%
  mutate(
    p=mapvalues(p, from=c(0.05, 0.5, 0.95), to=c("lo", "med", "hi")),
    data=mapvalues(data, from=c(TRUE, FALSE), to=c("data", "simulation"))
  ) %>%
  spread(p, q) %>%
  ggplot(aes(x=day, y=med, color=data, fill=data, ymin=lo, ymax=hi)) +
         geom_ribbon(alpha=0.2)

ggsave(file.path(plotting_folder, "4-simulated_percentiles.pdf"))

# =============================================================================
# Some utility functions to calculate cumulative case numbers

cumulative_curve <- function(dat, len) {
  total_sum <- 0
  daily_sum <- c()
  for (i in 1:len) {
    total_sum <- total_sum + dat$cases[i]
    daily_sum <- c(daily_sum, total_sum)
  }
  return(daily_sum)
}  
  
median_simulation <- function(sdat, n) {
  all_totals <- c()
  for (i in 1:n) {
    sim <- subset(sdat, .id == i)
    tot <- sum(sim$cases)
    all_totals <- c(all_totals, tot)
  }

  # Taking the median
  n2 <- 0.5 * n
  median_idx <- order(all_totals)[n2]
  median_sim <- subset(sdat, .id == median_idx)
  
  return(median_sim)
}

# =============================================================================
# Comparing cumulative actual data with real data

med_sim <- median_simulation(sim_data, num_sims)
csim <- cumulative_curve(med_sim, max_day)

# Getting the cumulative curve for observed data
cobs <- cumulative_curve(case_data, max_day)

df <- data.frame('day' = seq(1, max_day), 
                  'obs_data' = cobs,
                  'sim_data' = csim)

ggplot(df, aes(day)) + 
  geom_line(aes(y = obs_data, colour = "Real Data")) + 
  geom_line(aes(y = sim_data, colour = "Simulated")) +
  ylab('Cumulative Number of Cases')

ggsave(file.path(plotting_folder, "4-cumulative_cases.pdf"))

# =============================================================================
# CALCULATION OF THE CONFIDENCE INTERVAL FOR THE PARAMETERS (OPTIONAL)
# =============================================================================

# =============================================================================
# MCAP settings

mcap_confidence <- 0.95  # The desired confidence

# Profile design settings
mcap_profdes_len <- 15   # Number of subdivisions of the parameter range in the profile
mcap_nprof <- 10         # Number of starts per profile point
mcap_replicates <- 5
mcap_num_particles_pfilt <- 2000

# IF settings to generate the likelihood surface for each parameter
mcap_num_particles <- 2000
mcap_num_filter_iter <- 100

mcap_cool_type <- "geometric"
mcap_cool_frac <- 0.5
mcap_cool_frac_lastif <- 0.1

# Lambda closer to 1 increases the curvature of the likelihood approximation
mcap_lambda <- c(a0=0.75, a1=0.75, b0=0.75, b1=0.75, sigma=0.75, gamma=0.75)

mcap_ngrid <- c(a0=1000, a1=1000, b0=1000, b1=1000, sigma=1000, gamma=1000)

# =============================================================================
# Function to calculate the MCAP CIs 

mcap <- function(lp, parameter, confidence, lambda, Ngrid) {
  smooth_fit <- loess(lp ~ parameter, span=lambda)
  parameter_grid <- seq(min(parameter), max(parameter), length.out = Ngrid)
  smoothed_loglik <- predict(smooth_fit, newdata=parameter_grid)
  smooth_arg_max <- parameter_grid[which.max(smoothed_loglik)]
  dist <- abs(parameter - smooth_arg_max)
  included <- dist < sort(dist)[trunc(lambda*length(dist))]
  maxdist <- max(dist[included])
  weight <- rep(0, length(parameter))
  weight[included] <- (1 - (dist[included]/maxdist)^3)^3
  quadratic_fit <- lm(lp ~ a + b, weight=weight,
                      data = data.frame(lp=lp, b=parameter, a=-parameter^2))
  b <- unname(coef(quadratic_fit)["b"] )
  a <- unname(coef(quadratic_fit)["a"] )
  m <- vcov(quadratic_fit)
  
  var_b <- m["b", "b"]
  var_a <- m["a", "a"]
  cov_ab <- m["a", "b"]
  
  se_mc_squared <- (1 / (4 * a^2)) * (var_b - (2 * b/a) * cov_ab + (b^2 / a^2) * var_a)
  se_stat_squared <- 1/(2*a)
  se_total_squared <- se_mc_squared + se_stat_squared
  
  delta <- qchisq(confidence, df=1) * ( a * se_mc_squared + 0.5)
  loglik_diff <- max(smoothed_loglik) - smoothed_loglik
  ci <- range(parameter_grid[loglik_diff < delta])
  list(lp=lp, parameter=parameter, confidence=confidence,
       quadratic_fit=quadratic_fit, quadratic_max=b/(2*a),
       smooth_fit=smooth_fit,
       fit=data.frame(parameter=parameter_grid,
                      smoothed=smoothed_loglik,
                      quadratic=predict(quadratic_fit, list(b = parameter_grid, 
                                                            a = -parameter_grid^2))),
       mle=smooth_arg_max, ci=ci, delta=delta,
       se_stat=sqrt(se_stat_squared), se_mc=sqrt(se_mc_squared), 
       se=sqrt(se_total_squared)
  )
}

# =============================================================================
# Generates the model for the profile likelihood calculation

plik <- function(pname, pval0, pval1) {
  registerDoParallel()
  
  bake(file=file.path(cooking_folder, paste(pname, "_mcap.rds", sep="")), {  
    
    desel <- which(names(mle_params) %in% c("loglik", "loglik.se", pname))
    mle_params %>% 
      subset(
        loglik > max(loglik,na.rm=TRUE) - 20,
        select=-desel
      ) %>% 
      melt(id=NULL) %>% 
      daply(~variable, function(x)range(x$value)) -> box
    
    starts <- profileDesign(pname=seq(pval0, pval1, length=mcap_profdes_len),
                            lower=box[,1], upper=box[,2],
                            nprof=mcap_nprof)
    names(starts) <- sub("pname", pname,names(starts))
    
    psizes <- perturb_sizes[names(perturb_sizes) != pname]
    
    foreach(params=iter(starts, "row"),
            .combine=rbind,
            .packages="pomp",
            .options.multicore=list(set.seed=TRUE),
            .options.mpi=list(seed=mcap_plik_seed, chunkSize=1)
    ) %dopar% {
      mf <- mif2(model,
                 params=unlist(params),
                 Np=mcap_num_particles,
                 Nmif=mcap_num_filter_iter,
                 cooling.type=mcap_cool_type,
                 cooling.fraction.50=mcap_cool_frac,
                 rw.sd=do.call(rw.sd, psizes)
      )
      mf <- mif2(mf, 
                 Np=mcap_num_particles,
                 Nmif=mcap_num_filter_iter,
                 cooling.fraction.50=mcap_cool_frac_lastif)
      ll <- logmeanexp(replicate(mcap_replicates, 
                                 logLik(pfilter(mf, Np=mcap_num_particles_pfilt))), se=TRUE)
      data.frame(as.list(coef(mf)), loglik=ll[1], loglik.se=ll[2])
    }
  }) -> pmodel
  
  return(pmodel)
}

# =============================================================================
# Iterate over all the free parameters in the model to calculate their CIs... may take a while!

par_names <- row.names(free_param_box)
for (i in 1:nrow(free_param_box)) {
  name <- par_names[i]  
  print(sprintf("Calculating CI for %s...", name))
  
  row <- free_param_box[i,]
  mdl <- plik(pname=name, pval0=row[1], pval1=row[2])
  
  par_range <- seq(row[1], row[2], length=mcap_profdes_len)
  log_likelihoods <- c()
  for (val in par_range) {
    likelihoods <- subset(mdl, abs(mdl[[name]]-val)<1)$loglik
    if (length(likelihoods) == 0) next
    log_likelihoods <- c(log_likelihoods, max(likelihoods))
  }
  
  x <- mcap(log_likelihoods, par_range, mcap_confidence, mcap_lambda[[i]], mcap_ngrid[[i]])
  if (i == 1) {
    cis <- data.frame("name" = c(name), "x0" = c(x$ci[1]), "x1" = c(x$ci[2]), stringsAsFactors = FALSE)  
  } else {
    cis <- rbind(cis, c(name, x$ci[1], x$ci[2]))  
  }
  print(sprintf("%s %0.2f %0.2f", name, x$ci[1], x$ci[2]))

  ggplot(x$fit, aes(parameter, quadratic)) + geom_line() + 
    geom_vline(xintercept=c(x$ci[1], x$ci[2]), linetype=4, colour='red') +
    geom_point(data = data.frame('parameters'=par_range, 'loglik'=log_likelihoods), 
               aes(parameters, log_likelihoods)) 
  ggsave(file.path(plotting_folder, paste("5-", name, "_ci.pdf", sep="")))
}

write.csv(cis, file=file.path(output_folder, "param_confidence_intervals.csv"), row.names=FALSE, na="")