args = commandArgs(trailingOnly=TRUE)

main_folder <- "."
prop_file <- "default.properties"

if (0 < length(args)) {
  main_folder <- args[1]
}

if (1 < length(args)) {
  prop_file <- file.path(main_folder, args[2])
}

# =============================================================================
# Libraries

library(doRNG)
library(foreach)
library(doParallel)

library(properties)

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
# Folders and properties

prop <- read.properties(prop_file)
output_name <- prop$output_name
output_folder <- file.path(main_folder, output_name) 
if (!dir.exists(output_folder)) dir.create(output_folder)
cooking_folder <- file.path(main_folder, output_name, "bake")
if (!dir.exists(cooking_folder)) dir.create(cooking_folder)
plotting_folder <- file.path(main_folder, output_name, "plots")
if (!dir.exists(plotting_folder)) dir.create(plotting_folder)
code_folder <- file.path(main_folder, "code")
if (!dir.exists(code_folder)) dir.create(code_folder)
file_name <- "snippets"

# =============================================================================
# Observed data

csv_table <- read.csv(file.path(main_folder, "case_counts.csv"))
obs_data <- data.frame(day = csv_table$Time, cases = csv_table$Count)
head(obs_data)

ggplot(data=obs_data, aes(x=day, y=cases, group=1)) + geom_line()
ggsave(file.path(plotting_folder, "0-observed_data.pdf"))

# =============================================================================
# General Parameters

pop_size <- as.integer(prop$pop_size)
exp0 <- as.integer(prop$exp0)
inf0 <- as.integer(prop$inf0)
rec0 <- as.integer(prop$rec0)

time_unit <- prop$time_unit
time_start <- max(as.numeric(prop$time_start), obs_data$day[1])
time_end <- min(as.numeric(prop$time_end), max(obs_data$day))
time_step <- as.numeric(prop$time_step)

# Subset the observed data to the region of interest
obs_data <- subset(obs_data, (time_start <= day & day <= time_end))
#ggplot(data=obs_data, aes(x=day, y=cases, group=1)) + geom_line()

free_param_names <- c("a00", "a01", "a10", "a11")
free_param_box <- rbind(
  a00 = c(as.numeric(prop$bounds_a00_min), as.numeric(prop$bounds_a00_max)),
  a01 = c(as.numeric(prop$bounds_a01_min), as.numeric(prop$bounds_a01_max)),
  a10 = c(as.numeric(prop$bounds_a10_min), as.numeric(prop$bounds_a10_max)),
  a11 = c(as.numeric(prop$bounds_a11_min), as.numeric(prop$bounds_a11_max))
)

log_trans_params <- c("a00", "a01", "a10", "a11")
logit_trans_params <- c()

fixed_param_names <- c("pop", "S_0", "E_0", "I_0", "R_0", "rho", "sigma", "gamma")
fixed_param_values <- c(pop=pop_size, S_0=1-(exp0+inf0+rec0)/pop_size, E_0=exp0/pop_size, I_0=inf0/pop_size, R_0=rec0/pop_size, rho=1.0, sigma=0.26, gamma=0.6)
all_param_names <- c(free_param_names, fixed_param_names)

# Random seeds, keep unchanged to ensure reproducibilty of results
test_mle_seed <- as.integer(prop$test_mle_seed)
full_mle_seed <- as.integer(prop$full_mle_seed)
global_search_seed <- as.integer(prop$global_search_seed)
test_sim_seed <- as.integer(prop$test_sim_seed)
full_sim_seed <- as.integer(prop$full_sim_seed)
mcap_plik_seed <- as.integer(prop$mcap_plik_seed)

# =============================================================================
# Csnippets defining the SEIR model

# pomp C API:
# https://kingaa.github.io/pomp/vignettes/C_API.html

sir_step <- Csnippet("
  double beta;
  double foi;
  double rate[3], trans[3];

  beta = calc_beta(t, a00, a01, a10, a11);

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

extra <- Csnippet(gsub("MAIN_FOLDER", main_folder, "
double calc_beta(double td, double a00, double a01, double a10, double a11) {
  static int *indices = NULL;
  static double *contacts = NULL;
  static int max_t = 0;
  static int num_v = 0;

  if (indices == NULL) {
    FILE *file;

    file = fopen(\"MAIN_FOLDER/indices\", \"r\");

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

    file = fopen(\"MAIN_FOLDER/contacts\", \"r\");
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
      double p = a00 + a01 * y + a10 *x + a11 * x * y;
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
"))

# =============================================================================
# The documentation of the pomp object
# https://kingaa.github.io/pomp/manual/pomp.html

obs_data %>% 
  pomp(t0 = time_start,
       time = time_unit,
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
         log=log_trans_params,
         logit=logit_trans_params),
       paramnames = c(free_param_names, fixed_param_names),
       #compile=FALSE,
       verbose = TRUE
  ) -> model

# =============================================================================
# IF parameters, see more details in the manual and tutorial:
# https://kingaa.github.io/pomp/manual/mif2.html
# https://kingaa.github.io/sbied/mif/mif.html#choosing-the-algorithmic-settings-for-if2

num_test_runs <- as.integer(prop$num_test_runs)

num_guesses <- as.integer(prop$num_guesses)   
num_filter_iter <- as.integer(prop$num_filter_iter)
num_particles <- as.integer(prop$num_particles)
num_replicates <- as.integer(prop$num_replicates)

perturb_sizes <- list(a00=as.numeric(prop$perturb_size_a00),
                      a01=as.numeric(prop$perturb_size_a01),
                      a10=as.numeric(prop$perturb_size_a10),
                      a11=as.numeric(prop$perturb_size_a11))

cool_frac <- as.numeric(prop$cool_frac)
cool_type <- prop$cool_type

# Variables to use in the scatterplot matrix showing the result of the IF search
pair_vars <- ~loglik+a00+a01+a10+a11

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

pdf(file=file.path(plotting_folder, "2-mle_global_search.pdf"))
pairs(pair_vars, data=all, col=ifelse(all$type=="guess", grey(0.5), "red"), pch=16)
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
# Calculate mean and standard deviatio of MLEs over the top MLEs
mle_table <- read.csv(file.path(output_folder, "param_mle_global_search.csv"))
mle_sorted <- mle_table[order(-mle_table$loglik),]

n <- as.integer(prop$num_mean_sd)
mle_mean_sd <- NULL
for (p in free_param_names) {
  p_mean <- mean(mle_sorted[[p]][1:n])
  p_sd <- sd(mle_sorted[[p]][1:n])
  if (is.null(mle_mean_sd)) {
    mle_mean_sd <- data.frame("name" = c(p), "mean" = c(p_mean), "sd" = c(p_sd), stringsAsFactors = FALSE)  
  } else {
    mle_mean_sd <- rbind(mle_mean_sd, c(p, p_mean, p_sd))
  }  
  print(sprintf("%s %0.2f %0.2f", p, p_mean, p_sd))
}
write.csv(mle_mean_sd, file=file.path(output_folder, "param_mean_sdev.csv"), row.names=FALSE, na="")

# =============================================================================
# Simulation parameters

num_sims <- as.integer(prop$num_sims)

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

cumulative_curve <- function(dat, t0, t1) {
  total_sum <- 0
  daily_sum <- c()
  for (i in t0:t1) {
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
csim <- cumulative_curve(med_sim, time_start + 1, time_end)

# Getting the cumulative curve for observed data
cobs <- cumulative_curve(obs_data, time_start + 1, time_end)

df <- data.frame('day' = seq(time_start + 1, time_end), 
                 'obs_data' = cobs,
                 'sim_data' = csim)

ggplot(df, aes(day)) + 
  geom_line(aes(y = obs_data, colour = "Real Data")) + 
  geom_line(aes(y = sim_data, colour = "Simulated")) +
  ylab('Cumulative Number of Cases')

ggsave(file.path(plotting_folder, "4-cumulative_cases.pdf"))

# =============================================================================
# MCAP settings

calculate_cis <- as.logical(prop$calculate_cis)

mcap_confidence <- as.numeric(prop$mcap_confidence)

mcap_profdes_len <- as.integer(prop$mcap_profdes_len)
mcap_nprof <- as.integer(prop$mcap_nprof)
mcap_replicates <- as.integer(prop$mcap_replicates)
mcap_num_particles_pfilt <- as.integer(prop$mcap_num_particles_pfilt)

mcap_num_particles <- as.integer(prop$mcap_num_particles)
mcap_num_filter_iter <- as.integer(prop$mcap_num_filter_iter)

mcap_cool_type <- prop$mcap_cool_type
mcap_cool_frac <- as.numeric(prop$mcap_cool_frac)
mcap_cool_frac_lastif <- as.numeric(prop$mcap_cool_frac_lastif)

mcap_lambda <- c(a00=as.numeric(prop$mcap_lambda_a00),
                 a01=as.numeric(prop$mcap_lambda_a01),
                 a10=as.numeric(prop$mcap_lambda_a10),
                 a11=as.numeric(prop$mcap_lambda_a11))

mcap_ngrid <- c(a00=as.numeric(prop$mcap_ngrid_a00),
                a01=as.numeric(prop$mcap_ngrid_a01),
                a10=as.numeric(prop$mcap_ngrid_a10),
                a11=as.numeric(prop$mcap_ngrid_a11))

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
  
  bake(file=file.path(cooking_folder, paste("mcap_search_", pname, ".rds", sep="")), {  
    
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

if (calculate_cis) {
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
}