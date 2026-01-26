library(eDNAoccupancy)
library(dplyr)
library(tidyr)
data("fungusDetectionData")

fungusDet = occData(fungusDetectionData, siteCol = "site", sampleColName = "sample")
fit = occModel(detectionMats = fungusDet)
posteriorSummary(fit)



# Simulate survey results
n_sites = 1 # locations
n_samples_per_site = 30 # trees
n_tech_reps = 3

theta = runif(1000, 0.2, 0.4)  # Probability that DNA is present in a sample
p = rbeta(1000, 5, 5)      # Probability that a technical replicate detects the DNA if it is
# hist(rbeta(1000, 5, 5))


# df of site, sample (1,2,3 etc), tech_rep (1,2,3 etc) where values are 0/1 for detection/non-detection
df_list = list()
d = 1
for(i in 1:n_sites){
  for(j in 1:n_samples_per_site){
    for(k in 1:n_tech_reps){
      # simulate detection/non-detection data
      dna_present = rbinom(1, 1, sample(theta, size = 1))
      detection = rbinom(1, 1, dna_present * sample(p, size = 1))

      # save to list
      df_list[[d]] = tibble(site = paste0("site_", i),
                       sample = j,
                       tech_rep = paste0("rep_", k),
                       detection = detection)
      d = d + 1
    }
  }
}

surveyDetData = bind_rows(df_list)
surveyDetData = surveyDetData %>%
  pivot_wider(names_from = tech_rep, values_from = detection) %>%
  arrange(site, sample) %>%
  as.data.frame()
surveyDet = occData(surveyDetData, siteCol = "site", sampleColName = "sample")

fit2 = occModel(detectionMats = surveyDet, niter = 11000, niterInterval = 5000)
posteriorSummary(fit2, burnin = 4000, mcError = TRUE)

plotTrace(fit2, paramName = "beta.(Intercept)")
plotTrace(fit2, paramName = "alpha.(Intercept)")

posteriorSummaryOfSiteOccupancy(fit2, burnin = 1000, mcError = TRUE) # psi
posteriorSummaryOfSampleOccupancy(fit2, burnin = 1000, mcError = TRUE) # theta
posteriorSummaryOfDetection(fit2, burnin = 1000, mcError = TRUE) # p



# try msocc next and make simulated data based on N positives desired
library(msocc)

sim <- msocc_sim(M = 1, J = 100, K = 6, psi = 0.5, theta = 0.3, p = 0.6)
mod <- msocc_mod(wide_data = sim$resp,
                 site = list(model = ~ 1, cov_tbl = sim$site),
                 sample = list(model = ~ 1, cov_tbl = sim$sample),
                 rep = list(model = ~ 1, cov_tbl = sim$rep),
                 progress = F)
posterior_summary(mod, print = T)
posterior_summary(mod, level = 'site', print = T)
posterior_summary(mod, level = 'sample', print = T)
posterior_summary(mod, level = 'rep', print = T)



# create simulated df based on n positives and n positive tech reps
library(msocc)
library(dplyr)
library(tidyr)

n_sites = 1
n_samples = 30
n_pos_samples =6 # per site
n_tech_reps = 6
n_pos_tech_reps = 4 # per positive sample

# create df of zeros
df1 = msocc_sim(M = n_sites, J = n_samples, K = n_tech_reps,
                psi = 0.3, theta = 0.5, p = 0.5)

# df1$resp %>%
#   pivot_longer(cols = starts_with("pcr")) %>%
#   group_by(site) %>%
#   summarise(total_positives = sum(value))


# fill in desired n pos samples and n pos tech reps

## select n_pos_samples samples to be positive (per site)
df1$resp =
  df1$resp %>%
    # pivot long pcr replicates
    pivot_longer(cols = starts_with("pcr")) %>%
    group_by(site) %>%
    # select n_pos_samples samples to be positive for first n_pos_samples samples
    mutate(value =
             case_when(
      (sample %in% sample(1:n_samples, n_pos_samples)) &
        (name %in% paste0("pcr", sample(1:n_tech_reps, n_pos_tech_reps))) ~ 1,
      TRUE ~ 0
    )) %>%
    # to int
    mutate(value = as.integer(value)) %>%
    ungroup() %>%
    # pivot wider back to original format
    pivot_wider(names_from = name, values_from = value) %>%
    as.data.frame()


# run model
# sim <- msocc_sim(M = 1, J = 100, K = 6, psi = 0.5, theta = 0.3, p = 0.6)

mod <- msocc_mod(wide_data = df1$resp,
                 site = list(model = ~ 1, cov_tbl = df1$site),
                 sample = list(model = ~ 1, cov_tbl = df1$sample),
                 rep = list(model = ~ 1, cov_tbl = df1$rep),
                 progress = T, num.mcmc = 1e4)
posterior_summary(mod, print = T)
posterior_summary(mod, level = 'site', print = T)
posterior_summary(mod, level = 'sample', print = T)
posterior_summary(mod, level = 'rep', print = T)


# vary n total, n pos samples, n pos tech reps and see how CIs change

sim_data = function(n_sites,
                    n_samples, n_pos_samples,
                    n_tech_reps, n_pos_tech_reps
                    ){

  # generated simulated data template
  sim = msocc_sim(M = n_sites, J = n_samples, K = n_tech_reps,
                  psi = 0, theta = 0, p = 0)

  # fill in desired n pos samples and n pos tech reps
  sim$resp =
    sim$resp %>%
    # pivot long pcr replicates
    pivot_longer(cols = starts_with("pcr")) %>%
    group_by(site) %>%
    # select n_pos_samples samples to be positive for first n_pos_samples samples
    mutate(value =
             case_when(
               (sample %in% sample(1:n_samples, n_pos_samples)) &
                 (name %in% paste0("pcr", sample(1:n_tech_reps, n_pos_tech_reps))) ~ 1,
               TRUE ~ 0
             )) %>%
    # to int
    mutate(value = as.integer(value)) %>%
    ungroup() %>%
    # pivot wider back to original format
    pivot_wider(names_from = name, values_from = value) %>%
    as.data.frame()

  return(sim)

}

# create df grid of different scenarios
grid = expand.grid(n_sites = 1,
                   n_samples = c(30,50,100),
                   n_pos_samples = c(3,10),
                   n_tech_reps = c(6),
                   n_pos_tech_reps = c(1,3, 5))
## rm if more positive samples than total samples OR
## morepositive technical reps than total technical reps
grid = grid %>%
  filter(n_pos_samples <= n_samples) %>%
  filter(n_pos_tech_reps <= n_tech_reps)

# run sim_data for each scenario
sim_list = list()
for(i in 1:nrow(grid)){
  sim_list[[i]] = sim_data(n_sites = grid$n_sites[i],
                             n_samples = grid$n_samples[i],
                             n_pos_samples = grid$n_pos_samples[i],
                             n_tech_reps = grid$n_tech_reps[i],
                             n_pos_tech_reps = grid$n_pos_tech_reps[i])
}

# run model for each scenario and save posterior summaries
posterior_sample_list = list()
for(i in 1:length(sim_list)){
  mod <- msocc_mod(wide_data = sim_list[[i]]$resp,
                   site = list(model = ~ 1, cov_tbl = sim_list[[i]]$site),
                   sample = list(model = ~ 1, cov_tbl = sim_list[[i]]$sample),
                   rep = list(model = ~ 1, cov_tbl = sim_list[[i]]$rep),
                   progress = F, num.mcmc = 1e4)

  # get posterior summary for sample occupancy
  posterior_sample_list[[i]] = posterior_summary(mod, level = 'sample', print = F)
  ## add in scenario info
  posterior_sample_list[[i]]$n_samples = grid$n_samples[i]
  posterior_sample_list[[i]]$n_pos_samples = grid$n_pos_samples[i]
  posterior_sample_list[[i]]$n_tech_reps = grid$n_tech_reps[i]
  posterior_sample_list[[i]]$n_pos_tech_reps = grid$n_pos_tech_reps[i]

  message(paste0("Completed scenario ", i, " of ", nrow(grid)))
}

# to one df
posterior_sample_df = bind_rows(posterior_sample_list)

# get range
posterior_sample_df$range = posterior_sample_df$`0.975` - posterior_sample_df$`0.025`

library(ggplot2)

# plot where y axis is range, x is n total samples, color is n pos samples, facet is n pos tech reps
ggplot(posterior_sample_df, aes(x = as.factor(n_samples), y = range,
                                 color = as.factor(n_pos_samples))) +
  geom_point() +
  geom_line(aes(group = n_pos_samples)) +
  facet_wrap(~n_pos_tech_reps) +
  labs(x = "Total number of samples",
       y = "Range of 95% credible interval for sample occupancy",
       color = "Number of positive samples") +
  theme_minimal()
