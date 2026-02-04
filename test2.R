
# ---------------------------------------------------------------------------- #

# sim data fun
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



# ---------------------------------------------------------------------------- #

# fun to run power analysis
pa = function(n_sites,
              n_samples,
              n_pos_samples,
              n_tech_reps,
              n_pos_tech_reps){

  # create df grid of different scenarios
  grid = expand.grid(n_sites = n_sites,
                     n_samples = n_samples,
                     n_pos_samples = n_pos_samples,
                     n_tech_reps = n_tech_reps,
                     n_pos_tech_reps = n_pos_tech_reps)

  ## rm if more positive samples than total samples OR
  ## morepositive technical reps than total technical reps
  grid = grid %>%
    filter(n_pos_samples <= n_samples) %>%
    filter(n_pos_tech_reps <= n_tech_reps)


  ## get sim data for each scenario
  sim_list = list()
  for(i in 1:nrow(grid)){
    sim_list[[i]] = sim_data(n_sites = grid$n_sites[i],
                             n_samples = grid$n_samples[i],
                             n_pos_samples = grid$n_pos_samples[i],
                             n_tech_reps = grid$n_tech_reps[i],
                             n_pos_tech_reps = grid$n_pos_tech_reps[i])
  }


  # run model for each scenario and save posterior summaries
  mod_list = list()
  for(i in 1:length(sim_list)){
    mod_list[[i]] <- msocc_mod(wide_data = sim_list[[i]]$resp,
                               site = list(model = ~ 1, cov_tbl = sim_list[[i]]$site),
                               sample = list(model = ~ 1, cov_tbl = sim_list[[i]]$sample),
                               rep = list(model = ~ 1, cov_tbl = sim_list[[i]]$rep),
                               progress = F, num.mcmc = 5e4)
    message(paste0("Completed scenario ", i, " of ", nrow(grid)))
  }


  # extract info, takes a little depending number of scenarios
  posterior_sample_list = list()
  for(i in 1:length(mod_list)){

    # get posterior summaries
    posterior_sample_list[[i]] =
      bind_cols(
        # get posterior summary for sample detection
        posterior_summary(mod_list[[i]], level = 'sample', print = F)[1,] %>% # assuming no site difference (ie covariates)
          select(median, mean, `0.025`, `0.975`) %>%
          rename_with(.cols = everything(),
                      .fn = ~ paste0("sample_", .)) ,

        # get posterior summary for rep detection
        posterior_summary(mod_list[[i]], level = 'rep', print = F)[1,] %>% # assuming no site difference (ie covariates)
          select(median, mean, `0.025`, `0.975`) %>%
          rename_with(.cols = everything(),
                      .fn = ~ paste0("rep_", .))
      )

    ## add in scenario info
    posterior_sample_list[[i]]$n_sites = grid$n_sites[i]
    posterior_sample_list[[i]]$n_samples = grid$n_samples[i]
    posterior_sample_list[[i]]$n_pos_samples = grid$n_pos_samples[i]
    posterior_sample_list[[i]]$n_tech_reps = grid$n_tech_reps[i]
    posterior_sample_list[[i]]$n_pos_tech_reps = grid$n_pos_tech_reps[i]

    message(paste0("Extracted info for scenario ", i, " of ", nrow(grid)))

  }


  # to one df
  posterior_sample_df = bind_rows(posterior_sample_list)

  # get range
  posterior_sample_df$sample_range = posterior_sample_df$`sample_0.975` - posterior_sample_df$`sample_0.025`
  posterior_sample_df$rep_range = posterior_sample_df$`rep_0.975` - posterior_sample_df$`rep_0.025`

  return(posterior_sample_df)
}


# ---------------------------------------------------------------------------- #

library(msocc)
library(dplyr)
library(tidyr)
library(ggplot2)


# run pa
df1 = pa(n_sites = 1,
            n_samples = c(15,30,50, 100),
            n_pos_samples = c(1,3,5,10, 15, 22, 30, 40, 50, 75, 100),
            n_tech_reps = c(6),
            n_pos_tech_reps = c(1,3))


# saveRDS(df1, file = "sample_posterior_sample_df.rds")
# df1 = readRDS("sample_posterior_sample_df.rds")


ptr.labs <- paste0(unique(df1$n_pos_tech_reps), " pos TRs")
names(ptr.labs) <- unique(df1$n_pos_tech_reps)

tr.labs = paste0(unique(df1$n_tech_reps), " TRs")
names(tr.labs) <- unique(df1$n_tech_reps)

# line plot theta v2
df1 %>%
  mutate(prop_pos = n_pos_samples / n_samples) %>%
  ggplot(aes(x = prop_pos, y = sample_range, group = n_samples,
             color = as.factor(n_samples))) +
  geom_point() +
  geom_line() +
  # add h line at some ideal value
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50") +
  # scale_color_brewer(palette = "Reds") +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_color_viridis_d() +
  facet_grid(n_tech_reps ~ n_pos_tech_reps,
             labeller= labeller(n_pos_tech_reps = ptr.labs, n_tech_reps = tr.labs)) +
  labs(color = "Total number of samples",
       y = "Range of 95% credible interval for theta",
       x = "Proportion of positive samples") +
  theme_bw()


df1 %>%
  filter(n_tech_reps == 6) %>%
  ggplot(aes(x = as.factor(n_pos_samples),
             y = sample_median, group = n_samples, color = n_samples)) +
  geom_point(position=position_dodge(width=0.5)) +
  # geom_ribbon(aes(ymin = sample_0.025, ymax = sample_0.975, fill = n_samples), alpha = 0.2) +
  geom_errorbar(aes(ymin = sample_0.025, ymax = sample_0.975), width = 0.2,
                position=position_dodge(width=0.5)) +
  scale_color_viridis_b() +
  labs(x = "Number of positive samples",
       y = "Range of 95% CI for sample occupancy") +
  facet_grid(~ n_pos_tech_reps,
             labeller= labeller(n_pos_tech_reps = ptr.labs, n_tech_reps = tr.labs)
             ) +
  theme_bw()




# library(eDNAoccupancy)
# library(dplyr)
# library(tidyr)
# data("fungusDetectionData")
#
# fungusDet = occData(fungusDetectionData, siteCol = "site", sampleColName = "sample")
# fit = occModel(detectionMats = fungusDet)
# posteriorSummary(fit)
#
#
#
# # Simulate survey results
# n_sites = 1 # locations
# n_samples_per_site = 30 # trees
# n_tech_reps = 3
#
# theta = runif(1000, 0.2, 0.4)  # Probability that DNA is present in a sample
# p = rbeta(1000, 5, 5)      # Probability that a technical replicate detects the DNA if it is
# # hist(rbeta(1000, 5, 5))
#
#
# # df of site, sample (1,2,3 etc), tech_rep (1,2,3 etc) where values are 0/1 for detection/non-detection
# df_list = list()
# d = 1
# for(i in 1:n_sites){
#   for(j in 1:n_samples_per_site){
#     for(k in 1:n_tech_reps){
#       # simulate detection/non-detection data
#       dna_present = rbinom(1, 1, sample(theta, size = 1))
#       detection = rbinom(1, 1, dna_present * sample(p, size = 1))
#
#       # save to list
#       df_list[[d]] = tibble(site = paste0("site_", i),
#                        sample = j,
#                        tech_rep = paste0("rep_", k),
#                        detection = detection)
#       d = d + 1
#     }
#   }
# }
#
# surveyDetData = bind_rows(df_list)
# surveyDetData = surveyDetData %>%
#   pivot_wider(names_from = tech_rep, values_from = detection) %>%
#   arrange(site, sample) %>%
#   as.data.frame()
# surveyDet = occData(surveyDetData, siteCol = "site", sampleColName = "sample")
#
# fit2 = occModel(detectionMats = surveyDet, niter = 11000, niterInterval = 5000)
# posteriorSummary(fit2, burnin = 4000, mcError = TRUE)
#
# plotTrace(fit2, paramName = "beta.(Intercept)")
# plotTrace(fit2, paramName = "alpha.(Intercept)")
#
# posteriorSummaryOfSiteOccupancy(fit2, burnin = 1000, mcError = TRUE) # psi
# posteriorSummaryOfSampleOccupancy(fit2, burnin = 1000, mcError = TRUE) # theta
# posteriorSummaryOfDetection(fit2, burnin = 1000, mcError = TRUE) # p
#
#
#
# # try msocc next and make simulated data based on N positives desired
# library(msocc)
#
# sim <- msocc_sim(M = 1, J = 100, K = 6, psi = 0.5, theta = 0.3, p = 0.6)
# mod <- msocc_mod(wide_data = sim$resp,
#                  site = list(model = ~ 1, cov_tbl = sim$site),
#                  sample = list(model = ~ 1, cov_tbl = sim$sample),
#                  rep = list(model = ~ 1, cov_tbl = sim$rep),
#                  progress = F)
# posterior_summary(mod, print = T)
# posterior_summary(mod, level = 'site', print = T)
# posterior_summary(mod, level = 'sample', print = T)
# posterior_summary(mod, level = 'rep', print = T)
#
#
#
# # create simulated df based on n positives and n positive tech reps
# library(msocc)
# library(dplyr)
# library(tidyr)
#
# n_sites = 2
# n_samples = 30
# n_pos_samples =6 # per site
# n_tech_reps = 6
# n_pos_tech_reps = 4 # per positive sample
#
# # create df of zeros
# df1 = msocc_sim(M = n_sites, J = n_samples, K = n_tech_reps,
#                 psi = 0.3, theta = 0.5, p = 0.5)
#
# # df1$resp %>%
# #   pivot_longer(cols = starts_with("pcr")) %>%
# #   group_by(site) %>%
# #   summarise(total_positives = sum(value))
#
#
# # fill in desired n pos samples and n pos tech reps
#
# ## select n_pos_samples samples to be positive (per site)
# df1$resp =
#   df1$resp %>%
#     # pivot long pcr replicates
#     pivot_longer(cols = starts_with("pcr")) %>%
#     group_by(site) %>%
#     # select n_pos_samples samples to be positive for first n_pos_samples samples
#     mutate(value =
#              case_when(
#       (sample %in% sample(1:n_samples, n_pos_samples)) &
#         (name %in% paste0("pcr", sample(1:n_tech_reps, n_pos_tech_reps))) ~ 1,
#       TRUE ~ 0
#     )) %>%
#     # to int
#     mutate(value = as.integer(value)) %>%
#     ungroup() %>%
#     # pivot wider back to original format
#     pivot_wider(names_from = name, values_from = value) %>%
#     as.data.frame()
#
#
# # run model
# # sim <- msocc_sim(M = 1, J = 100, K = 6, psi = 0.5, theta = 0.3, p = 0.6)
#
# mod <- msocc_mod(wide_data = df1$resp,
#                  site = list(model = ~ 1, cov_tbl = df1$site),
#                  sample = list(model = ~ 1, cov_tbl = df1$sample),
#                  rep = list(model = ~ 1, cov_tbl = df1$rep),
#                  progress = T, num.mcmc = 1e4)
# posterior_summary(mod, print = T)
# posterior_summary(mod, level = 'site', print = T)
# posterior_summary(mod, level = 'sample', print = T)
# posterior_summary(mod, level = 'rep', print = T)
#
#
#
#
#
# # ---------------------------------------------------------------------------- #
#
#
#
# # vary n total, n pos samples, n pos tech reps and see how CIs change
# sim_data = function(n_sites,
#                     n_samples, n_pos_samples,
#                     n_tech_reps, n_pos_tech_reps
#                     ){
#
#   # generated simulated data template
#   sim = msocc_sim(M = n_sites, J = n_samples, K = n_tech_reps,
#                   psi = 0, theta = 0, p = 0)
#
#   # fill in desired n pos samples and n pos tech reps
#   sim$resp =
#     sim$resp %>%
#     # pivot long pcr replicates
#     pivot_longer(cols = starts_with("pcr")) %>%
#     group_by(site) %>%
#     # select n_pos_samples samples to be positive for first n_pos_samples samples
#     mutate(value =
#              case_when(
#                (sample %in% sample(1:n_samples, n_pos_samples)) &
#                  (name %in% paste0("pcr", sample(1:n_tech_reps, n_pos_tech_reps))) ~ 1,
#                TRUE ~ 0
#              )) %>%
#     # to int
#     mutate(value = as.integer(value)) %>%
#     ungroup() %>%
#     # pivot wider back to original format
#     pivot_wider(names_from = name, values_from = value) %>%
#     as.data.frame()
#
#   return(sim)
#
# }
#
# # create df grid of different scenarios
# grid = expand.grid(n_sites = c(1,2),
#                    n_samples = c(15,30,50),
#                    n_pos_samples = c(1,3,5,10),
#                    n_tech_reps = c(6),
#                    n_pos_tech_reps = c(1,3, 5))
#
# ## rm if more positive samples than total samples OR
# ## morepositive technical reps than total technical reps
# grid = grid %>%
#   filter(n_pos_samples <= n_samples) %>%
#   filter(n_pos_tech_reps <= n_tech_reps)
#
# # run sim_data for each scenario
# sim_list = list()
# for(i in 1:nrow(grid)){
#   sim_list[[i]] = sim_data(n_sites = grid$n_sites[i],
#                              n_samples = grid$n_samples[i],
#                              n_pos_samples = grid$n_pos_samples[i],
#                              n_tech_reps = grid$n_tech_reps[i],
#                              n_pos_tech_reps = grid$n_pos_tech_reps[i])
# }
#
# # run model for each scenario and save posterior summaries
# mod_list = list()
# for(i in 1:length(sim_list)){
#   mod_list[[i]] <- msocc_mod(wide_data = sim_list[[i]]$resp,
#                      site = list(model = ~ 1, cov_tbl = sim_list[[i]]$site),
#                      sample = list(model = ~ 1, cov_tbl = sim_list[[i]]$sample),
#                      rep = list(model = ~ 1, cov_tbl = sim_list[[i]]$rep),
#                      progress = F, num.mcmc = 5e4)
#   message(paste0("Completed scenario ", i, " of ", nrow(grid)))
# }
#
#
# # extract info, takes a little depending number of scenarios
# posterior_sample_list = list()
# for(i in 1:length(mod_list)){
#
#   # get posterior summaries
#   posterior_sample_list[[i]] =
#     bind_cols(
#       # get posterior summary for sample detection
#       posterior_summary(mod_list[[i]], level = 'sample', print = F)[1,] %>% # assuming no site difference (ie covariates)
#         select(median, mean, `0.025`, `0.975`) %>%
#         rename_with(.cols = everything(),
#                     .fn = ~ paste0("sample_", .)) ,
#
#       # get posterior summary for rep detection
#       posterior_summary(mod_list[[i]], level = 'rep', print = F)[1,] %>% # assuming no site difference (ie covariates)
#         select(median, mean, `0.025`, `0.975`) %>%
#         rename_with(.cols = everything(),
#                     .fn = ~ paste0("rep_", .))
#     )
#
#   ## add in scenario info
#   posterior_sample_list[[i]]$n_sites = grid$n_sites[i]
#   posterior_sample_list[[i]]$n_samples = grid$n_samples[i]
#   posterior_sample_list[[i]]$n_pos_samples = grid$n_pos_samples[i]
#   posterior_sample_list[[i]]$n_tech_reps = grid$n_tech_reps[i]
#   posterior_sample_list[[i]]$n_pos_tech_reps = grid$n_pos_tech_reps[i]
#
#   message(paste0("Extracted info for scenario ", i, " of ", nrow(grid)))
#
# }
#
#
# # to one df
# posterior_sample_df = bind_rows(posterior_sample_list)
#
# # get range
# posterior_sample_df$sample_range = posterior_sample_df$`sample_0.975` - posterior_sample_df$`sample_0.025`
# posterior_sample_df$rep_range = posterior_sample_df$`rep_0.975` - posterior_sample_df$`rep_0.025`
#
# # # save for future use
# # saveRDS(posterior_sample_df, file = "sample_posterior_sample_df.rds")
# # posterior_sample_df = readRDS("sample_posterior_sample_df.rds")
#
# library(ggplot2)
#
#
# # Define new facet labels
# ptr.labs <- paste0("# pos TRs = ", unique(posterior_sample_df$n_pos_tech_reps))
# names(ptr.labs) <- unique(posterior_sample_df$n_pos_tech_reps)
#
# tr.labs = paste0("# TRs = ", unique(posterior_sample_df$n_tech_reps))
# names(tr.labs) <- unique(posterior_sample_df$n_tech_reps)
#
# # # line plot theta
# # posterior_sample_df %>%
# #   ggplot(aes(x = n_samples, y = sample_range, group = n_pos_samples,
# #              color = as.factor(n_pos_samples))) +
# #   geom_point() +
# #   geom_line() +
# #   # add h line at some ideal value
# #   geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50") +
# #   # scale_color_brewer(palette = "Reds") +
# #   scale_color_viridis_d() +
# #   facet_grid(n_pos_tech_reps ~ n_tech_reps,
# #              labeller= labeller(n_pos_tech_reps = ptr.labs, n_tech_reps = tr.labs)) +
# #   labs(color = "Number of positive samples",
# #        y = "Range of 95% credible interval for theta",
# #        x = "Total number of samples") +
# #   theme_bw()
#
# # line plot theta v2
# posterior_sample_df %>%
#   ggplot(aes(x = n_pos_samples, y = sample_range, group = n_samples,
#              color = as.factor(n_samples))) +
#   geom_point() +
#   geom_line() +
#   # add h line at some ideal value
#   geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50") +
#   # scale_color_brewer(palette = "Reds") +
#   scale_x_continuous(breaks = scales::pretty_breaks()) +
#   scale_color_viridis_d() +
#   facet_grid(n_pos_tech_reps ~ n_tech_reps,
#              labeller= labeller(n_pos_tech_reps = ptr.labs, n_tech_reps = tr.labs)) +
#   labs(color = "Total number of samples",
#        y = "Range of 95% credible interval for theta",
#        x = "Number of positive samples") +
#   theme_bw()
#
#
#
#
# # line plot p
# posterior_sample_df %>%
#   ggplot(aes(x = n_samples, y = rep_range, group = n_pos_samples,
#              color = as.factor(n_pos_samples))) +
#   geom_point() +
#   geom_line() +
#   # add h line at some ideal value
#   geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50") +
#   scale_color_brewer(palette = "Reds") +
#   facet_grid(n_pos_tech_reps ~ n_tech_reps,
#              # labeller= labeller(n_pos_tech_reps = ptr.labs, n_tech_reps = tr.labs)
#              ) +
#   labs(color = "Number of positive samples",
#        y = "Range of 95% credible interval for p",
#        x = "Total number of samples") +
#   theme_bw()
#
#
# ### notes,
# # n sites don't matter really, it's equal to n positives
# # why is less positives have lower ci range for theta? ... perhaps bc theta itself is smaller so smllare bands
# # next step to plot theta CIs at 3posTR, vs nPosSamples
#
# posterior_sample_df %>%
#   filter(n_sites == 1, n_tech_reps == 6) %>%
#   ggplot(aes(x = as.factor(n_pos_samples),
#                  y = sample_median)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = sample_0.025, ymax = sample_0.975), width = 0.2) +
#   labs(x = "Number of positive samples",
#        y = "Range of 95% CI for sample occupancy") +
#   facet_grid(n_samples ~ n_pos_tech_reps) +
#   theme_bw()
#
# posterior_sample_df %>%
#   filter(n_tech_reps == 6) %>%
#   ggplot(aes(x = as.factor(n_pos_samples),
#              y = sample_median, group = n_samples, color = n_samples)) +
#   geom_point(position=position_dodge(width=0.5)) +
#   # geom_ribbon(aes(ymin = sample_0.025, ymax = sample_0.975, fill = n_samples), alpha = 0.2) +
#   geom_errorbar(aes(ymin = sample_0.025, ymax = sample_0.975), width = 0.2,
#                 position=position_dodge(width=0.5)) +
#   scale_color_viridis_b() +
#   labs(x = "Number of positive samples",
#        y = "Range of 95% CI for sample occupancy") +
#   facet_grid(~ n_pos_tech_reps) +
#   theme_bw()
#
#
# posterior_sample_df %>%
#   filter(n_tech_reps == 6, n_pos_tech_reps == 3, n_pos_samples == 1) %>%
#   View()
#
#
#
# # # heatmap
# # posterior_sample_df %>%
# #   ggplot(aes(x = n_pos_samples, y = n_samples, fill = sample_range)) +
# #   geom_tile() +
# #   facet_grid(n_pos_tech_reps ~ n_tech_reps) +
# #   scale_fill_fermenter(breaks = c(0.1, 0.2, 0.3, 0.5, 0.75)) +
# #   labs(x = "Number of positive samples",
# #        y = "Total number of samples",
# #        fill = "Range of 95% CI for sample occupancy") +
# #   theme_bw()
