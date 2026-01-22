


?eDNA_power

t1 = eDNA_power(n = 10, k = 3, theta = 0.6, p = 0.8, s_min = 5)
t2 = eDNA_power(n = 15, k = 2, theta = runif(1000, 0.4, 0.6), p = rbeta(1000, 1.8, 4.2), s_min = 5)




?eDNA_power_grid

g1 = eDNA_power_grid(n_vals = 5:10, k_vals = 1:3, theta = 0.6, p = 0.8, s_min = 5)
g2 = eDNA_power_grid(n_vals = 5:100, k_vals = 1:9, theta = runif(1000, 0.4, 0.6), p = rbeta(1000, 1.8, 4.2), s_min = 5)



?eDNA_power_extract
e2 = eDNA_power_extract(g2, target_power = 0.90)


?eDNA_power_heatmap

# eDNA_power_heatmap(g1, breaks = c(0.25,0.5,0.75))
eDNA_power_heatmap(g2, 5, breaks = c(0.25,0.5,0.75, 0.9, 0.95))




#
# freq_power(10, 2, 0.6, 0.8)
# g1 = freq_power_grid(1:100, 1:6, 0.6, 0.5, 10)
# plot_freq_heatmap(g1)
#
# # beta distrib centered on 0.3 with sd of ~0.2
# theta_prior = rbeta(10000, 1.8, 4.2)
# hist(theta_prior)
# p_prior = rbeta(10000, 1.8, 4.2)
# hist(p_prior)
#
# theta_prior = runif(10000, 0.2, 0.4)
# hist(theta_prior)
# p_prior = runif(10000, 0.4, 0.6)
# hist(p_prior)
#
#
# g2 = bayes_power_grid(1:100, 1:9, 0.6, 0.5, 10,
#                       s_min = 5,
#                      theta_prior = theta_prior,
#                      p_prior = p_prior)
#
# plot_bayes_heatmap(g2)
