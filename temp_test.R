

freq_power(10, 2, 0.6, 0.8)
g1 = freq_power_grid(1:100, 1:6, 0.6, 0.5, 10)
plot_freq_heatmap(g1)

# beta distrib centered on 0.3 with sd of ~0.2
theta_prior = rbeta(10000, 1.8, 4.2)
hist(theta_prior)
p_prior = rbeta(10000, 1.8, 4.2)
hist(p_prior)

theta_prior = runif(10000, 0.2, 0.4)
hist(theta_prior)
p_prior = runif(10000, 0.4, 0.6)
hist(p_prior)


g2 = bayes_power_grid(1:100, 1:9, 0.6, 0.5, 10,
                      s_min = 5,
                     theta_prior = theta_prior,
                     p_prior = p_prior)

plot_bayes_heatmap(g2)
