library(BayesianDesignAna)
ratios <- c(Placebo = 1, Dose1 = 2, Dose2 = 2, Dose3 = 1)
n <- allocate_arms(180, ratios)
cap <- borrowing_limits(n_control = n[['Placebo']], historical_ess = 65,
                        max_fraction = 0.4, max_absolute = 25)
print(cap)
prior <- borrowing_prior(historical_rate = .3, ess = cap$allowed_ess)
print(analyze_binary(responses = c(Placebo = 9, Dose1 = 25, Dose2 = 34, Dose3 = 17),
                     n = n, control_prior = c(prior$alpha, prior$beta), seed = 1))
rates <- c(Placebo = .3, Dose1 = .4, Dose2 = .5, Dose3 = .55)
print(simulate_power(n, rates, borrow_ess = cap$allowed_ess,
                     historical_rate = .3, n_sim = 100, n_post = 1000))
print(find_sample_size(seq(120, 240, 40), ratios, rates, target = 'any',
                       borrow_ess = cap$allowed_ess, historical_rate = .3,
                       n_sim = 100, n_post = 1000)$results)
