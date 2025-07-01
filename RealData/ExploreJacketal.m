

DAT = load("Jackson_etal_data.mat");

histogram2(DAT.ts_dat(:,2),DAT.ts_dat(:,3))


figure;
subplot(1,3,1)
scatter(DAT.fit_sim(:,1),DAT.fit_sim(:,4),1,'blue','filled')
title("Estimated Fertility")

subplot(1,3,2)
scatter(DAT.fit_sim(:,1),DAT.fit_sim(:,5),1,'blue','filled')
title("Estimated Survival")


% Calculate total population
years = min(DAT.ts_dat(:,3)):max(DAT.ts_dat(:,3));
tot_pop = zeros(length(years),1);
for k = 1:length(years)
    tot_pop(k) = sum(DAT.ts_dat(:,3) == years(k));
end

subplot(1,3,3)
plot(years,tot_pop)
title("Total Population")



%% === Process and Fit Gompertz Model to Cumulative Survival ===

% Step 1: Extract and average survival to next age
ages_raw = DAT.fit_sim(:,1);
surv_probs_raw = DAT.fit_sim(:,5);
Mort_probs = 1-surv_probs_raw;

% Get unique ages and mean survival per age class
unique_ages = unique(ages_raw);
mean_surv = arrayfun(@(a) mean(surv_probs_raw(ages_raw == a)), unique_ages);

% Step 2: Convert to cumulative survival: S(age) = prod_{i=0}^{age-1} s(i)
cum_surv = ones(size(mean_surv));
for i = 2:length(mean_surv)
    cum_surv(i) = cum_surv(i-1) * mean_surv(i-1);
end

% Step 3: Fit Gompertz model to cumulative survival
% Gompertz survival function: S(a) = exp(-(b1/b2) * (exp(b2*a) - 1) - b3*a)
Model = @(b, x) exp(-(b(1)/b(2)) * (exp(b(2) * x) - 1) - b(3)*x);
b0 = [0.02, 0.05, 0.02];  % Initial guess
b_opt = lsqcurvefit(Model, b0, unique_ages, cum_surv);

% Step 4: Generate smooth survival curve
age_fine = linspace(min(unique_ages), max(unique_ages), 100);
S_fit = Model(b_opt, age_fine);

% Step 5: Estimate mortality rate: mu(a) = -d/dx log(S(a))
dS_dx = gradient(S_fit, age_fine);
mu_fit = -dS_dx ./ S_fit;

% ---- Fit Power Law: mu(a) = c1 * a^c2 ----
power_fun = @(c, x) c(1) * x.^c(2);
c0 = [0.001, 2];
valid_idx = age_fine > 0;
c_opt = lsqcurvefit(power_fun, c0, age_fine(valid_idx), mu_fit(valid_idx));
mu_power = power_fun(c_opt, age_fine);

% ---- Fit Exponential Mortality: mu(a) = d1 * exp(d2 * (a - d3)) ----
exp_fun = @(d, x) d(1) * exp(d(2) * (x - d(3))) + d(4)*exp(d(5)*(x-d(3)));
d0 = [0.001, 0.02, 90,0.002,-0.001];
d_opt = lsqcurvefit(exp_fun, d0, ages_raw, Mort_probs);
mu_exp = exp_fun(d_opt, age_fine);

% === Plotting ===
figure;
subplot(2,1,1)
scatter(unique_ages, cum_surv, 'ro', 'DisplayName', 'Empirical Cumulative'); hold on;
plot(age_fine, S_fit, 'b-', 'LineWidth', 2, 'DisplayName', 'Gompertz Fit');
plot(age_fine, exp(-cumtrapz(age_fine, mu_exp)), 'k--', 'DisplayName', 'Exp. Mortality Integral');
xlabel('Age'); ylabel('Survival Probability');
title('Fitted Survival Curve');
legend('Location', 'best');

subplot(2,1,2)
plot(age_fine, mu_fit, 'k-', 'LineWidth', 2, 'DisplayName', 'Gompertz Mortality'); hold on;
plot(age_fine, mu_power, 'r--', 'LineWidth', 2, 'DisplayName', 'Power Law Fit');
plot(age_fine, mu_exp, 'b--', 'LineWidth', 2, 'DisplayName', 'Exponential Fit');
xlabel('Age'); ylabel('Mortality Rate \mu(a)');
title('Estimated Mortality Function');
legend(); grid on;


%%
unique_ages = unique(ages_raw);
n_ages = length(unique_ages);

% Infer number of draws from the length
n_draws = length(surv_probs_raw) / n_ages;

% Reshape: each row = one draw; each column = one age class
surv_matrix = reshape(surv_probs_raw, [n_ages, n_draws])';  % size = [n_draws x n_ages]

cum_surv_matrix = cumprod(surv_matrix, 2);  % cumulative survival along age axis
cum_surv_mean = mean(cum_surv_matrix, 1);
cum_surv_median = median(cum_surv_matrix, 1);
cum_surv_low = prctile(cum_surv_matrix, 2.5, 1);   % 2.5% quantile
cum_surv_high = prctile(cum_surv_matrix, 97.5, 1); % 97.5% quantile

figure;
hold on;
fill([unique_ages; flipud(unique_ages)], ...
     [cum_surv_low, fliplr(cum_surv_high)], ...
     [0.9 0.9 1], 'EdgeColor', 'none', 'DisplayName', '95% CI');
plot(unique_ages, cum_surv_median, 'b-', 'LineWidth', 2, 'DisplayName', 'Median');
xlabel('Age'); ylabel('Cumulative Survival');
title('Cumulative Survival with Uncertainty');
legend(); grid on;


%%
% === Define arbitrary hazard function: e.g., 2nd-degree polynomial
haz_fun = @(c, a) 1e-2*(c(1)*a.^2 + c(2)*a + c(3));

% === Survival from integrated hazard
survival_model = @(c, a) exp(-cumtrapz(a, haz_fun(c, a)));

% === Initial guess for parameters
c0 = [1e-4, 1e-2, 0.1];

% === Fit to cumulative survival
loss_fun = @(c) survival_model(c, unique_ages) - cum_surv;  % vector residuals
c_fit = lsqnonlin(loss_fun, c0);

% === Fit to mortality directly
c_fit2 = lsqcurvefit(haz_fun, c0, ages_raw, Mort_probs);

% === Evaluate fit
S_custom = survival_model(c_fit, age_fine);
mu_custom = haz_fun(c_fit, age_fine);

% Add to plots
figure;
subplot(2,1,1)
scatter(unique_ages, cum_surv, 'ro'); hold on;
plot(age_fine, S_custom, 'g-', 'LineWidth', 2, 'DisplayName', 'Survival Fit');
plot(age_fine, survival_model(c_fit2,age_fine), 'r--', 'LineWidth', 2, 'DisplayName', 'Hazard Fit');
xlabel('Age'); ylabel('Survival');
legend();

subplot(2,1,2); hold on
scatter(ages_raw,Mort_probs,1,'blue','filled');
plot(age_fine, mu_custom, 'g-', 'LineWidth', 2, 'DisplayName', 'Survival Fit');
plot(age_fine, haz_fun(c_fit2,age_fine),'r--','LineWidth',2,'DisplayName',"Hazard Fit");
xlabel('Age'); ylabel('Hazard \mu(a)');
legend(); grid on;


%% ---- Fit Fertility as Sum of Two Gaussians ----
age_data = DAT.fit_sim(:,1);
fert_data = DAT.fit_sim(:,4);

% Define age bins (same as used for population density or choose here)
age_edges = [10 15 20 25 30 35 40 45 50]; % Example bins
age_bin_centers = 0.5 * (age_edges(1:end-1) + age_edges(2:end));

% Initialize array to hold mean fertility per bin
mean_fert = zeros(size(age_bin_centers));

% Compute mean fertility per age class
for i = 1:length(age_bin_centers)
    in_bin = (age_data >= age_edges(i)) & (age_data < age_edges(i+1));
    mean_fert(i) = mean(fert_data(in_bin));
end

% Two-Gaussian model
fert_model = @(p, a) ...
    p(1) * exp(-0.5 * ((a - p(2)) / p(3)).^2) + ...
    p(4) * exp(-0.5 * ((a - p(5)) / p(6)).^2);

% Initial guess: [A1, mu1, sigma1, A2, mu2, sigma2]
p0 = [0.05, 20, 5, 0.03, 35, 7];

% Lower and upper bounds to constrain fitting
lb = [0, 10, 1, 0, 25, 1];
ub = [1, 40, 20, 1, 50, 20];

p_opt = lsqcurvefit(fert_model, p0, age_data, fert_data, lb, ub);

% One-Gaussian model
fert_model2 = @(p, a) ...
    p(1) * exp(-0.5 * ((a - p(2)) / p(3)).^2);

% Initial guess: [A1, mu1, sigma1, A2, mu2, sigma2]
p0 = [0.05, 20, 5];

% Lower and upper bounds to constrain fitting
lb = [0, 10, 1];
ub = [1, 40, 20];

p_opt2 = lsqcurvefit(fert_model2, p0, age_data, fert_data, lb, ub);


% Generate smooth fit curve
age_fine = linspace(min(age_data), max(age_data), 200);
fert_fit = fert_model(p_opt, age_fine);
fert_fit2 = fert_model2(p_opt2, age_fine);

% Plot fit
figure;
scatter(age_data, fert_data, 10, 'filled', 'DisplayName', 'Fertility Data'); hold on;
plot(age_fine, fert_fit2, 'g-', 'LineWidth', 2, 'DisplayName', 'Gaussian Fit');
plot(age_fine, fert_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Two-Gaussian Fit');
plot(age_bin_centers, mean_fert, 'k-', 'LineWidth', 2, 'DisplayName', 'Mean');
xlabel('Age'); ylabel('Fertility Rate');
title('Fitted Fertility Curve (Two Gaussians)');
legend(); grid on;