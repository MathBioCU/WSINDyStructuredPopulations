%% =========================================================================
%  WSINDy for Real Population Data: Jackson et al.
% This script uses public domain data (shared in accordance with its license).
% =========================================================================
clear; clc;

%% -------------------------------------------------------------------------
% 1. Data and Library Selection
% -------------------------------------------------------------------------
DataBase = {"./RealData/fv_population_density.mat"};
Libraries = {"./RealData/Library_Jackson_etal.m"};

DataChoice = 1;

addpath(genpath('./utils'));

%% -------------------------------------------------------------------------
% 2. Experiment Settings
% -------------------------------------------------------------------------

% Subsampling and cutoff
dt_skip = 1;
dx_skip = 1;
TrainingCutoff = 1.0;

% Boundary cross-validation
BoundaryCV = true;  

% Noise estimation 
SF = 0.5;            % Smoothing factor for estimating variance

% WSINDy test function parameters
TestFunc_supp_ratio_t = 0.2;
TestFunc_supp_ratio_x = 0.2;
TestFunc_power = 14;
qp_sub = 1;

% Threshold ranges and sparsity weights
Threshhold_pde = logspace(-4, 0, 1e4);
Threshhold_ode = logspace(-4, 0, 1e4);
Sparsity_weight = [0.3, 0.4]; % can be adjusted depending on library size

% Display and plotting options
verbose = 1;
Want_Plots = 1;

% Domain flux boundary
fluxflag = false;  % 0 flux at right endpoint?

warning('off', 'MATLAB:rankDeficientMatrix');

%% -------------------------------------------------------------------------
% 3. Load and Preprocess Data
% -------------------------------------------------------------------------
load(DataBase{DataChoice}, "u", "t", "x");
U_exact = u;

% Load library functions
run(Libraries{DataChoice});

% Ensure proper vector orientation
if ~iscolumn(x); x = x'; end
if iscolumn(t); t = t'; end

% Compute original dx and downsample
dx1 = x(2) - x(1);
t = t(1:dt_skip:end);
x = x(1:dx_skip:end);
dx = mean(diff(x));

% Downsample U_exact in time
U_exact = U_exact(:,1:dt_skip:end);

% Downsample spatially by averaging if needed
if dx_skip == 1
    U_exact = U_exact(1:end, :);
else
    U_temp = zeros(length(x), length(t));
    for i = 1:length(x)-1
        U_temp(i,:) = dx1/dx * sum(U_exact((i-1)*dx_skip+1:i*dx_skip, :), 1);
    end
    U_temp(end,:) = dx1/dx * sum(U_exact((length(x)-1)*dx_skip+1:end, :), 1);
    U_exact = U_temp;
end

% Normalize population
U_exact = U_exact / (dx * sum(U_exact(:,1)));

%% -------------------------------------------------------------------------
% 4. Trim Data and Estimate Variance
% -------------------------------------------------------------------------
t_full = t;
U_full = U_exact;

% Truncate training data
t = t(1:floor(TrainingCutoff * end));
U_exact = U_exact(:, 1:length(t));
U_noisy = U_exact;

% Estimate observational variance from lognormal model
U_total = dx * sum(U_noisy, 1);
[var_est, ~] = estimateVariance(U_noisy, x, t, "lognormal", SF);
U_total = (1 / exp(var_est/2)) * U_total;

%% -------------------------------------------------------------------------
% 5. Run WSINDy
% -------------------------------------------------------------------------
tic;
[w_pde, w_ode, TransportTrials, SourceTrials, BoundaryTrials, ...
 Ttags, Stags, Btags, G, b, G_ode, b_ode, phix] = ...
    wsindyStructuredPop1D2(x, t, U_noisy, U_total, ...
        {TransportParams, SourceParams, BoundaryParams}, ...
        TestFunc_supp_ratio_x, TestFunc_supp_ratio_t, ...
        TestFunc_power, qp_sub, ...
        Threshhold_pde, Threshhold_ode, ...
        Sparsity_weight, var_est, BoundaryCV, verbose, 0);
times = toc;

w = [w_pde; w_ode];

%% -------------------------------------------------------------------------
% 6. Display Summary Statistics
% -------------------------------------------------------------------------
if verbose
    disp("===== Stats =====");
    disp("Classes: " + num2str(length(x)));
    disp("Estimated variance [σ²]: " + num2str(var_est));
    disp("Computation time: " + num2str(times) + " sec");

    disp("===== Performance Metrics =====");
    disp("size(G): " + mat2str(size(G)));
    disp("rank(G): " + num2str(rank(G)));
    disp("cond(G): " + sprintf('%10e', cond(G)));
end

%% -------------------------------------------------------------------------
% 7. Plot Results and Compare to Data
% -------------------------------------------------------------------------
if Want_Plots
    disp("Plotting...");

    % Ensure scale consistency
    w_pde(length(TransportTrials)+1:end) = w_pde(length(TransportTrials)+1:end);
    w_ode = w_ode * 1;

    % Plot model-predicted transport, source, boundary terms
    % In this plot, `true' dynamics refers to the first entry in the
    % plotResults function and was removed from the paper.
    [g_learned, f_learned, b_learned, u_dd] = ...
        plotResults(U_full, U_noisy, U_total, {t_full, t}, x, ...
                    TransportTrials, SourceTrials, BoundaryTrials, ...
                    w_pde(1:length(TransportTrials)), ...
                    w_pde(length(TransportTrials)+1:end), ...
                    w_ode, fluxflag);

    % Load empirical demographic data
    load("./RealData/Jackson_etal_ages.mat");
    load("./RealData/Jackson_etal_fert.mat");
    load("./RealData/Jackson_etal_surv.mat");

    % ---------------------------------------------------------------------
    % Plot 1: Fertility Function
    figure; subplot(1,2,1); hold on;
    scatter(ages, fert.*surv, 1, 'blue', 'filled');
    plot(x, b_learned(x,1,1), 'k', 'LineWidth', 2, 'DisplayName', 'Learned');
    title("Estimated Fertility");

    % ---------------------------------------------------------------------
    % Plot 2: Mortality Function
    subplot(1,2,2); hold on;
    scatter(ages, 1 - surv, 1, 'blue', 'filled');
    plot(x, -f_learned(x,1,1), 'k', 'LineWidth', 2, 'DisplayName', 'Learned');
    title("Estimated Mortality");

    % ---------------------------------------------------------------------
    % Plot 3: Survival per Age Interval
    figure; hold on;
    scatter(ages, surv, 1, 'blue', 'filled');
    S_learned = zeros(length(x)-1, 1);
    for i = 1:length(x)-1
        S_learned(i) = exp(integral(@(a) f_learned(a,1,1), x(i), x(i+1)));
    end
    plot(x(1:end-1), S_learned, 'k', 'LineWidth', 2, 'DisplayName', 'Learned');
    title("Estimated Survival");


end
