%% ========================================================================
%  WSINDy for Structured Population Models
% ========================================================================
clear; clc;

% Feel free to adjust parameters in sections 1 and 2 to reproduce figures
% from the paper.

%% ------------------------------------------------------------------------
% 1. Data and Library Selection
% ------------------------------------------------------------------------
DataChoice = 6;   % <-- Choose dataset here (1 to 6)

DataBase = {
    "./ExampleData/AgeStructure_1_03_04.mat",               % L.1 / A.1
    "./ExampleData/AgeStructure_1_exp08_gauss10-5.mat",     % L.2 / A.2
    "./ExampleData/SizeStructure_vb_x_sig2-1.mat",          % L.3 / S.1
    "./ExampleData/SizeStructure_xvb_exp_x.mat",            % L.4 / S.2 
    "./ExampleData/AgeStructure_1_constL_quadR",            % NL.1/ A.3
    "./ExampleData/SizeStructure_vbR_linP_sig2-1L.mat"      % NL.2/ S.3
};

Libraries = {...
    "./ExampleData/Library_A1.m", ...
    "./ExampleData/Library_A2.m",...
    "./ExampleData/Library_S1.m",...
    "./ExampleData/Library_S2.m", ...
    "./ExampleData/Library_NL1.m",...
    "./ExampleData/Library_S3.m"...
};


run(Libraries{DataChoice});
addpath(genpath('./utils'));

%% ------------------------------------------------------------------------
% 2. Experiment Settings
% ------------------------------------------------------------------------
rng(5);  % For 'Typical Examples' figures 

% Data trimming, subsampling, and noise
TrainingCutoff = 0.5;
dt_skip = 10;
dx_skip = 1;
Noise_Ratio = 0.6;       % Standard deviation for lognormal noise
Noise_type = "lognormal";
BoundaryCV = true;      % Not recommended for nonlinear examples
SF = 0.2; % Smoothing factor for variance estimate (must be tuned for each example)

% WSINDy test function parameters
TestFunc_supp_ratio_t = 0.5;
TestFunc_supp_ratio_x = 0.4;
TestFunc_power = 14;
qp_sub = 3; % subsample test functions for speed

% Threshold ranges & Sparsity weight
Threshhold_pde = logspace(-4, 0, 1e3);
Threshhold_ode = logspace(-4, 0, 1e3);
Sparsity_weight = [0.5, 0.5];

% Display & plotting
verbose = 1;
Want_Plots = 0;

% For size-structured simulations
fluxflag = contains(DataBase{DataChoice}, "Size");

%% ------------------------------------------------------------------------
% 3. Load and Preprocess Data
% ------------------------------------------------------------------------
load(DataBase{DataChoice}, "U_exact", "t", "x", "True_funcs", "True_w");

% Check column orientation
if ~iscolumn(x)
    x = x';
end
if iscolumn(t)
    t = t';
end

% Extract ground-truth weights
True_trans = True_w{1};
True_source = True_w{2};
True_pde = [True_trans; True_source];
True_boundary = True_w{3};

% Subsample time and structure
dx1 = x(2) - x(1);
t = t(1:dt_skip:end);
x = x(1:dx_skip:end);
dx = mean(diff(x));
U_exact = U_exact(1:dx_skip:end, 1:dt_skip:end);

% Downsample spatially with averaging if dx_skip > 1
if dx_skip > 1
    U_temp = zeros(length(x), length(t));
    for i = 1:length(x)-1
        U_temp(i,:) = dx1/dx * sum(U_exact((i-1)*dx_skip+1:i*dx_skip, :), 1);
    end
    U_temp(end,:) = dx1/dx * sum(U_exact((length(x)-1)*dx_skip+1:end, :), 1);
    U_exact = U_temp;
end

%% ------------------------------------------------------------------------
% 4. Add Noise and Normalize Population
% --------------------------------------------------------------------------
t_full = t; U_full = U_exact;
t = t(1:floor(TrainingCutoff * end));
U_exact = U_exact(:, 1:length(t));

[U_noisy, Noise_var] = addNoise(U_exact, Noise_Ratio, Noise_type);

U_total = dx * sum(U_noisy, 1); var_est = Noise_Ratio^2;

% % estimating noise (comment out if using exact variance)
% [var_est,~] = estimateVariance(U_noisy,x,t,Noise_type,SF);



U_total = 1/exp(var_est/2) *  U_total;


%% ------------------------------------------------------------------------
% 5. Run WSINDy
% ------------------------------------------------------------------------
tic;
[w_pde, w_ode, TransportTrials, SourceTrials, BoundaryTrials, ...
 Ttags, Stags, Btags, G, b, G_ode, b_ode, phix] = ...
    wsindyStructuredPop1D2(x, t, U_noisy, U_total, ...
        {TransportParams, SourceParams, BoundaryParams}, ...
        TestFunc_supp_ratio_x, TestFunc_supp_ratio_t, ...
        TestFunc_power, qp_sub, ...
        Threshhold_pde, Threshhold_ode, ...
        Sparsity_weight, var_est, BoundaryCV, verbose, 1);
times = toc;

w = [w_pde; w_ode];

%% -----------------------------------------------------------------------
% 6. Evaluate Performance
% ------------------------------------------------------------------------
[E2_pde, Einfty_pde, TPR_pde, true_w] = computePerformance2( ...
    w_pde, w_ode, TransportTrials, SourceTrials, BoundaryTrials, ...
    True_pde, True_boundary, Ttags, Stags, Btags, TrueTags);



%% -----------------------------------------------------------------------
% 7. Display Results
% ---------------------------------------------------------------------------
if verbose
    disp("===== Stats =====");
    disp("Residual (full): " + num2str(norm(G*true_w - b) / norm(b)));
    disp("Residual (ODE):  " + num2str(norm(G_ode * true_w(end-length(w_ode)+1:end) - b_ode) / norm(b_ode)));
    disp("Classes: " + num2str(length(x)));
    disp("Noise type: " + Noise_type);
    disp("Noise ratio: " + Noise_var);
    disp("[σ², estimated]: " + num2str([Noise_Ratio^2, var_est]));
    disp("Compute time: " + num2str(times) + " sec");

    disp("===== Performance =====");
    disp("size(G): " + mat2str(size(G)));
    disp("rank(G): " + num2str(rank(G)));
    disp("cond(G): " + sprintf('%10e', cond(G)));
    disp("E2: " + num2str(E2_pde));
    disp("E∞: " + num2str(Einfty_pde));
    disp("TPR: " + num2str(TPR_pde));
end



%% -----------------------------------------------------------------------
% 8. Plot Results 
% ------------------------------------------------------------------------
% Skip poor fits
if E2_pde >= 3
    return
end

if Want_Plots
    disp("Plotting...");

    [g_learned, f_learned, b_learned, u_dd] = ...
        plotResults(U_full, U_noisy, U_total, {t_full, t}, x, ...
                    TransportTrials, SourceTrials, BoundaryTrials, ...
                    w_pde(1:length(TransportTrials)), ...
                    w_pde(length(TransportTrials)+1:end), ...
                    w_ode, fluxflag);

    if TrainingCutoff < 1
        t_idx = floor(TrainingCutoff * size(U_full, 2)) + 1;
        pred_err = lp_lq_norm(U_full(:,t_idx:end) - u_dd(:,t_idx:end), ...
                              ones(size(u_dd(:,t_idx:end))), 2, 2) / ...
                   lp_lq_norm(U_full(:,t_idx:end), ...
                              ones(size(u_dd(:,t_idx:end))), 2, 2);
        disp("Prediction Error: " + pred_err);
    end

    % Relative errors for each term
    N = U_total;
    disp("==== Relative Weighted Errors ====");
    disp("Transport: " + lp_lq_norm(g_learned(x,U_noisy,N)-True_funcs{1}(x,U_noisy,N).*U_noisy, ...
                                 ones(size(U_noisy)), 1, 1) / ...
                     lp_lq_norm(True_funcs{1}(x,U_noisy,N).*U_noisy, ...
                                 ones(size(U_noisy)), 1, 1));
    disp("Source:    " + lp_lq_norm(f_learned(x,U_noisy,N)-True_funcs{2}(x,U_noisy,N).*U_noisy, ...
                                 ones(size(U_noisy)), 1, 1) / ...
                     lp_lq_norm(True_funcs{2}(x,U_noisy,N).*U_noisy, ...
                                 ones(size(U_noisy)), 1, 1));
    disp("Boundary:  " + lp_lq_norm(b_learned(x,U_noisy,N)-True_funcs{3}(x,U_noisy,N).*U_noisy, ...
                                 ones(size(U_noisy)), 1, 1) / ...
                     lp_lq_norm(True_funcs{3}(x,U_noisy,N).*U_noisy, ...
                                 ones(size(U_noisy)), 1, 1));
end
