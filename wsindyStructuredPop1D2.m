function [w_pder, w_oder, TransportTrials, SourceTrials, BoundaryTrials, ...
          TransportTags, SourceTags, BoundaryTags, G, b, G_ode2, b_ode2, phix] = ...
    wsindyStructuredPop1D2(x, t, u_data, u_total, libraryparams, ...
                           tf_rad_x, tf_rad_t, tf_power, qp, ...
                           Threshhold_pde, Threshhold_ode, Sparsity_w, var, ...
                           BCV, verbose, want_plots)
% WSINDYSTRUCTUREDPOP1D Weak SINDy for Structured Population Models in 1D
%
% This function applies the weak-form Sparse Identification of Nonlinear 
% Dynamics (WSINDy) to structured population models in one spatial dimension.
%
% INPUTS:
%   x              - Spatial grid (Nx-by-1)
%   t              - Temporal grid (1-by-Nt)
%   u_data         - State data u(x,t) (Nx-by-Nt)
%   u_total        - Total population over time (1-by-Nt)
%   libraryparams  - Cell array of structs for {Transport, Source, Boundary} trial library parameters
%   tf_rad_x       - Support radius (fraction of x domain) for test functions
%   tf_rad_t       - Support radius (fraction of t domain) for test functions
%   tf_power       - Degree of polynomial test functions
%   qp             - Subsampling factor for test function placement
%   Threshhold_pde - Regularization path or threshold for PDE sparsification
%   Threshhold_ode - Same as above, but for ODE part
%   Sparsity_w     - Desired sparsity level for PDE and ODE weights
%   var            - Variance for noise model (used in normalization)
%   BCV            - Boolean: enable ODE-based cross-validation refinement
%   verbose        - Boolean: print learned model details
%   want_plots     - Boolean: show residual/sparsity plots
%
% OUTPUTS:
%   w_pder, w_oder - Identified coefficients for PDE (spatial terms) and ODE (boundary terms)
%   ___Trials      - Cell arrays of trial functions
%   ___Tags        - Cell arrays of corresponding tag strings
%   G, b           - Assembled linear system (for diagnostic use)
%   G_ode2, b_ode2 - ODE system with learned source substituted (for BCV)
%   phix           - Test functions over space (for visualization or reuse)
%
% Author: Rainey Lyons, University of Colorado Boulder
% Contact: rainey.lyons@colorado.edu

% Initialize threshold method (ignored — suggest fixing logic later)
Thresholding = 'MSTLS'; % Ignored 'LB'

% Spatial and temporal discretizations (assumed approximately uniform)
dx = mean(diff(x));
dt = mean(diff(t));

% Unpack trial function library parameters
TransportParams = libraryparams{1};
SourceParams    = libraryparams{2};
BoundaryParams  = libraryparams{3};

%% ------------------------------------------------------------------------
% 1. Generate Test Functions (time & space)
% -------------------------------------------------------------------------

TestFuncs_t = {}; TestDerivs_t = {};
TestFuncs_x = {}; TestDerivs_x = {};

for tf_t = tf_rad_t
    for tf_x = tf_rad_x
        % Define support regions
        try
            Supports_t = [t(1:ceil(end*(1-tf_t)))', t(ceil(end*tf_t):end)'];
        catch
            Supports_t = [t(1:ceil(end*(1-tf_t)))', t(ceil(end*tf_t):end-1)'];
        end
        try
            Supports_x = [x(1:ceil(end*(1-tf_x))), x(ceil(end*tf_x):end)];
        catch
            Supports_x = [x(1:ceil(end*(1-tf_x))), x(ceil(end*tf_x):end-1)];
        end
        
        % Subsample test function locations
        Supports_x = Supports_x(1:qp:end, :);
        Supports_t = Supports_t(1:qp:end, :);

        % Generate polynomial test functions
        [tf_t_set, d_tf_t_set, ~] = generatePolynomialTestFuncs(Supports_t(:,1), Supports_t(:,2), tf_power, tf_power);
        [tf_x_set, d_tf_x_set, ~] = generatePolynomialTestFuncs(Supports_x(:,1), Supports_x(:,2), tf_power, tf_power);

        % Append to full lists
        TestFuncs_t = [TestFuncs_t; tf_t_set];
        TestDerivs_t = [TestDerivs_t; d_tf_t_set];
        TestFuncs_x = [TestFuncs_x; tf_x_set];
        TestDerivs_x = [TestDerivs_x; d_tf_x_set];
    end
end

% Vectorize test functions
phit  = vectorizeTests(TestFuncs_t, t);
Dphit = vectorizeTests(TestDerivs_t, t);
phix  = vectorizeTests(TestFuncs_x, x);
Dphix = vectorizeTests(TestDerivs_x, x);

%% ------------------------------------------------------------------------
% 2. Generate Trial Functions (Transport, Source, Boundary)
% ----------------------------------------------------------------------------

[BoundaryTrials, BoundaryTags] = generateTrialFunctions(BoundaryParams);
[SourceTrials,   SourceTags]   = generateTrialFunctions(SourceParams);
[TransportTrials,TransportTags]= generateTrialFunctions(TransportParams);

% Vectorize the trial functions
BTrials = vectorizeTrials(BoundaryTrials, t, x, u_data, u_total);
STrials = vectorizeTrials(SourceTrials,   t, x, u_data, u_total);
TTrials = vectorizeTrials(TransportTrials,t, x, u_data, u_total);

%% -------------------------------------------------------------------------
% 3. Construct Linear System: Integrate trials/data against test functions
% --------------------------------------------------------------------------

% Integrate against spatial test derivatives 
U  = integrateDataAgainstTests(u_data, phix, 1);
T  = integrateTrialsAgainstTests(TTrials, Dphix, 1);
S  = integrateTrialsAgainstTests(STrials, phix, 1);
TS = cat(3, T, S); % Combined PDE terms
PDE_Tags = [TransportTags, SourceTags];

% ODE-level integrations for total population dynamics
S_ode = integrateTrialsAgainstTests(STrials, {dx * ones(size(x))}, 1);
B_ode = integrateTrialsAgainstTests(BTrials, {dx * ones(size(x))}, 1);

% RHS 
b_pde = integrateDataAgainstTests(U, Dphit, 2);
b_pde = reshape(-b_pde, [], 1);

% LHS
G = integrateTrialsAgainstTests(TS, phit, 2);
G = reshape(G, [], size(TS, 3));
G = paddata(G, [size(G, 1), size(G, 2) + length(BTrials)]);  % pad for boundary terms

% Assemble ODE system
b_ode = integrateDataAgainstTests(u_total, Dphit, 2);
b_ode = reshape(-b_ode, [], 1);

G_ode = [ ...
    integrateTrialsAgainstTests(S_ode, phit, 3), ...
    integrateTrialsAgainstTests(B_ode, phit, 3) ...
];
G_ode = reshape(G_ode, [], length(SourceTrials) + length(BoundaryTrials));
G_ode = paddata(G_ode, [size(G_ode,1), size(TS,3) + length(BTrials)], "Side", "leading");

% Normalize by noise level
G_ode = G_ode ./ exp(var/2);

% Final full system
G = [G; G_ode];
b = [b_pde; b_ode];

%% ------------------------------------------------------------------------
% 4. Sparse Regression (MSTLS or LB)
% ------------------------------------------------------------------------

if isscalar(Threshhold_pde)
    % --- LB-style ensemble sparsification (optional)
    LBp.nE1 = 0.5;    % Fraction of library used per ensemble
    LBp.nE2 = 200;    % Number of ensembles
    LBp.ensT = 0.9;   % Retention threshold for ensemble features
    LBp.nE3 = 200;    % Bootstrap samples for uncertainty
    LBp.ensT2 = 0.5;  % Threshold for UQ phase
    LBp.DB = true;    % Distributional bagging enabled

    w_pde = sparsifyDynamicsLB(G, b, Threshhold_pde, 0, ones(size(G, 2), 1), LBp);
    w_ls  = G \ b;  % Least-squares for comparison

elseif strcmp(Thresholding, 'MSTLS')
    % --- Multi-stage thresholded least squares
    [w_pde, w_ls, ~, tb, loss_vals, fv, sv] = ...
        sequential_thresholding_ls(G, b, Threshhold_pde, Sparsity_w(1), ...
                                   length(b_pde), length(PDE_Tags));

    if verbose && want_plots
        figure;  
        subplot(1,2,1); hold on; 
        title('Loss vs \lambda');
        plot(Threshhold_pde, loss_vals); 
        scatter(tb, loss_vals(Threshhold_pde == tb), 'o'); 
        xlabel('\lambda'); xscale('log');

        subplot(1,2,2); hold on;
        title('Model Diagnostics');
        yyaxis left;  plot(Threshhold_pde, fv, 'DisplayName', 'Residual');
        yyaxis right; plot(Threshhold_pde, sv, 'DisplayName', 'Sparsity');
        xlabel('\lambda'); xscale('log'); legend(); 
    end
end

% --- Partition learned weights
w_T = w_pde(1:length(TransportTrials));
w_S = w_pde(length(TransportTrials)+1 : length(TransportTrials)+length(SourceTrials));
w_B = w_pde(length(TransportTrials)+length(SourceTrials)+1:end);

% --- Display learned PDE components
if verbose
    disp("Learned Transport:");
    disp([w_T(w_T ~= 0), [TransportTags{w_T ~= 0}]']);
    
    disp("Learned Source:");
    disp([w_S(w_S ~= 0), [SourceTags{w_S ~= 0}]']);
    
    disp("Boundary Terms:");
    disp([w_B(w_B ~= 0), [BoundaryTags{w_B ~= 0}]']);
    
    disp("Residual (PDE): " + num2str(norm(G*w_pde - b) / norm(b)));
    disp("Residual (LS):  " + num2str(norm(G*w_ls - b) / norm(b)));
end

% Output initial weight partitions
w_pder = [w_T; w_S];
w_oder = w_B;

% (may be overwritten by BCV)
G_ode2 = 1; b_ode2 = 1;


%% ----------------------------------------------------------------------
% 5. Boundary Correction via Cross-Validation (if enabled)
% --------------------------------------------------------------------

if BCV
    % Construct learned source function 
    Learned_F = zeros(size(STrials{1}));
    for i = 1:length(w_S)
        Learned_F(:,:,i) = w_S(i) * STrials{i};
    end
    Learned_F = sum(Learned_F, 3);  
    Learned_F = dx * sum(Learned_F, 1);  % Integrate over space

    % Build ODE-only system using corrected source term
    B = zeros(length(t), length(BTrials));
    for i = 1:length(BTrials)
        B(:, i) = dx * sum(BTrials{i}, 1);
    end

    % RHS: d/dt(u_total) + source(u)
    b_ode2 = integrateDataAgainstTests(u_total, Dphit, 2) + ...
             integrateDataAgainstTests(Learned_F, phit, 2);
    b_ode2 = reshape(-b_ode2, [], 1);

    G_ode2 = integrateTrialsAgainstTests(B_ode, phit, 3);
    G_ode2 = reshape(G_ode2, [], length(BoundaryTrials));
    G_ode2 = G_ode2 / exp(var/2);  % Normalize for noise

    % Solve again with sparsity-promoting regression
    if isscalar(Threshhold_ode)
        w_ode = sparsifyDynamics(G_ode2, b_ode2, Threshhold_ode, 1, 0, ...
                                 ones(size(G_ode2,2),1));
        w_ls2 = G_ode2 \ b_ode2;
    else
        [w_ode, w_ls2, ~, ~, loss_vals, ~, ~] = ...
            sequential_thresholding_ls(G_ode2, b_ode2, Threshhold_ode, ...
                                       Sparsity_w(2), length(b_pde), length(PDE_Tags));
    end

    % Display ODE component results
    if verbose
        disp("---- ODE Cross-Validation ----");
        disp("Learned Boundary:");
        disp([w_ode(w_ode ~= 0), [BoundaryTags{w_ode ~= 0}]']);
        disp("Residual (ODE):    " + num2str(norm(G_ode2*w_ode - b_ode2)/norm(b_ode2)));
        disp("Residual (ODE LS): " + num2str(norm(G_ode2*w_ls2 - b_ode2)/norm(b_ode2)));
    end

    % Check if new ODE weights are more consistent with full system
    if ~isequal(~w_ode, ~w_B) && ~isempty(w_ode(w_ode ~= 0))
        supp1 = logical([ones(size(w_T)); ones(size(w_S)); w_B ~= 0]);
        supp2 = logical([ones(size(w_T)); ones(size(w_S)); w_ode ~= 0]);

        if sum(w_ode ~= 0) == sum(w_B ~= 0)
            common_idx = or(supp1, supp2);
        else
            common_idx = and(supp1, supp2);
        end

        w_corrected = zeros(size(w_pde));
        w_corrected_temp = sequential_thresholding_ls( ...
            G(:,common_idx), b, Threshhold_pde, Sparsity_w(1), ...
            length(b_pde), length(PDE_Tags));
        w_corrected(common_idx) = w_corrected_temp;

        % Split corrected weights
        w_T2 = w_corrected(1:length(TransportTrials));
        w_S2 = w_corrected(length(TransportTrials)+1 : length(TransportTrials)+length(SourceTrials));
        w_B2 = w_corrected(length(TransportTrials)+length(SourceTrials)+1:end);

        if verbose
            disp('---- ODE-Corrected Model ----');
            disp("Transport:");
            disp([w_T2(w_T2 ~= 0), [TransportTags{w_T2 ~= 0}]']);
            disp("Source:");
            disp([w_S2(w_S2 ~= 0), [SourceTags{w_S2 ~= 0}]']);
            disp("Boundary:");
            disp([w_B2(w_B2 ~= 0), [BoundaryTags{w_B2 ~= 0}]']);
            disp("Residual (corrected): " + num2str(norm(G*w_corrected - b)/norm(b)));
        end

        % Update outputs
        w_pder = [w_T2; w_S2];
        w_oder = w_B2;

    else
        % If original + new ODE combined gives better residual, prefer that
        w_combined = [w_T; w_S; w_ode];
        if norm(G*w_combined - b) < norm(G*w_pde - b)
            if verbose
                disp("Replaced with ODE-corrected model. Residual: " + ...
                     num2str(norm(G*w_combined - b)/norm(b)));
            end
            w_pder = [w_T; w_S];
            w_oder = w_ode;
        end
    end
end

end
