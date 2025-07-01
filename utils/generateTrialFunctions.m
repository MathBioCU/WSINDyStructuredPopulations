function [trialFunctions,tags] = generateTrialFunctions(params)
% Generates structured trial functions from given params.
% Outputs:
%   trialFunctions - Cell array of function handles for each birth rate
%       function constructed from every combination of parameters.
%   tags - Cell containing string labels for generated functions


% Initialize empty cell array to store function handles
trialFunctions = {};
tags = {};

%% ---------------------GAUSSIAN-------------------------------------------
if isfield(params,"gaussian")
if ~isnan(params.gaussian.x_peak)
    % Extract Gaussian parameters
    a_peak_values = params.gaussian.x_peak;
    sigma_values = params.gaussian.sigma;
    
    % Generate all combinations of Gaussian parameters
    for a_peak = a_peak_values
        for sigma = sigma_values
            % Define Gaussian distribution for current parameter combination
            trialFunctions{end+1} = @(a,u,N)  normpdf(a,a_peak,sigma) .*u;
            tags{end+1} = "Gauss("+num2str(a_peak)+","+num2str(sigma)+")*u";
        end
    end
end
if isfield(params.gaussian,'scale')
if ~isnan(params.gaussian.scale)
    for i = 1:length(trialFunctions)
        trialFunctions{i} = @(a,u,N)trialFunctions{i}(a,u,N).*params.gaussian.scale;
    end
end
end
end


%% ---------------------Exponential-------------------------------------------
if isfield(params,"exponential")
if ~isnan(params.exponential.k)  
    if isempty(params.exponential.x0)
        % Generate all combinations of Gaussian parameters
        for k = params.exponential.k
            trialFunctions{end+1} = @(a,u,N)  exp(k.*a) .*u;
            tags{end+1} = "exp("+num2str(k)+"x)*u";
        end
    else
        % Generate all combinations of Gaussian parameters
        for k = params.exponential.k
            for x0 = params.exponential.x0
                trialFunctions{end+1} = @(a,u,N)  exp(k.*(a-x0)) .*u;
                tags{end+1} = "exp("+num2str(k)+"(x-"+num2str(x0)+"))*u";
            end
        end
    end
end
end

%% ------------------------GAMMA-------------------------------------------
if isfield(params,"gamma")
if ~isnan(params.gamma.shape)
    % Exctract Gamma parameters
    a_values = params.gamma.shape;
    b_values = params.gamma.scale;
    
    % Generate all combinations of gamma parameters
    for a = a_values
        for b = b_values
            % Define gamma distribution for current parameter combination
            trialFunctions{end+1} = @(x,u,N) gampdf(x,a,b).*u;
            tags{end+1} = "gampdf(a,"+num2str(a)+","+num2str(b)+")*u";
        end
    end
end
end

%% ------------------------SIGMOID-----------------------------------------
if isfield(params,"sigmoid")
if ~isnan(params.sigmoid.x_mid)
    % Extract Sigmoid parameters
    a_maturity_values = params.sigmoid.x_mid;
    k_values = params.sigmoid.k;
    
    % Generate all combinations of Sigmoid parameters
    for a_maturity = a_maturity_values
        for k = k_values
            % Define Sigmoid function for current parameter combination
            trialFunctions{end+1} = @(a,u,N) 1 ./ (1 + exp(-k * (a - a_maturity))).*u;
            tags{end+1} = "sigmoid("+num2str(a_maturity)+","+num2str(k)+")*u"; 
        end
    end
end
end

%% ------------------------POLYNOMIALS--------------------------------------
if isfield(params,"polynomial")
if ~isnan(params.polynomial.p)
    if isfield(params.polynomial,"scale")
        if ~isnan(params.polynomial.scale)
            scale = params.polynomial.scale;
        else
            scale = 1;
        end
    else 
        scale = 1;
    end
    % Extract Polynomial parameters
    p_values = params.polynomial.p;
    for p = p_values
        trialFunctions{end+1} = @(a,u,N) scale.* a.^p .*u;
        if scale == 1
            tags{end+1} = "a^"+num2str(p)+"*u";
        else
            tags{end+1} = "a^"+num2str(p)+"*u"+"*"+num2str(scale);
        end
    end
end
end
% if ~isnan(params.Hermite.p)
%     p_vals = params.Hermite.p;
%     for p = p_vals
%         trialFunctions{end+1} = @(a,N) exp(-(a).^2/2) .* hermiteH(p,a);
%         tags{end+1} = "Hermite "+num2str(p);
%     end
% end

%--------------------Nonlinearities--------------------------------------
if isfield(params,"NL")
if ~isnan(params.NL.Logistic.k)
    L = length(trialFunctions);
    for k = params.NL.Logistic.k
    for i = 1:L
        trialFunctions{end+1} = @(a,u,N) trialFunctions{i}(a,u,N) .* N.^k;
        tags{end+1} = "N^"+num2str(k)+"*"+ tags{i}; 
    end
    end
end

if ~isnan(params.NL.Ricker.k)
    L = length(trialFunctions);
    for i = 1:L
        for k =params.NL.Ricker.k
            trialFunctions{end+1} = @(a,u,N) trialFunctions{i}(a,u,N) .* exp(-N/k);
            tags{end+1} = "exp(-N/"+num2str(k)+")*"+tags{i};
        end
    end
end

if ~isnan(params.NL.Sigmoid.k)
    L = length(trialFunctions);
    for i = 1:L
        for k = params.NL.Sigmoid.k
            for N_m = params.NL.Sigmoid.N_mid
                trialFunctions{end+1}= @(a,u,N) trialFunctions{i}(a,u,N) .* 1./(1+exp(k*(N-N_m)));
                tags{end+1} = "SigNL(k = "+num2str(k)+", N_mid = "+num2str(N_m)+")*"+tags{i};
            end
        end
    end
end
end



%=====Density dependence=================
if isfield(params,"DD")
if ~isnan(params.DD.polynomial.p)
    L = length(trialFunctions);
    for k = params.DD.polynomial.p -1
        if k ==0
            continue
        end
    for i = 1:L
        trialFunctions{end+1} = @(a,u,N) trialFunctions{i}(a,u,N) .* u.^k;
        tags{end+1} =  tags{i}+"^"+num2str(k+1); 
    end
    end
end

end