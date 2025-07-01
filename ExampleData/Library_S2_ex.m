
%% Set Boundary Trial function params
% 1. Polynomial: b(x,u) = x^p *u
BoundaryParams.polynomial.p = [];
% 2. Sigmoidal:b(x,u) = 1/ (1+ exp(-k(x-x_maturity)) *u
BoundaryParams.sigmoid.x_mid =[1];
BoundaryParams.sigmoid.k = [2];
% 3. Gaussian: b(x,u) = exp( -(x-x_peak)^2/(2sigma^2)) *u
BoundaryParams.gaussian.x_peak = [];
BoundaryParams.gaussian.sigma = [];
% 4. Gamma: b(x,u) = 1/(b^a Γ(a)) x^(a-1) exp(-x/b)*u. shape = a, scale = b
BoundaryParams.gamma.shape = [];
BoundaryParams.gamma.scale = [];

% Nonlinear terms (Multiplies above functions by these nonlinearities)
% 5. Logistic: b(x,u,N) = b(x,u) * (N^k) (for Logistic include k=1)
BoundaryParams.NL.Logistic.k = [1]; 
% 6. Ricker: b(x,u,N) = b(x,u) * exp(-N/k)
BoundaryParams.NL.Ricker.k = [];
% 7. SigmoidNL: b(x,u,N) = b(x,u) * 1/(1+ exp(k(N-N_mid))) 
BoundaryParams.NL.Sigmoid.k = [];
BoundaryParams.NL.Sigmoid.N_mid = [];%0.5;

%% Set Source Trial function params
% Reaction Terms (f(x,N)u)
% 1. Exponential: f(x,u) = exp(kx)*u
SourceParams.exponential.k =[]; 0.08;% [0.05];0.08:0.02:0.16;
% 2. Sigmoidal:f(x,u) = 1/ (1+ exp(-k(x-x_mid))*u
SourceParams.sigmoid.x_mid = [];
SourceParams.sigmoid.k = [];
% 3. Polynomial: f(x,u) = x^p *u
SourceParams.polynomial.p = [1];

% Nonlinear terms (Multiplies above functions by these nonlinearities)
% 5. Polynomial: f(x,u,N) = f(x,u) * (N^k) (for Logistic include k=1)
SourceParams.NL.Logistic.k = [1];
% 6. Ricker: f(x,u,N) = f(x,u) * exp(-N/k)
SourceParams.NL.Ricker.k = [];
% 7. SigmoidNL: f(x,u,N) = f(x,u) * 1/(1+ exp(k(N-N_mid))) 
SourceParams.NL.Sigmoid.k = [];
SourceParams.NL.Sigmoid.N_mid = [];


%% Set Transport Trail function params
% 1. Exponential: g(x,u) = exp(kx)*u
TransportParams.exponential.k = [];
% 2. Polynomial: g(x,u) = x^p*u (for Von Bertalanffy growth, include p = [0,1])
TransportParams.polynomial.p = [0:1];

% Nonlinear terms
% 3. Polynomial: g(x,u,N) = g(x,u) * (N^k) (for Logistic include k=1)
TransportParams.NL.Logistic.k = [];[1];
% 4. Ricker: g(x,u,N) = g(x,u) * exp(-N/k)
TransportParams.NL.Ricker.k = [4];
% 5. SigmoidNL: g(x,u,N) = g(x,u) * 1/(1+ exp(k(N-N_mid)) 
TransportParams.NL.Sigmoid.k = [];10;
TransportParams.NL.Sigmoid.N_mid =[]; 0.5;
