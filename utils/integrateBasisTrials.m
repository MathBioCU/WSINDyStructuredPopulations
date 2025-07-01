function RetMat = integrateBasisTrials(x, t, u_data, u_total, test_t, test_x, x_grid, N_grid, basis_type)
if nargin < 9
    basis_type = 'hat'; % default
end

dx = mean(diff(x_grid));
dN = mean(diff(N_grid));
nx = length(x);
nt = length(t);
nbx = length(x_grid);
nbN = length(N_grid);
n_basis = nbx * nbN;
n_test = length(test_t) * length(test_x);

BasisTensor = zeros(nx, nt, n_basis);
count = 1;
for i = 1:nbx
    phi_x = basisFunctionIntegration(x, x_grid, i, dx, basis_type); % [nx × 1]
    for j = 1:nbN
        phi_N = basisFunctionIntegration(u_total.', N_grid, j, dN, basis_type); % [nt × 1], u_total is [1 × nt] so transpose
        phi_N = phi_N.'; % back to [1 × nt] for broadcasting
        BasisTensor(:, :, count) = phi_x * phi_N; % [nx × 1] × [1 × nt] = [nx × nt]
        count = count + 1;
    end
end

RetMat = zeros(n_test, n_basis);
row = 1;
for ix = 1:length(test_x)
    phi_tx = test_x{ix}; % [nx × 1]
    for it = 1:length(test_t)
        phi_tt = test_t{it}; % [nt × 1]
        test_func = phi_tx * phi_tt; % [nx × nt] test function
        for b = 1:n_basis
            integrand = test_func .* BasisTensor(:, :, b) .* u_data; % [nx × nt]
            RetMat(row, b) = sum(integrand(:)); % total integral
        end
        row = row + 1;
    end
end
end

function phi = basisFunctionIntegration(x, grid, idx, d, basis_type)
switch basis_type
    case 'hat'
        center = grid(idx);
        phi = max(1 - abs((x - center) / d), 0);
    case 'spline'
        if iscolumn(grid)
            knots = [grid(1)-d; grid; grid(end)+d]; % clamped knots, column
            knots = knots';
        else
            knots = [grid(1)-d, grid, grid(end)+d];
        end
        order = 4; % linear B-spline
        coefs = zeros(1, length(knots) - order);
        coefs(idx) = 1;
        pp = spmak(knots, coefs); % spmak expects row knots
        phi = fnval(pp, x);
        phi(phi < 0) = 0;
    case 'poly'
        center = grid(idx);
        phi = (x - center).^2;
    otherwise
        error('Unknown basis type.');
end
phi = phi / (sum(phi) * d); % normalize
end
