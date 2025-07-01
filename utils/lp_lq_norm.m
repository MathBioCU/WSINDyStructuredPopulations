function norm_val = lp_lq_norm(A, B, p, q)
    % Computes the Lp-Lq norm of the elementwise product of matrices A and B
    % A, B: NxM matrices
    % p: Lp norm in the N direction (rows)
    % q: Lq norm in the M direction (columns)
    
    % Compute the element-wise product
    D = abs(A .* B);
    
    % Compute the Lp norm along the rows (N direction)
    Lp_norms = sum(D.^p, 1).^(1/p);
    
    % Compute the Lq norm along the columns (M direction)
    norm_val = (sum(Lp_norms.^q)^(1/q));
end