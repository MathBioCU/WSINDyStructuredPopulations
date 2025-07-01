function [w_best,w_Ls,coefficients, threshold_best,loss_vals,fit_vals,spars_vals] = sequential_thresholding_ls(A, b, threshold_array,sparsity,row_pde,col_pde)
    % Implements a modified Sequential Thresholding Least-Squares algorithm.
    % A: Regression matrix
    % b: Target vector
    % threshold_array: Array of threshold values to test
    % sparsity: weight on sparsity constraint of loss function
    % Output:
    %   - w_best: learned weights which minimize the loss function
    %   - coefficients: Cell array containing results for different thresholds
    %   - loss_vals: Array of values from los function
    

    num_thresholds = length(threshold_array);
    coefficients = cell(1, num_thresholds);
    loss_vals = inf*ones(1,num_thresholds);
    length_vec = [];
    


    
    if size(A,1)>row_pde
    % row_ode = size(b,1)-row_pde;
    % alpha = sqrt(row_pde/(row_ode*1));
    % A(row_pde+1:end,:) = alpha * A(row_pde+1:end,:);
    % b(row_pde+1:end,:) = alpha * b(row_pde+1:end,:);

    % colnorm = vecnorm(A);
    % A = A./ colnorm;

    w_Ls = A \ b;
    % Perform STLS for each threshold
    for t_idx = 1:num_thresholds
        lambda = threshold_array(t_idx); % Current threshold
        
        w = sparsifyDynamics(A,b,lambda,1,0,ones(size(A,2),1));
        coefficients{t_idx} = w;

        %sparsity = 1-0.5*1/length(w)/estimate_delta_residual(A,w,w_Ls);
        
        fit_vals(t_idx) = norm(A*(w - w_Ls),2)/norm(A*w_Ls,2);norm(A*w - b)/norm(b); 

        % fit_vals(t_idx) = norm(A(1:row_pde,:)*w - b(1:row_pde))/norm(b(1:row_pde))...
            % + norm(A(row_pde+1:end,:)*w - b(row_pde+1:end))/norm(b(row_pde+1:end));
        
        % %  norm(A*(w - w_Ls),2)/norm(A*w_Ls,2);  norm(A*(w - w_Ls),2)/norm(A*w_Ls,2);  
        spars_vals(t_idx) = sum(w~=0)/length(w);

        % spars_vals(t_idx) = sum(w(1:col_pde)~=0)/length(w(1:col_pde)) + sum(w(col_pde+1:end)~=0)/length(w(col_pde+1:end)); log(1+sum(w~=0)/length(w_Ls));
        miss_vals(t_idx) = 100 * (sum(w~=0)<3);

        loss_vals(t_idx) = (1-sparsity) * fit_vals(t_idx) ...
                            +  sparsity * spars_vals(t_idx) ;
        length_vec = [length_vec sum(w~=0)];
        % loss_vals(t_idx) = 2*sparsity*norm(A*(w - w_Ls),2) ...
        %                     + 2*(1-sparsity) * norm(w,1);

        % loss_vals(t_idx) = norm(A*w - b,2)/norm(b) + sparsity * sum(w~=0)/length(w);

    end
    else
            % Perform STLS for each threshold
    w_Ls = A \ b;
    for t_idx = 1:num_thresholds
        lambda = threshold_array(t_idx); % Current threshold
        
        w = sparsifyDynamics(A,b,lambda,1,0,ones(size(A,2),1));
        coefficients{t_idx} = w;

        %sparsity = 1-0.5*1/length(w)/estimate_delta_residual(A,w,w_Ls);
        
        fit_vals(t_idx) =  norm(A*w - b)/norm(b);  

        % fit_vals(t_idx) = norm(A(1:row_pde,:)*w - b(1:row_pde))/norm(b(1:row_pde))...
            % + norm(A(row_pde+1:end,:)*w - b(row_pde+1:end))/norm(b(row_pde+1:end));
        
        % %  norm(A*(w - w_Ls),2)/norm(A*w_Ls,2);  

        spars_vals(t_idx) = sum(w~=0)/length(w); log(1+sum(w~=0)/length(w_Ls));
        miss_vals(t_idx) = 100 * (sum(w~=0)<3);

        loss_vals(t_idx) = 2*(1-sparsity)* fit_vals(t_idx) ...
                            + 2 * sparsity*spars_vals(t_idx) ;
        length_vec = [length_vec sum(w~=0)];
        % loss_vals(t_idx) = 2*sparsity*norm(A*(w - w_Ls),2) ...
        %                     + 2*(1-sparsity) * norm(w,1);

        % loss_vals(t_idx) = norm(A*w - b,2)/norm(b) + sparsity * sum(w~=0)/length(w);

    end
    end


    [~,m_idx] = min(loss_vals);
    w_best = coefficients{m_idx};
    threshold_best = threshold_array(m_idx);
    
end

function delta_res = estimate_delta_residual(G, w, w_ls)
    % Compute the residual vector
    r = G * (w - w_ls);
    norm_r = norm(r, 2);

    % Handle the zero-residual case
    if norm_r == 0
        grad = zeros(size(w));
    else
        grad = (G' * r) / norm_r;  % Gradient of residual norm
    end

    
    delta_res = median(abs(grad));
end


