function [Xi,its,thrs_EL] = sparsifyDynamics(Theta,dXdt,lambda,n,gamma,M)
    % Copyright 2015, All Rights Reserved
    % Code by Steven L. Brunton
    % For Paper, "Discovering Governing Equations from Data: 
    %        Sparse Identification of Nonlinear Dynamical Systems"
    % by S. L. Brunton, J. L. Proctor, and J. N. Kutz
    %
    % Modified by Daniel A. Messenger, 2020 
    
    [~,nn] =size(Theta);
    if  gamma ~= 0
        Theta = [Theta;gamma*eye(nn)];
        dXdt = [dXdt;zeros(nn,n)];
    end
    
    Xi = Theta \ dXdt;  % initial guess: Least-squares

    % noise_type = 'multiplicative';
    % % Apply weighted least squares to account for multiplicative lognormal noise
    % if any(strcmp(noise_type, 'multiplicative'))
    %     R = dXdt - Theta * Xi;
    %     w = 1 ./ max(sum(R.^2, 2), 1e-6);  % per-row residual variance
    %     W = diag(w);
    %     Theta = W * Theta;
    %     dXdt = W * dXdt;
    %     Xi = Theta \ dXdt;
    % end

    if ~isempty(M)
        Xi = M.*Xi;
        bnds = norm(dXdt)./vecnorm(Theta)'.*M; 
        LBs = lambda*max(1,bnds);ones(size(bnds));
        UBs = 1/lambda*min(1,bnds);ones(size(bnds));
        thrs_EL = [LBs bnds UBs];
    else
        thrs_EL = [];
    end
    
    smallinds = 0*Xi;
    for j=1:nn
        if ~isempty(M)
            smallinds_new = or(abs(Xi)<LBs,abs(Xi)>UBs);
            if all(smallinds_new(:)==smallinds(:))
                its = j;
                return
            else
                smallinds = smallinds_new;
                Xi(smallinds)=0;    
                for ind=1:n
                    Xi(~smallinds,ind) = M(~smallinds).*(Theta(:,~smallinds)\dXdt(:,ind));
                end
            end
        else
            smallinds_new = (abs(Xi)<lambda);
            if all(smallinds_new(:)==smallinds(:))
                its = j;
                return
            else
                smallinds = smallinds_new;
                Xi(smallinds)=0;
                for ind = 1:n        
                    biginds = ~smallinds(:,ind);
                    Xi(biginds,ind) = Theta(:,biginds)\dXdt(:,ind);
                end
            end
        end
    end
    its = nn;
end
