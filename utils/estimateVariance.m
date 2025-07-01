function [sig_sq_est,residuals,log_u,log_u_smooth] = estimateVariance(u_data,x,t,noise_type,SF)

if isequal(noise_type,'lognormal')
    threshold = 0.01;
    residuals = [];
    for k = 1:1:size(u_data,2)
        valid_idx = u_data(:,k)>threshold;
        log_u = log(u_data(valid_idx,k));
        log_u_smooth = smoothdata(log_u, 1, 'loess','SmoothingFactor', SF);
        residuals = [residuals;log_u-log_u_smooth];
    end


    sig_sq_est = var(residuals(:));
end