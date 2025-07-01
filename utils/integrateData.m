function u_total = integrateData(x,t,u_data,type)
dx = mean(diff(x));
if isequal(type,"normal")    
    u_total = dx * sum(u_data,1);
elseif isequal(type,"lognormal")
    u_total = zeros(size(t));
    for i = 1:length(t)
    valid_idx = find(u_data(:,i)~=0);
    log_u = log(u_data(valid_idx,i));
    sp_log = fit(x(valid_idx), log_u, 'smoothingspline');
    u_smooth_log = exp(feval(sp_log, x(valid_idx)));
    u_total(i) = dx* sum(u_smooth_log);
    end
end