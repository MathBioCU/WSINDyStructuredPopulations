function [dg,dd,db] = learnedFunctionError(x,t,U_noisy,U_total,Growth_cell,Death_cell,Birth_cell,Growth_w,Death_w,Birth_w,true_g,true_d,true_b)

dx = mean(diff(x)); dt = mean(diff(t));

birthrate = @(x,u,N) 0;
if ~isempty(Birth_cell)
for i = 1:length(Birth_cell)
    birthrate =@(x,u,N) birthrate(x,u,N) + Birth_cell{i}(x,u,N) .* Birth_w(i);
end
end
birthrate = @(x,u,N) max(birthrate(x,u,N),0);
% figure()
% fsurf(birthrate,[0 30 0 10])

growthrate = @(x,u,N) 0;
if ~isempty(Growth_cell)
for i = 1:length(Growth_cell)
    growthrate =@(x,u,N) growthrate(x,u,N) + Growth_cell{i}(x,u,N) .*Growth_w(i);
end
end

deathrate = @(x,u,N) 0;
if ~isempty(Death_cell)
for i = 1:length(Death_cell)
    deathrate =@(x,u,N) deathrate(x,u,N) + Death_cell{i}(x,u,N) .* Death_w(i);
end
end

dg = lp_lq_norm(growthrate(x,U_noisy,U_total)-true_g(x,U_noisy,U_total).*U_noisy,ones(size(U_noisy)),2,2)./lp_lq_norm(true_g(x,U_noisy,U_total).*U_noisy,ones(size(U_noisy)),2,2);
dd = lp_lq_norm(deathrate(x,U_noisy,U_total)-true_d(x,U_noisy,U_total).*U_noisy,ones(size(U_noisy)),2,2)./lp_lq_norm(true_d(x,U_noisy,U_total).*U_noisy,ones(size(U_noisy)),2,2);
db = lp_lq_norm(birthrate(x,U_noisy,U_total)-true_b(x,U_noisy,U_total).*U_noisy,ones(size(U_noisy)),2,2)./lp_lq_norm(true_b(x,U_noisy,U_total).*U_noisy,ones(size(U_noisy)),2,2);

end