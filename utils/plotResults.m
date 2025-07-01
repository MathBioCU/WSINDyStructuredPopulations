function [growthrate,deathrate,birthrate,u_dd] =plotResults(U,U_n,U_t,t_cell,x,Growth_cell,Death_cell,Birth_cell,Growth_w,Death_w,Birth_w,fluxflag)
%Solves the learned structured population model usings an upwind scheme


t = t_cell{1};
t_data = t_cell{2};


dt = t(2)-t(1);
dx = 1e-2;
dx_data = x(2) - x(1);


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

if fluxflag
    try
    x_top = fzero(@(x)growthrate(x,1,1),x(end));
    if x_top > x(2)
        x_sim = (x(1):dx:x_top)';
    else
        x_sim = (x(1):dx:x(end))';
    end
    u_0 = interp1(x,U(:,1),x_sim,"linear","extrap"); 
    catch
    x_sim = (x(1):dx:x(end))';
    u_0 = interp1(x,U(:,1),x_sim,"linear","extrap");
    end
else
    x_sim = (x(1):dx:x(end))';
    u_0 = interp1(x,U(:,1),x_sim,"linear","extrap");
end
u_0(isnan(u_0)) = 0;

function dudt = RHS(t,u,x,dx,growthrate,birthrate,deathrate,flag)
    dudt = zeros(size(u));
    N = dx*sum(u);
    % N = trapz(a,u);
    dudt(1) = -(1/dx)*(growthrate(x(1),u(1),N) ...
        - dx*sum(birthrate(x,u,N)))...
        +deathrate(x(1),u(1),N);


    if flag
        dudt(2:end-1) = -(1/dx)*(growthrate(x(2:end-1),u(2:end-1),N) ...
            - growthrate(x(1:end-2),u(1:end-2),N))...
            +deathrate(x(2:end-1),u(2:end-1),N);
    
        dudt(end) = -(1/dx)*(0 - growthrate(x(end-1),u(end-1),N))...
            +deathrate(x(end),u(end),N);
    else
        dudt(2:end) = -(1/dx)*(growthrate(x(2:end),u(2:end),N) ...
            - growthrate(x(1:end-1),u(1:end-1),N))...
            +deathrate(x(2:end),u(2:end),N);
    end

end

tic;
[t,u] = ode45(@(t,u) RHS(t,u,x_sim,dx,growthrate,birthrate,deathrate,fluxflag), t,u_0);
toc

u_dd = u';

u_dd = interp1(x_sim,u_dd,x); u_dd(isnan(u_dd)) = 0;

figure()
subplot(1,3,1)
surf(x,t_data,U_n','EdgeColor','none')
ylabel('Time')
xlabel('Structure')
title("Data")
M = max(u,[],"all")*(3/3);
colormap("jet"); clim([0,M]); colorbar();
% zlim([0,M])
ylim([0,t(end)])
view(2)

subplot(1,3,2)
surf(x,t,u_dd','edgecolor','none')
xlabel('Structure')
ylabel('Time')
colormap("jet"); clim([0,M]);%colorbar();
title("Learned model")
% zlim([0,M])
ylim([0,t(end)])
view(2)

% subplot(2,2,3)
% surf(x,t,min(abs(u-U')./abs(u),1),'edgecolor','none')
% view(2)
% colorbar()
% xlabel('Structure')
% ylabel('Time')
% colormap("jet")
% title("Relative error with true dynamics")

subplot(1,3,3)
hold on
plot(t,dx*sum(u,2),'LineWidth',2,'DisplayName','Learned')
plot(t,dx_data*sum(U,1),'--','LineWidth',2,'DisplayName','True Dynamics')
scatter(t_data,U_t, 4,'DisplayName','Data')
xlabel("Time")
title("Total Population")
lg = legend();
lg.Location = 'Best';




end