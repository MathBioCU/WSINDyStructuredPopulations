%% Load Data
DAT = load("Jackson_etal_data.mat");

%% ---- USER SETTINGS ----
% Define age classes (edges, variable width)
age_edges = 0:1:72; % Cell edges
% Time range filter (optional)
year_min = 1980;
year_max = 2010;

% Age filter (optional)
age_min = min(age_edges);
age_max = max(age_edges);

%% ---- Filter Data ----
ages = DAT.ts_dat(:,2);
years = DAT.ts_dat(:,3);

valid = (years >= year_min) & (years <= year_max) & ...
        (ages >= age_min) & (ages <= age_max);
ages = ages(valid);
years = years(valid);

% Unique time points and normalized time vector (t = 0 at start)
unique_years = sort(unique(years));
t = unique_years - min(unique_years); % time midpoints (in years since start)
n_time = length(t);

% Finite volume cell midpoints for age
x = 0.5 * (age_edges(1:end-1) + age_edges(2:end)); % age midpoints
x = x';
bin_widths = diff(age_edges);
n_bins = length(x);

%% ---- Compute Population Density u(t,x) ----
u = zeros(n_time, n_bins); % rows: time, cols: age bins

for i = 1:n_time
    yr = unique_years(i);
    ages_this_year = ages(years == yr);
    counts = histcounts(ages_this_year, age_edges);
    u(i, :) = counts ./ bin_widths; % density = count / bin width
end
u = u';
%% ---- Save in FV format ----
output_filename = 'fv_population_density.mat';
save(output_filename, 'x', 't', 'u');

disp(['Saved FV-style population density data to ', output_filename]);

%% ---- Optional Plot ----
figure;
imagesc(x, t, u'); axis xy;
ylabel('Time (years since start)');
xlabel('Age (midpoints)');
title('Population Density');
colorbar;
