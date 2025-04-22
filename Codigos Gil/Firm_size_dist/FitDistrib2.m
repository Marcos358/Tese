% Fit the distribution for n_workers

% choose sector
services = 1;

% load data
if services==1
    load draw_n_workers_services.mat;
    X = double(n_workers_services);
else
    load draw_n_workers_industry.mat;
    X = double(n_workers_industry);
end

N = length(X);
X = X';

% Compute emprirical cdf (ecdf) for scatter plot
[f,x] = ecdf(X);

% grid for theoretical cdf
grid = min(X):1:max(X);

% Theoretical cdf
% 1) By fitdist
gpj = fitdist(X,'gp','theta',0.9999)

Y1 = gpcdf(grid,gpj.k,gpj.sigma,gpj.theta);

% 2) Maximum likelihood (pareto type 1)
xi_hat = N/sum(log(X));

Y2 = gpcdf(grid,xi_hat,xi_hat,1);

% plot
scatter(x(1:20),f(1:20)); hold on
plot(grid(1:20),Y1(1:20)); hold on
plot(grid(1:20),Y2(1:20)); hold on

legend('ecdf','fitdist','maximum likelihood');


