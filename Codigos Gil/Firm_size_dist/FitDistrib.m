clear clc
% Compute distribution that fits the data on firm-size distribution
formal = 1;
industry = 0;
% Load data
load dist_n_workers_formal.mat n_firms_industry n_firms_services;

x_s = log(double(n_firms_industry));
x_i = log(double(n_firms_services));
clear n_firms_industry n_firms_services;

load dist_n_workers_informal.mat n_firms_industry n_firms_services

xi_s = log(double(n_firms_industry));
xi_i = log(double(n_firms_services));

% Location Parameter (Pareto)
Xm = 1;



% Histogram
% formal
if formal ==1
    if industry ==1
        x = x_i; % industry
    else
        x = x_s; % services
    end
else
    if industry ==1
        x = xi_i; % industry
    else
        x = xi_s; % services
    end
end

% Sample size
N = length(x);

% Generate uniform draws
X1 = rand(N,1);
% Shape parameter (maximum likelihood)
xi_hat = N/(sum(x) - N*log(Xm));

sample = GeneratePareto(X1,1,xi_hat,x,formal);


