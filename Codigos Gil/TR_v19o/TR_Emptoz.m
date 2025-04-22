%%% Employment to productivity

function Fz = TR_Emptoz(FixPar,Prices,tauW,tauY_j)

% Prices
w = Prices.w;
p_j = Prices.p_j;

% Parameters
J = length(p_j);
alpha_j = FixPar.alpha_j;
theta_j = FixPar.theta_j;
lambda_jk = FixPar.lambda_jk;
%
omegap_j = log( (alpha_j.*(1-tauY_j)./((1+tauW)*w)).^(1-theta_j) .* theta_j.^theta_j);
Omegap_j = omegap_j + theta_j.*sum(log(lambda_jk.^lambda_jk)');
Lambda0 = theta_j'.*lambda_jk;
Theta = Omegap_j' + (eye(J) - Lambda0)* log(p_j');

% Employment do productivity function
% Fz:(sector j,firm employment) -> firm productivity (formal)
Fz = @(j,emp) exp((1-alpha_j(j)-theta_j(j))*log(emp) - Theta(j));
    
end