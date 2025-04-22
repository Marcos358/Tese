%//////////////////////////////////////////////////////////////////////////
% Indirect Tax Reform
%//////////////////////////////////////////////////////////////////////////
% TR_SetParam
% Fixed parameters for all simulations 
%//////////////////////////////////////////////////////////////////////////

% Load data and some parameters
load(FN.modelinputs)

% Structures 
Prices.name = "Equilibrium prices";
FixPar.name = "Fixed parameters";
EndoPar.name = "Endogenously set parameters";
EstPar.name = "Estimated parameters";

%% Set some fixed parameters
%%%
LmaxInf = 5;                % Maximum number of employees in informal firms
jE = 2;                     % Sector 'jE' = numeraire good
%
TaxW = 0.375;               % Standard payroll tax rate [Ulyssea 2018]
PolicyPar.tauW = TaxW;                              % Payroll tax rate [Ulyssea 2018]
PolicyPar.tauWSimp = 0.175;                         % SIMPLES payroll tax rate (Alvarez etal, 2023)
%PolicyPar.RmaxSimp_j = 114*w*[1,1,1];               % Maximum revenue allowed for SIMPLES tax regime [discussion in PolicyPar]
%PolicyPar.RmaxSimp_j = 30*w*[1,1,1];               % Maximum revenue allowed for SIMPLES tax regime [discussion in PolicyPar]
%%%
PolicyPar.RmaxSimp_j = RmaxSimp_j;               % Maximum revenue allowed for SIMPLES tax regime [discussion in PolicyPar]

%% Tax rates
if SO.FirmOptions > 2 % Add SIMPLES
    %%%
    Ts_Th_j = [0,0.018,0.089];%GIL, automatizar         % SIMPLES taxes / total taxes 
    %%%
    VAsh_j = Data.VAj_GDP;                              % GDP shares
    shv2 = VAsh_j(2)/(VAsh_j(2)+VAsh_j(3));
    tauS_tauY3 = (0.65 - shv2*tauS_tauY2)/(1-shv2);
    tauS_tauY_j = [1,tauS_tauY2,tauS_tauY3];
    VAs_VAh_j = Ts_Th_j ./ tauS_tauY_j;
    VAs_VAF_j = 1./(1+1./VAs_VAh_j);
    VAh_VAF_j = 1 - VAs_VAF_j;
    % tY_j = (1./(1-Data.InfVA_j)).* Data.TaxYVA_j;     
    Th_VA_j = Data.TaxYVA_j./(1+Ts_Th_j);
    
    %TaxY_j = Th_VA_j./(VAh_VAF_j.*(1-Data.InfVA_j));    % Standard tax rate
    %TaxY_j(1) = Data.TaxYVA_j(1);
    PolicyPar.tauYSimp_j = TaxY_j.*tauS_tauY_j;      % SIMPLES Value added tax rates
    %TaxY_j = [tY_j(1),0.5152,0.1140];%[tY_j(1),0.4362,0.1139];%[tY_j(1),0.4290,0.1719];%[6.2,25.9,25]/100;       % Stantard value added tax rate [fixed for j=1, initial guess = tY_j]
    %tauYSimpForm_j = [1,0.4,0.7336];%[1,1,1]*0.65;%% Value added tax rates [discussion below]
else
    %TaxY_j = (1./(1-Data.InfVA_j)).* Data.TaxYVA_j;    % Standard tax rate
    %TaxY_j(1) = Data.TaxYVA_j(1);    
end
%%
J = SO.NumberofSectors;
onesJ = ones(1,J);
zerosJ = zeros(1,J);
%
FixPar.J = J;
FixPar.jE = jE; 
FixPar.EForm_j = EForm_j;
FixPar.LmaxInf_j = LmaxInf*onesJ;
EstPar.xi_j = xi_j;
if SO.prodknown == 1
    sigmaz_j = [1,1,1];
end
FixPar.sigmaz_j = sigmaz_j;
%

%% Production 
FixPar.retscale_j = 1 - (1-Data.WVA_j).*vash';
FixPar.theta_j = FixPar.retscale_j .* (onesJ - vash');
FixPar.alpha_j = FixPar.retscale_j .* vash';
FixPar.lambda_jk = ioiish';
FixPar.kappaForm_j = Data.ExitRate_j;                       % Discount rate, formal standard taxation firms
% Informal
FixPar.kappaInf_j = kappaInfForm_j.*FixPar.kappaForm_j;     % Discount rate, informal firms
if SO.FirmOptions > 2 % Add SIMPLES
    % SIMPLES
    FixPar.kappaSimp_j = kappaSimpForm_j.*FixPar.kappaForm_j;   % Discount rate, SIMPLES firms 
    FixPar.ESimp_j = ESimpForm_j.*FixPar.EForm_j; 
    %
    FixPar.EInf_j = EInfForm_j.*((1-VAs_VAF_j).*FixPar.EForm_j+VAs_VAF_j.*FixPar.ESimp_j);                 % Entry sunk cost, informal firms
else
    FixPar.EInf_j = EInfForm_j.*FixPar.EForm_j;                 % Entry sunk cost, informal firms
end

%% Prices: initial steady state 
Prices.w = w;
Prices.p_j = p_j;
%Prices.pE = Prices.p_j(jE);
Prices.pE = w;
Prices = TR_Piota(Prices,FixPar);                       % piota_j, intermediate inputs price index

%% Simulated firms
load(FN.firmsimuls, 'X1', 'X2', 'N')
EstPar.Nf = N;                                          % Number of simulated firms
EstPar.X1_f = X1(1:N,11);                               % Simulated draws for pre-entry productivity
EstPar.X2_f = X2(1:N,11);                               % Simulated draws for post-entry productivity

%% Productivity distribution   

%%% Productivity lower bound
% Formal firms
tauY_j = TaxY_j;
tauW = TaxW;
Fz_Form = TR_Emptoz(FixPar,Prices,tauW,tauY_j);
zmin_Form_j = onesJ;
%zmax_Form_j = onesJ;
% Informal firms
tauY_j = zeros(1,J);
tauW = 0;
Fz_Inf = TR_Emptoz(FixPar,Prices,tauW,tauY_j);
zmin_Inf_j = onesJ;
lzlow_j = zerosJ;
lzup_j = zerosJ;
Xm_j = zerosJ;


% Bound
lzlow_j(1) = 1;
Xm_j(1) = 1;
for j=2:J
    zmin_Form_j(j) = Fz_Form(j,1);
    %zmax_Form_j(j) = Fz_Form(j,Data.Ne);
    zmin_Inf_j(j) = Fz_Inf(j,1);
    lzminj = log(min(zmin_Form_j(j),zmin_Inf_j(j)));
    lzlow_j(j) = floor(lzminj*10)/10;         % Productivity grid, lower bound (normalization)
    %Xm_j(j) = exp(lzlow_j(j));
end
lzlow = min(lzlow_j);
lzlow_j = [1,lzlow,lzlow];
Xm_j = [1,exp(lzlow),exp(lzlow)];
%
FixPar.Xm_j = Xm_j;
%lzmin = log(min(min(zmin_Form_j(2:end)),min(zmin_Inf_j(2:end))));
%lzmax = log(max(zmax_Form_j(2:end)));
%
%lzlow = floor(lzmin*10)/10;         % Productivity grid, lower bound (normalization)
%lzup = ceil(lzmax*10000)/10000;           % Productivity grid, upper bound COMPLETAR AUTOMATIZAR
%lzup = floor(lzmax*10)/10;

%%% Simulated firms by sector
Nf = EstPar.Nf;
Nf_j = floor(Nf/J) * onesJ;
Nf_j(1) = Nf - sum(Nf_j(2:end));
EndoPar.Nf_j = Nf_j;                                    % Number of simulated firms by sector
X1_fj = cell(1,J);
X2_fj = cell(1,J);
jn = 0;
EndoPar.lz0_fj = cell(1,J);
EndoPar.lz_fj = cell(1,J);
for j = 1:J
    j0 = jn + 1;
    jn = Nf_j(j) + j0 - 1;
    X1_fj{j} = EstPar.X1_f(j0:jn);
    X2_fj{j} = EstPar.X2_f(j0:jn);

    % Pre-entry productivity (log scale), simulated firms
    z0 = Xm_j(j) *((ones(Nf_j(j),1) - X1_fj{j}).^(-1/EstPar.xi_j(j)));
    lz0_f = log(z0);
    lzup_j(j) = max(lz0_f)*1.05;

    % Keep parameters
    EndoPar.lz0_fj{j} = lz0_f;
end

%%% Productivity grid
FixPar.lz_grid_j = cell(1,J);
j=1;FixPar.lz_grid_j{j} = lzlow_j(j);   % Agricultura
for j=2:J
    dz = 0.05;                          % Distance between z grid points
    lz_grid = (lzlow_j(j):dz:lzup_j(j))';
    while length(lz_grid) < 101         % Productivity grid with more than 100 points
        dz = dz/2;
        lz_grid = (lzlow_j(j):dz:lzup_j(j))';
    end    
    FixPar.lz_grid_j{j} = lz_grid;    
end


% Post-entry x pre-entry probabilities
lz_grid_j = FixPar.lz_grid_j;
sigmaz_j = FixPar.sigmaz_j;
%
Probzz_j = cell(1,J);
for j=2:J
    lz_grido = lz_grid_j{j};
    Nz = length(lz_grido);                                  % Number of points in productivity grid
    if SO.prodknown == 0
        Probzz = zeros(Nz,Nz);                                  % Probability of productivity values after entry, conditional on the signal before entry
        for iz=1:Nz   
            q = lz_grido(iz); 
            Probzz(iz,:) = PreEntryProb(q,sigmaz_j(j),dz,lz_grido);
        end
    else
        Probzz = eye(Nz);  % IMPONDO PRE ENTRY = POST ENTRY    
    end
    EndoPar.Probzz_j{j} = Probzz;                                % 'Prob' in Ulyssea (2018) code
end

for j = 1:J
    % Post-entry productivity (log scale), simulated firms
    if SO.prodknown == 0
        Epsilon = sigmaz_j(j)*X2_fj{j};                         % The post-entry shock, mean 'mu_j' = 0
    else
        Epsilon = 0*sigmaz_j(j)*X2_fj{j};                         % The post-entry shock, mean 'mu_j' = 0    
    end
    lz_f = EndoPar.lz0_fj{j} + Epsilon;
    
    % Keep parameters
    EndoPar.lz_fj{j} = lz_f;
end

%%

EndoPar.lA = log(A); % log da produtividade da agricultura

return
