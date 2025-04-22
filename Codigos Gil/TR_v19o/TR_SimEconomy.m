%//////////////////////////////////////////////////////////////////////////
% Indirect Tax Reform
%//////////////////////////////////////////////////////////////////////////
% TR_SimEconomy
% Simulates the economy for a given wage and parameter set
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [Totals_j, Msim_j]=TR_SimEconomy(SO,FixPar,PolicyPar,EndoPar,Prices)

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% PARAMETERS, GRIDS, AND PRICES
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% All firms
wage        = Prices.w;
p_j         = Prices.p_j;
piota_j     = Prices.piota_j;
pE          = Prices.pE;
J           = FixPar.J;
lz_grid_j   = FixPar.lz_grid_j;
lambda_jk   = FixPar.lambda_jk;
theta_j     = FixPar.theta_j;
alpha_j     = FixPar.alpha_j;
Probzz_j    = EndoPar.Probzz_j;
Nf_j        = EndoPar.Nf_j;
lz0_fj      = EndoPar.lz0_fj;
lz_fj       = EndoPar.lz_fj;

% Formal firms
kappaForm_j     = FixPar.kappaForm_j;
EForm_j         = pE*FixPar.EForm_j;
tauY_j          = PolicyPar.tauY_j;
tauW            = PolicyPar.tauW;
if SO.FirmOptions > 2 % Including SIMPLES firms
    kappaSimp_j = FixPar.kappaSimp_j;
    ESimp_j     = pE*FixPar.ESimp_j;
    tauYSimp_j  = PolicyPar.tauYSimp_j;
    tauWSimp    = PolicyPar.tauWSimp;
    RmaxSimp_j  = PolicyPar.RmaxSimp_j;
end

if SO.FirmOptions > 1 % Including informal firms
    % Informal firms
    kappaInf_j      = FixPar.kappaInf_j;
    EInf_j          = pE*FixPar.EInf_j;
    LmaxInf_j       = FixPar.LmaxInf_j;
end

% Sectoral loop
Totals_j = cell(1,J);
Msim_j = cell(1,J);

for j = 1   % Agriculture
    
    %%% Productivity 
    lz_grid = EndoPar.lA;

    %%% Sector-specific parameters
    p = p_j(j);
    piota = piota_j(j);
    alpha = alpha_j(j);
    theta = theta_j(j);
    p_k = p_j;
    lambda_k = lambda_jk(j,:);
    tauY = tauY_j(j);

    % variables not used for j=1
    kappaForm = kappaForm_j(j);
    Probzz = Probzz_j{j};
    lz0_f = lz_grid;
    lz_f = lz_grid;
    Nf = Nf_j(j);
    EForm = EForm_j(j);
    
    % LABOR DEMAND    
    LaborPolicy = TR_LaborFunc(SO,alpha,theta,tauY,tauW,lz_grid,wage,p,piota,j,0,1,1);
    
    % Value functions   
    [Firms,~] = TR_ValueFunc(SO,wage,p,piota,p_k,lambda_k,alpha,theta,...
            tauY,tauW,kappaForm,Probzz,lz_grid,lz0_f,lz_f,LaborPolicy,j,0,0,0,0,0);
    
    % Aggregates:
    Totals = TR_AggregateFunc(SO,Firms,lz_f,Nf,EForm,kappaForm,0,0,0,0,j);    
    
    % Moments
    Msim = 0;
    
    %% Output
    Totals_j{j} = Totals;
    Msim_j{j} = Msim;
    
    clear  alpha theta kappaForm tauY 
            
end

for j = 2:J
    
    %%% Productivity 
    lz_grid = lz_grid_j{j};
    
    %%% Sector-specific parameters
    p = p_j(j);
    piota = piota_j(j);
    alpha = alpha_j(j);
    theta = theta_j(j);
    p_k = p_j;
    lambda_k = lambda_jk(j,:);
    kappaForm = kappaForm_j(j);
    EForm = EForm_j(j);
    tauY = tauY_j(j);
    if SO.FirmOptions > 2 % Including SIMPLES firms
        kappaSimp = kappaSimp_j(j);
        ESimp = ESimp_j(j);
        tauYSimp = tauYSimp_j(j);
        RmaxSimp = RmaxSimp_j(j);
    end
    if SO.FirmOptions > 1 % Including informal firms
        kappaInf = kappaInf_j(j);
        EInf = EInf_j(j);
        LmaxInf = LmaxInf_j(j);
    end
    Nf = Nf_j(j);
    lz0_f = lz0_fj{j};
    lz_f = lz_fj{j};
    Probzz = Probzz_j{j};

    %%
    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    % LABOR DEMAND
    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    if SO.FirmOptions == 3       % Including SIMPLES firms
        LaborPolicy = TR_LaborFunc(SO,alpha,theta,tauY,tauW,lz_grid,wage,p,piota,j,LmaxInf,tauYSimp,tauWSimp);
    elseif SO.FirmOptions == 2   % Including informal firms
        LaborPolicy = TR_LaborFunc(SO,alpha,theta,tauY,tauW,lz_grid,wage,p,piota,j,LmaxInf,1,1);
    elseif SO.FirmOptions == 1
        LaborPolicy = TR_LaborFunc(SO,alpha,theta,tauY,tauW,lz_grid,wage,p,piota,j,0,1,1);
    end
    
    %%
    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    % Value functions
    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    %%%%%%%%%%%%%%%% SIMPLES, arrumar daqui em diante
    if SO.FirmOptions == 1          % Only formal firms, main tax code
        [Firms,~] = TR_ValueFunc(SO,wage,p,piota,p_k,lambda_k,alpha,theta,...
            tauY,tauW,kappaForm,Probzz,lz_grid,lz0_f,lz_f,LaborPolicy,j,0,0,0,0,0);
    elseif SO.FirmOptions == 2      % Including informal firms
        [Firms,~] = TR_ValueFunc(SO,wage,p,piota,p_k,lambda_k,alpha,theta,...
            tauY,tauW,kappaForm,Probzz,lz_grid,lz0_f,lz_f,LaborPolicy,j,kappaInf,0,0,0,0);        
    elseif SO.FirmOptions == 3      % Including SIMPLES firms
        [Firms,LaborPolicy] = TR_ValueFunc(SO,wage,p,piota,p_k,lambda_k,alpha,theta,...
            tauY,tauW,kappaForm,Probzz,lz_grid,lz0_f,lz_f,LaborPolicy,j,kappaInf,tauYSimp,tauWSimp,kappaSimp,RmaxSimp);
    end
    PEND=0;
    PEND=PEND+1; % a) Completar gráficos do SIMPLES; 
    PEND=PEND+1; % ARRUMAR gráficos profits (principalmente SIMPLES), Colocar pontos de corte nos gráficos separados de tipos de firmas
    
    %%
    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    % Aggregates:
    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    if SO.FirmOptions == 1          % Only formal firms, main tax code
        Totals = TR_AggregateFunc(SO,Firms,lz_f,Nf,EForm,kappaForm,0,0,0,0,j);
    elseif SO.FirmOptions == 2      % Including informal firms
        Totals = TR_AggregateFunc(SO,Firms,lz_f,Nf,EForm,kappaForm,EInf,kappaInf,0,0,j);
    elseif SO.FirmOptions == 3      % Including SIMPLES firms
        Totals = TR_AggregateFunc(SO,Firms,lz_f,Nf,EForm,kappaForm,EInf,kappaInf,ESimp,kappaSimp,j);
    end
    PEND=PEND+1; % i) ARRUMAR toda a parte do SIMPLES, desativada aqui; ii) fazer limpeza geral aqui, muito texto desativado
        
    %% Moments
    Msim = TR_GenMoments(SO,Totals,j);
    PEND=PEND+1; % i) AJUSTAR momentos à estimação; ii) acrescentar SIMPLES
    
    %% Output
    Totals_j{j} = Totals;
    Msim_j{j} = Msim;
    
    clear  alpha theta kappaForm tauY 
    
end

end
