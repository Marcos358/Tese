
%%%
%% Variables determined in SMM
var_single  = {A,tauS_tauY2};
var_2sec    = {EForm_j,xi_j,sigmaz_j,kappaInfForm_j,kappaSimpForm_j,TaxY_j};
%var_single  = {A,tauS_tauY2,xi_j(2),sigmaz_j(2),RmaxSimp_j(2)};
%var_2sec    = {EForm_j,TaxY_j,EInfForm_j,ESimpForm_j};
%%%

%% EstPar0 (initial guess)
% Must enter TR_SM_Estimate as vector (not structure)

N1v = length(var_single);
N2v = length(var_2sec);

EstPar0 = [];
for jv = 1:N1v
    EstPar0 = [EstPar0;log(var_single{jv})];
end
for jv = 1:N2v
    EstPar0 = [EstPar0;log(var_2sec{jv}(2:3))'];
end

% Call TR_SMM_Estimate (Objective Function)
% Function that returns (MS-MD)W(MS-MD)'

%Q = TR_SMM_Estimate(EstPar0,Prices,FixPar,PolicyPar,EndoPar,Data,SO);%,SO,FN,FixPar,PolicyPar,EndoPar,Prices,checkfig,N,X1,X2,D,J,MD);

% Minimize Objective Function
options = optimset('Display','on','PlotFcns',@optimplotfval,...
    'MaxFunEvals', 1e4, 'MaxIter', 1e4);

Estimates = fminsearch(@(EstParam)TR_SMM_Estimate(EstParam,Prices,FixPar,...
                        PolicyPar,EndoPar,Data,SO,FN),EstPar0,options);
%%%
k=1;
%% var_single
A = exp(Estimates(k));k=k+1;
tauS_tauY2 = exp(Estimates(k));k=k+1;
%xi_j(2) = exp(Estimates(k));k=k+1;
%sigmaz_j(2) = exp(Estimates(k));k=k+1;
%RmaxSimp_j(2) = exp(Estimates(k));k=k+1;
%% var_2sec    
EForm_j(2:3) = exp(Estimates(k:k+1));k=k+2;
xi_j(2:3) = exp(Estimates(k:k+1));k=k+2;
sigmaz_j(2:3) = exp(Estimates(k:k+1));k=k+2;
kappaInfForm_j(2:3) = exp(Estimates(k:k+1));k=k+2;
kappaSimpForm_j(2:3) = exp(Estimates(k:k+1));k=k+2;
TaxY_j(2:3) = exp(Estimates(k:k+1));k=k+2;
%EInfForm_j(2:3) = exp(Estimates(k:k+1));k=k+2;
%ESimpForm_j(2:3) = exp(Estimates(k:k+1));k=k+2;
%RmaxSimp_j(2:3) = exp(Estimates(k:k+1));k=k+2;
%%%

%% Remaining parameters 
TR_SetParam
TR_PolicyParam


% Firms
[Totals_j,Msim_j] = TR_SimEconomy(SO,FixPar,PolicyPar,EndoPar,Prices);

% Labor market clearing
LD = 0;
for j=1:J
    LD = LD + Totals_j{j}.L;                        
end
Lbar = LD;
FixPar.Lbar = Lbar;

    % Goods market clearing
    c_j = zerosJ;
    Q_j = zerosJ;       % gross output of good j
    Iota_j = zerosJ;    % demand for good j as intermediate input
    for k=1:J
        for j=1:J
            Iota_j(k) = Iota_j(k) + Totals_j{j}.Iota_k(k);
        end
        Q_j(k) = Totals_j{k}.Q;
    end
    
    Ecost = 0;   % entry costs
    if SO.FirmOptions == 1
        for k=2:J
            Ecost = Ecost + Totals_j{k}.Ecost_Form;
        end    
    elseif SO.FirmOptions == 2      % add informal firms    
        for k=2:J
            Ecost = Ecost + Totals_j{k}.Ecost_Form + Totals_j{k}.Ecost_Inf;
        end    
    elseif SO.FirmOptions == 3      % add SIMPLES firms    
        for k=2:J
            Ecost = Ecost + Totals_j{k}.Ecost_Form + Totals_j{k}.Ecost_Inf + Totals_j{k}.Ecost_Simp;
        end    
    end
    %DH_j = p_j .* Q_j;      % Domestic demand for good j [VERSÃO ECO FECHADA]
    
    for j = 1:J
        c_j(j) = Q_j(j) - Iota_j(j);
    end
    c_j(jE) = c_j(jE) - Ecost/Prices.pE;
    
    % Government
    Gov.name = "Government, revenues and transfers";
    TauW_j = zerosJ;
    TauY_j = zerosJ;
    TauWSimp_j = zerosJ;
    TauYSimp_j = zerosJ;
    TaxWVA_j = zerosJ;
    TaxYVA_j = zerosJ;
    InfGDP_j = zerosJ;
    SimpGDP_j = zerosJ;
    GDP = 0;
    GDPinf = 0;
    GDPsimp = 0;
    for j = 1:J
        TauW_j(j) = PolicyPar.tauW * w * Totals_j{j}.LH;
        TauY_j(j) = PolicyPar.tauY_j(j) * Totals_j{j}.VAh; 
        TaxWVA_j(j) = TauW_j(j) / Totals_j{j}.VA;
        TaxYVA_j(j) = TauY_j(j) / Totals_j{j}.VA;
        GDP = GDP + Totals_j{j}.VA;
        if j > 1
            if SO.FirmOptions > 1       % add informal firms
                GDPinf = GDPinf + Totals_j{j}.VAi;
                InfGDP_j(j) = Totals_j{j}.VAi / Totals_j{j}.VA;
                if SO.FirmOptions > 2   % add SIMPLES firms
                    TauWSimp_j(j) = PolicyPar.tauWSimp * w * Totals_j{j}.LS;
                    TauYSimp_j(j) = PolicyPar.tauYSimp_j(j) * Totals_j{j}.VAs; 
                    TaxWVA_j(j) = (TauW_j(j)+TauWSimp_j(j)) / Totals_j{j}.VA;
                    TaxYVA_j(j) = (TauY_j(j)+TauYSimp_j(j)) / Totals_j{j}.VA;
                    GDPsimp = GDPsimp + Totals_j{j}.VAs;
                    SimpGDP_j(j) = Totals_j{j}.VAs / Totals_j{j}.VA;
                end 
            end 
        end  
    end
    TauW = sum(TauW_j);
    TauY = sum(TauY_j);
    if SO.FirmOptions > 2   % add simples firms
        TauW = TauW + sum(TauWSimp_j);
        TauY = TauY + sum(TauYSimp_j); 
    end
    Gov.T = TauW + TauY;
    %Gov.TauWSimp = TauWSimp;
    %Gov.TauYSimp = TauYSimp;
    Gov.TaxWVA_j = TaxWVA_j;
    Gov.TaxYVA_j = TaxYVA_j;
    Gov.TaxWGDP = TauW / GDP;
    Gov.TaxYGDP = TauY / GDP;
    Gov.InfGDP_j = InfGDP_j;
    Gov.SimpGDP_j = SimpGDP_j;
    % Real GDP
    Gov.GDP = GDP;
    Gov.InfGDP = GDPinf / GDP;
    Gov.SimpGDP = GDPsimp / GDP;
    
    % Households
    Hous.name = "Households, GDP shares";
    Profits = 0;
    if SO.FirmOptions == 1      % Only formal firms
        for j=1:J
            Profits = Profits + Totals_j{j}.Profits_Form; 
        end
    elseif SO.FirmOptions == 2  % Formal and informal firms
        Profits = Profits + Totals_j{1}.Profits_Form;
        for j=2:J
            Profits = Profits + Totals_j{j}.Profits_Form + Totals_j{j}.Profits_Inf; 
        end
    elseif SO.FirmOptions == 3  % Formal, informal, and SIMPLES firms
        Profits = Profits + Totals_j{1}.Profits_Form;
        for j=2:J
            Profits = Profits + Totals_j{j}.Profits_Form + Totals_j{j}.Profits_Inf + Totals_j{j}.Profits_Simp; 
        end   
    end
    %
    Wages = w * Lbar;
    Income = Wages + Profits + Gov.T;
    zeta_j = p_j .* c_j ./ Income 
    FixPar.zeta_j = zeta_j;
    %
    Hous.wagesGDP = Wages / GDP;
    Hous.profitsGDP = Profits / GDP;
    Hous.trGDP = Gov.T / GDP;
    Hous.incomeGDP = Income / GDP;
    Hous.pc_jGDP = p_j .*c_j / GDP;
    Hous.pcGDP = sum(p_j .* c_j) / GDP;
    Hous.c_j = c_j;
    %Hous.wages_wagesGOS = Wages / (Wages + Profits + Ecost);
    Hous.Ecost_GOS = Ecost / (Profits + Ecost);
    Hous.Consumption_FinalDemand_jE = c_j(jE) / (Q_j(jE) - Iota_j(jE));
    %
    assert(abs(Ecost/GDP + Hous.incomeGDP -1) < 1e-3)
    %
