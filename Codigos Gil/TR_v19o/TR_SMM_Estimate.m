function [Q]=TR_SMM_Estimate(EstParam,Prices,FixPar,PolicyPar,EndoPar,Data,SO,FN)
% Simulate the economy for a given parameter
% returns (MS-MD)W(MS-MD) -> objective function

% Parameters
J = FixPar.J;
jE = FixPar.jE;
zerosJ = zeros(1,J);
is = 1;

% Prices
w = Prices.w;
p_j = Prices.p_j;

%% Unstack EstPar into FixPar,EndoPar,etc
%%%
k=1;
% var_single
A = exp(EstParam(k));k=k+1;
tauS_tauY2 = exp(EstParam(k));k=k+1;
%xi_j(2) = exp(EstParam(k));k=k+1;
%sigmaz_j(2) = exp(EstParam(k));k=k+1;
%RmaxSimp_j(2) = exp(EstParam(k));k=k+1;
%% var_2sec    
EForm_j(2:3) = exp(EstParam(k:k+1));k=k+2;
xi_j(2:3) = exp(EstParam(k:k+1));k=k+2;
sigmaz_j(2:3) = exp(EstParam(k:k+1));k=k+2;
kappaInfForm_j(2:3) = exp(EstParam(k:k+1));k=k+2;
kappaSimpForm_j(2:3) = exp(EstParam(k:k+1));k=k+2;
TaxY_j(2:3) = exp(EstParam(k:k+1));k=k+2;
%EInfForm_j(2:3) = exp(EstParam(k:k+1));k=k+2;
%ESimpForm_j(2:3) = exp(EstParam(k:k+1));k=k+2;
%RmaxSimp_j(2:3) = exp(EstParam(k:k+1));k=k+2;
%% Excluded from estimation
%kappaInfForm_j = [1,1,1];
%kappaSimpForm_j = [1,1,1];
%xi_j(3) = 4.3941;
%sigmaz_j(3) = 0.5334;
%RmaxSimp_j(3) = 110.8102;
EInfForm_j = [1,1,1]*0.47;
ESimpForm_j = [1,1,1]*0.5;
RmaxSimp_j = [1,1,1]*114;
%%%
%% Bounds
%RmaxSimp_j(RmaxSimp_j>114) = 114;
%RmaxSimp_j(RmaxSimp_j<114/2) = 114/2;
EForm_j(EForm_j>10) = 10;
EForm_j(EForm_j<5) = 5;
TaxY_j(TaxY_j>0.49)=0.49;
TaxY_j(TaxY_j<0.13)=0.13;
%ESimpForm_j(ESimpForm_j>0.8) = 0.8;
%EInfForm_j = min(ESimpForm_j,EInfForm_j);


%%%
%% Remaining parameters 
TR_SetParam
TR_PolicyParam

% e.g.
% FixPar.EForm_j(2:3) = abs(EstPar(1:2));
% FixPar.kappaInf_j(2:3) = abs(EstPar(3:4));
% FixPar.EInf_j(2:3) = abs(EstPar(5:6));
% FixPar.kappaSimp_j(2:3) = abs(EstPar(7:8));
% FixPar.ESimp_j(2:3) = abs(EstPar(9:10));
% FixPar.sigmaz_j(2:3) = abs(EstPar(11:12));
% xi???
%EndoPar.lA = abs(EstPar(15));
%PolicyPar.tauY_j(2:3) = abs(EstPar(16:17));

% Simulate Economy
[Totals_j,Msim_j] = TR_SimEconomy(SO,FixPar,PolicyPar,EndoPar,Prices);

% Labor market clearing
LD = 0;
L_j = zerosJ;
for j=1:J
    L_j(j) = Totals_j{j}.L;
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
    VA_j = zerosJ;
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
        VA_j(j) = Totals_j{j}.VA; 
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
    Gov.VA_j = VA_j;
    Gov.VAj_Lj = VA_j./L_j;
    
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
     %% Employment distribution
    % Formal firms population
    empForm_f_j = cell(1,J);
    NTForm_j = ones(1,J);
    NSimp_j = zeros(1,J);
    for j=2:J
        NFormj = Totals_j{j}.NumEntryH;             % Number of entrants in sector j
        Emphj0 = round(Totals_j{j}.Lh);
        emph_fj = sort(Emphj0(Emphj0>0));           % Employees by firm in sector j 
        if SO.FirmOptions == 3                      % Add SIMPLES firms
            NSimpj = Totals_j{j}.NumEntryS;
            NTFormj = NFormj + NSimpj;
            Empsj0 = round(Totals_j{j}.Ls);
            emps_fj = sort(Empsj0(Empsj0>0));        
            empForm_fj = sort([emps_fj;emph_fj]);   % Formal: standard taxation (h) + SIMPLES (s)
            %Data.empForm_f_j;
            NSimp_j(j) = NSimpj;
        else
            NTFormj = NFormj;
            empForm_fj = emph_fj;
        end
        NTForm_j(j) = NTFormj;
        empForm_f_j{j} = empForm_fj;            % Formal firm size population, by secotr j
    end
    
    % Formal employment distribution
    nue = Data.nu0e;                            % Fit Pareto distribution for firms with number of employees >= nue
    Ne = Data.Ne;    
    xie_j = zeros(1,J);
    %
    theta_gp = nue-0.0001;                      % 'fitdist(...,'gp') requires observed values strictly greater than 'theta'
    shTForm_0tonue = zeros(1,J); 
    for j=2:J
        shTForm_0tonue(j) = length(empForm_f_j{j}(empForm_f_j{j}<nue)) / length(empForm_f_j{j});
        %gpj = fitdist(empForm_f_j{j}(empForm_f_j{j}>theta_gp), 'gp','theta',theta_gp);
        %K=gpj.k;
        %xie_j(j) = 1/K;
        empForm_f = empForm_f_j{j}(empForm_f_j{j}>theta_gp);
        N = length(empForm_f);
        grid = min(empForm_f):1:max(empForm_f);
        xie_j(j) = N/sum(log(empForm_f/nue));
    end


% Data  and Model Moments
% each one must be a vector

k=1;

% Mean employees by firm, all formal and SIMPLES firms [target: w]
%MD(k) = sum(Data.empForm_j(2:3)) / sum(Data.Nf_j(2:3));
%MS(k) = (Totals_j{2}.LH+Totals_j{3}.LH) / sum(NTForm_j(2:3));
%%MS(k) = (Totals_j{2}.LH+Totals_j{2}.LS+Totals_j{3}.LH+Totals_j{3}.LS) / sum(NTForm_j(2:3));
%WW(k) = 1;
%k=k+1;

%VA_L_TFormal_j = VA_TFormal_j ./ L_TFormal_j;
MS(k) = Gov.VAj_Lj(3)/Gov.VAj_Lj(2);
MD(k) = Data.VAj_Lj(3)/Data.VAj_Lj(2);
WW(k) = 100;
k=k+1;

% GDP shares by sector [target: p_j/p_jE]
MD(k:k+2) = Data.VAj_GDP;
MS(k:k+2) = [Totals_j{1}.VA,Totals_j{2}.VA,Totals_j{3}.VA] / GDP;
WW(k:k+2) = 20;
k=k+3;

% Mean employees by firm, services / manufacturing ratio
%MS(k) = (Totals_j{3}.LH / NTForm_j(3)) / (Totals_j{2}.LH / NTForm_j(2));
%%MS(k) = ((Totals_j{3}.LH + Totals_j{3}.LS) / NTForm_j(3)) / ((Totals_j{2}.LH + Totals_j{2}.LS) / NTForm_j(2));
%MD(k) = (Data.empForm_j(3) / Data.Nf_j(3)) / (Data.empForm_j(2) / Data.Nf_j(2));
%WW(k) = 1;
%k=k+1;

% Share of firms with less than 'nue' employees [target: EForm_j/EForm_jE]
%MD(k:k+1) = Data.empForm_cdf_j(nue-1,2:3);
%MS(k:k+1) = shTForm_0tonue(2:3);
%WW(k:k+1) = 1;
%k=k+2;

% Indirect taxesj / VAj
MD(k:k+1) = Data.TaxYVA_j(2:3);
MS(k:k+1) = Gov.TaxYVA_j(2:3);
WW(k:k+1) = 40;
k=k+2;

% Share of informal value added in sector value added [target: EInf_j/EForm_j]
MD(k:k+1) = Data.InfVA_j(2:3);
MS(k:k+1) = Gov.InfGDP_j(2:3);
WW(k:k+1) = 70;
k=k+2;

% Share of SIMPLES in tax records by sector [target: ESimp_j/EForm_j]
MD(k:k+1) = [0.018,0.089];
MS(k:k+1) = TauYSimp_j(2:3)./(TauY_j(2:3)+TauYSimp_j(2:3));
WW(k:k+1) = 70;
k=k+2;

% Share of SIMPLES in number of total formal firms by sector [target: sigma]
%MD(k:k+1) = [0.68,0.55];
%MS(k:k+1) = NSimp_j(2:3)./NTForm_j(2:3);

% xi_employment, 'nue' or more employees [target: xi_j] 
%MD(k:k+1) = Data.xi0e_j(2:3);
%MS(k:k+1) = xie_j(2:3);
%WW(k:k+1) = 1;
%k=k+2;


% Weight matrix
%W = eye(length(MD));
%W(1,1) = 30;
W = diag(WW);

Q=(MS-MD)*W*(MS-MD)';

% Display results
ParamA = exp(EstParam(1))'
Param = exp(EstParam(2:end))'
Mdat = MD
Mmod = MS

end