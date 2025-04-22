%//////////////////////////////////////////////////////////////////////////
% Indirect Tax Reform
%//////////////////////////////////////////////////////////////////////////
% TR_Equilibrium
% Find equilibrium
%//////////////////////////////////////////////////////////////////////////

function [Prices,Totals_j,Msim_j,Gov,Hous,walras] = TR_FindEquilibrium(SO,FixPar,PolicyPar,EndoPar,Prices0,shnew,toler,jmm,is)

%%% Parameters
Lbar = FixPar.Lbar;
J = FixPar.J;
jE = FixPar.jE;
zerosJ = zeros(1,J);

%%% Prices loop
%'Prices0' = initial guess for prices
Prices = Prices0;
COUNT=0;
tol=1;

while (tol>toler)
    %% Prices
    w = Prices.w;
    p_j = Prices.p_j;
    piota_j = Prices.piota_j;
    pE = Prices.pE;
    
    %% Firms
    [Totals_j,Msim_j] = TR_SimEconomy(SO,FixPar,PolicyPar,EndoPar,Prices);

    %% Government and GDP

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
    Gov.TaxYSimp_TaxY_j = TauYSimp_j./(TauY_j+TauYSimp_j);

    % Real GDP
    Gov.GDP = GDP;
    Gov.InfGDP = GDPinf / GDP;
    Gov.SimpGDP = GDPsimp / GDP;

    
    %% Households
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
    Wages = w * Lbar;    
    Income = Wages + Profits + Gov.T;    
    c_j = Income * FixPar.zeta_j ./ p_j;
    Hous.wagesGDP = Wages / GDP;
    Hous.profitsGDP = Profits / GDP;
    Hous.trGDP = Gov.T / GDP;
    Hous.incomeGDP = Income / GDP;
    Hous.pc_jGDP = p_j .*c_j / GDP;
    Hous.pcGDP = sum(p_j .* c_j) / GDP;
    Hous.c_j = c_j;
        
    %% Market clearing
    
    %%% Labor mkt clearing       
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

    LD = 0;
    for j=1:J
        LD = LD + Totals_j{j}.L;                        
    end
    LD = LD + Ecost/pE;
    Hous.LD = LD;
    tf = (LD - Lbar) / Lbar;
    wn = (1 + tf/2) * w;
    Prices.w = (1-shnew)* w + shnew* wn;
        
        
    %%% Goods mkt clearing, good J = numeraire
    mktgood_j = zeros(1,J);
    Q_j = zeros(1,J);       % gross output of good j
    Iota_j = zeros(1,J);    % demand for good j as intermediate input
    for k=1:J
        for j=1:J
            Iota_j(k) = Iota_j(k) + Totals_j{j}.Iota_k(k);
        end
        Q_j(k) = Totals_j{k}.Q;
    end
    
        
    for j = 1:J
        mktgood_j(j) = c_j(j) + Iota_j(j) - Q_j(j);
    end
    %mktgood_j(jmm) = mktgood_j(jmm) + Ecost/pE;
        
    jj = 1:j;    
    jj0 = jj(jj ~= jmm);    
    tf_j = mktgood_j(jj0) ./ Q_j(jj0);
    pn_j = (ones(size(tf_j)) + tf_j/2) .* p_j(jj0);
    Prices.p_j(jj0) = (1-shnew)* p_j(jj0) + shnew* pn_j;    
    %Prices.pE = Prices.p_j(jE);
    Prices.pE = Prices.w;
    Prices = TR_Piota(Prices,FixPar);                    %piota_j, intermediate inputs price index
        
    % Tolerance
    tol = abs(tf) + sum(abs(tf_j))
    ttf = tf
    ttf_j = tf_j
    COUNT=COUNT+1

    if is > 2
        parad=1;
    end
    
    if COUNT > 100
        return
    end

       
end

%% Check Walras' law
walras = mktgood_j(jmm) / Q_j(jmm)
assert(abs(walras) < 1e-2, "Good J market does not clear")
end


