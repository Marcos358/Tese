
%%%
pareto_fig = 0;                             % plot figures for cheking employment distribution: yes [1] or not [0]
%%%

% Parameters
J = FixPar.J;
jE = FixPar.jE;
zerosJ = zeros(1,J);

% Prices
w = Prices.w;
p_j = Prices.p_j;

    % Firms
    [Totals_j,Msim_j] = TR_SimEconomy(SO,FixPar,PolicyPar,EndoPar,Prices);

    % Entry costs
    Ecost = 0;   
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

    % Labor market clearing
    LD = 0;
    L_j = zerosJ;
    for j=1:J
        L_j(j) = Totals_j{j}.L;
        LD = LD + L_j(j);                      
    end
    Lbar = LD + Ecost/Prices.pE;
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
        
    for j = 1:J
        c_j(j) = Q_j(j) - Iota_j(j);
    end
    %c_j(jE) = c_j(jE) - Ecost/Prices.pE;
    
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
    Gov.TaxYSimp_TaxY_j = TauYSimp_j./(TauY_j+TauYSimp_j);
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
    %Hous.Ecost_GOS = Ecost / (Profits + Ecost);
    %Hous.Consumption_FinalDemand_jE = c_j(jE) / (Q_j(jE) - Iota_j(jE));
    Hous.Consumption_FinalDemand_jE = 0;
    %
    %assert(abs(Ecost/GDP + Hous.incomeGDP -1) < 1e-3)
    %
    Gov
    Hous
    %Data_wages_wagesGOS = Data.w_GDP
    
    
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

        if pareto_fig == 1
            X=nue:Ne;
            THETA=theta_gp;
            K_data = 1/Data.xi0e_j(j);
            SIGMA_data = K_data*THETA;
            Y_data = gppdf(X,K_data,SIGMA_data,THETA);
            %
            K_model = 1/xie_j(j);
            SIGMA_model = K_model*THETA;
            Y_model = gppdf(X,K_model,SIGMA_model,THETA);
            nne = nue:100; % grid for figures
            figure;scatter(X(1:nne(end)),Y_data(1:nne(end)));hold on;scatter(X(1:nne(end)),Y_model(1:nne(end))); 
        end
    end
    
    
    %% DATA X MODEL
    % OBS: TForm = formal taxação padrão + SIMPLES [Total Formal]
    
    % Consumption share in industry final demand [target: EForm_jE]
    %Model_Consumption_FinalDemand_jE = Hous.Consumption_FinalDemand_jE*0
    %Data_Consumption_FinalDemand_jE = Data.Consumption_FinalDemand_jE*0
    
    % Mean employees by firm, all formal and SIMPLES firms [target: w]
    if SO.FirmOptions > 2
        Model_meanLfirm = (Totals_j{2}.LH+Totals_j{2}.LS+Totals_j{3}.LH+Totals_j{3}.LS) / sum(NTForm_j(2:3))
        Data_meanLfirm = sum(Data.empForm_j(2:3)) / sum(Data.Nf_j(2:3))
        % by industry
        Model_meanLfirm_j = [0,Totals_j{2}.LH+Totals_j{2}.LS,Totals_j{3}.LH+Totals_j{3}.LS]./ NTForm_j
        Data_meanLfirm_j = Data.empForm_j ./ Data.Nf_j
    else
        Model_meanLfirm = (Totals_j{2}.LH+Totals_j{3}.LH) / sum(NTForm_j(2:3))
        Data_meanLfirm = sum(Data.empForm_j(2:3)) / sum(Data.Nf_j(2:3))
        % by industry
        Model_meanLfirm_j = [0,Totals_j{2}.LH,Totals_j{3}.LH]./ NTForm_j
        Data_meanLfirm_j = Data.empForm_j ./ Data.Nf_j
    end        
    % GDP shares by sector [target: p_j/p_jE]
    Model_GDPsh_j = [Totals_j{1}.VA,Totals_j{2}.VA,Totals_j{3}.VA] / GDP
    Data_GDPsh_j = Data.VAj_GDP
    
    % Share of firms with less than 'nue' employees [target: EForm_j/EForm_jE]
    Model_shTForm_lessthannue = shTForm_0tonue
    Data_shTForm_lessthannue = Data.empForm_cdf_j(nue-1,:)
    
    % xi_employment, 'nue' or more employees [target: xi_j] 
    Model_xie_j = xie_j
    Data_xi0e_j = Data.xi0e_j
    
    % Share of number of total formal firms by sector
    % [comparar com RAIS, usar como target sobreidentificado ou como momento not targeted]
    Model_shNTForm_j = NTForm_j/sum(NTForm_j)
    Data_shNTForm_j = Data.Nf_j/sum(Data.Nf_j(2:end))
    
    %%%% Informal
    
    % Share of informal value added in sector value added [target: EInf_j/EForm_j]
    Model_InfVA_j = Gov.InfGDP_j
    Data_InfVA_j = Data.InfVA_j
        
    % share of informal employment in total employment, firms with 1 to 5 employees  ---> sigma_j
    %Data_shLinf_j
    %Model_shLinf_j
    % OU share da informalidade no PIB
    
    %%%% SIMPLES
    if SO.FirmOptions > 2
        % Share of SIMPLES in number of total formal firms by sector [target: sigma]
        Model_shfirms_SimpTForm_j = NSimp_j./NTForm_j
        Data_shfirms_SimpTForm_j = [0,0.6847,0.552]%[0,0.63,0.59]        % RAIS, AUTOMATIZAR ---->GIL
        
        % Share of SIMPLES in tax records by sector [target: ESimp_j/EForm_j]
        Model_TaxYSimp_TaxY_j = Gov.TaxYSimp_TaxY_j
        Data_TaxYSimp_TaxY_j = [0,0.018,0.089]%[0,0.021,0.09]           % RAIS, AUTOMATIZAR ---->GIL
        Data.TaxYSimp_TaxY_j = Data_TaxYSimp_TaxY_j;
    end
    
    % formal firms, relative labor productivity (labor / value added), industry / services  ---> 
    VA_TFormal_j = zeros(1,J);
    L_TFormal_j = zeros(1,J);
    for j=1:J
        VA_TFormal_j(j) = Totals_j{j}.VAh;
        L_TFormal_j(j) = Totals_j{j}.LH;            
        if (SO.FirmOptions > 2) && (j>1)      % add SIMPLES firms
            VA_TFormal_j(j) = VA_TFormal_j(j) + Totals_j{j}.VAs;
            L_TFormal_j(j) = L_TFormal_j(j) + Totals_j{j}.LS;            
        end
    end
    VA_L_TFormal_j = VA_TFormal_j ./ L_TFormal_j;
    %
    Model_FormalLaborProductivity2_3 = VA_L_TFormal_j(2)/VA_L_TFormal_j(3)
    
    
    %%% MOMENTS
    % Ecost médio: parear "Ecost_GOS" com "Data_Investment_GOS"
    % A_1: parear com share da agricultura no PIB
    % xi_j: parear xie, pareto do emprego para employees >=6, no modelo e nos dados
    % z_min: normalização
        
    
    % Table 
    is = 1;
    disp(['______________________________________________________________'])
    disp(['            Moments: model (before tax reform) x data         '])
    disp(['______________________________________________________________'])
    disp([' '])
    disp(['                                            Model  Data                '])
    disp([' '])
    %var=Model_Consumption_FinalDemand_jE;
    %dat=Data_Consumption_FinalDemand_jE ;cc=2;
    %disp(['Manufacturing consumption / final demand:' '    ' num2str(100*var,cc) '%    ' num2str(100*dat,cc) '%     '])
    %disp([' '])
    var=round(Model_meanLfirm);
    dat=round(Data_meanLfirm) ;cc=5;
    disp(['Mean employees by firm, formal and SIMPLES:' '  ' num2str(var,cc) '    ' num2str(dat,cc) '     '])
    disp([' '])
    disp(['______________________________________________________________'])

    disp([' '])
    disp([' '])
    disp([' '])
    
    is = 1;
    disp(['______________________________________________________________'])
    disp(['     Moments, all firms: model (before tax reform) x data     '])
    disp(['______________________________________________________________'])
    disp([' '])
    disp(['            Sector    1.Agro  2.Manuf  3.Serv                  '])
    disp([' '])
    disp(['VAj / GDP'])
    var=Model_GDPsh_j ;cc=2;disp(['             Model' '       ' num2str(100*var(1),cc) '%   ' num2str(100*var(2),cc) '%     ' num2str(100*var(3),cc)  '%     '  ])
    dat=Data_GDPsh_j  ;cc=2;disp(['             Data ' '       ' num2str(100*dat(1),cc) '%     ' num2str(100*dat(2),cc) '%     ' num2str(100*dat(3),cc)  '%     '  ])
    disp([' '])
    disp(['VAj / Lj, relative to j=2'])
    var=Gov.VAj_Lj/Gov.VAj_Lj(jE)    ;cc=3;disp(['             Model' '       ' num2str(100*var(1)*0,cc) '%     ' num2str(100*var(2),cc) '%     ' num2str(100*var(3),cc)  '%     '  ])  
    dat=Data.VAj_Lj/Data.VAj_Lj(jE)  ;cc=3;disp(['             Data ' '       ' num2str(100*dat(1)*0,cc) '%     ' num2str(100*dat(2),cc) '%     ' num2str(100*dat(3),cc)  '%     '  ])
    disp([' '])
    disp(['______________________________________________________________'])
%return
    disp([' '])
    disp([' '])
    disp([' '])

    disp(['______________________________________________________________'])
    disp([' Moments, all formal firms: model (before tax reform) x data  '])
    disp(['______________________________________________________________'])
    disp([' '])
    disp(['            Sector    1.Agro  2.Manuf  3.Serv                  '])
    disp([' '])
    disp(['Number of firms with 1 to 5 employess / number of firms'])
    var=Model_shTForm_lessthannue;cc=2;disp(['             Model' '       ' num2str(100*var(1),cc) '%     ' num2str(100*var(2),cc) '%     ' num2str(100*var(3),cc)  '%     '  ])
    dat=Data_shTForm_lessthannue ;cc=2;disp(['             Data ' '       ' num2str(100*dat(1),cc) '%     ' num2str(100*dat(2),cc) '%     ' num2str(100*dat(3),cc)  '%     '  ])
    disp([' '])
    disp(['Firm size Pareto distribution, shape parameter, 6 or more employees'])
    var=Model_xie_j;cc=3;disp(['             Model' '       ' num2str(var(1),cc) '     ' num2str(var(2),cc) '     ' num2str(var(3),cc)  '     '  ])
    dat=Data_xi0e_j;cc=3;disp(['             Data ' '       ' num2str(dat(1),cc) '     ' num2str(dat(2),cc) '     ' num2str(dat(3),cc)  '     '  ])
    disp([' '])
    disp(['Number of firms sector j/ number of firms'])
    var=Model_shNTForm_j;cc=2;disp(['             Model' '       ' num2str(100*var(1)*0,cc) '%     ' num2str(100*var(2),cc) '%     ' num2str(100*var(3),cc)  '%     '  ])
    dat=Data_shNTForm_j; cc=2;disp(['             Data ' '       ' num2str(100*dat(1)*0,cc) '%     ' num2str(100*dat(2),cc) '%     ' num2str(100*dat(3),cc)  '%     '  ])
    disp([' '])
    disp(['Mean employees by firm, formal and SIMPLES:'])
    var=Model_meanLfirm_j;cc=2;disp(['             Model' '       ' num2str(var(1)*0,cc) '     ' num2str(var(2),cc) '     ' num2str(var(3),cc)  '     '  ])
    dat=Data_meanLfirm_j; cc=2;disp(['             Data ' '       ' num2str(dat(1)*0,cc) '     ' num2str(dat(2),cc) '     ' num2str(dat(3),cc)  '     '  ])
    disp([' '])
    disp(['______________________________________________________________'])

    disp([' '])
    disp([' '])
    disp([' '])

    disp(['______________________________________________________________'])
    disp(['Moments, informal and SIMPLES: model (before tax reform) x data '])
    disp(['______________________________________________________________'])
    disp([' '])
    disp(['            Sector    1.Agro  2.Manuf  3.Serv                  '])
    disp([' '])
    disp(['Informal firms'])
    disp(['________________'])
    disp([' '])
    disp(['Indirect taxesj / VAj'])
    var=Gov.TaxYVA_j ;cc=2;disp(['             Model' '        ' num2str(100*var(1),cc) '%    ' num2str(100*var(2),cc) '%    ' num2str(100*var(3),cc)  '% '])
    dat=Data.TaxYVA_j;cc=2;disp(['             Data ' '        ' num2str(100*dat(1),cc) '%    ' num2str(100*dat(2),cc) '%    ' num2str(100*dat(3),cc)  '% '])
    disp([' '])% TOTAL:     '  num2str(100*Gov.TaxYGDP,cc) '%; '  num2str(100*Data.TaxGDP,cc) '%
    disp(['Informal VAj / VAj'])
    var=Model_InfVA_j;cc=2;disp(['             Model' '        ' num2str(100*var(1),cc) '%    ' num2str(100*var(2),cc) '%    ' num2str(100*var(3),cc)  '%'])
    dat=Data_InfVA_j ;cc=2;disp(['             Data ' '        ' num2str(100*dat(1),cc) '%    ' num2str(100*dat(2),cc) '%    ' num2str(100*dat(3),cc)  '%'])
    disp([' '])% TOTAL:      '  num2str(100*Gov.InfGDP,cc) '%; '  num2str(100*Data.InfGDP,cc) '%
   if SO.FirmOptions > 2
    disp(['SIMPLES firms'])
    disp(['________________'])
    disp([' '])
    disp(['Number of SIMPLES firms in j / Number of formal firms in j'])
    var=Model_shfirms_SimpTForm_j;cc=3;disp(['             Model' '        ' num2str(100*var(1),cc) '%    ' num2str(100*var(2),cc) '%    ' num2str(100*var(3),cc)  '%     '])
    dat=Data_shfirms_SimpTForm_j ;cc=2;disp(['             Data ' '        ' num2str(100*dat(1),cc) '%    ' num2str(100*dat(2),cc) '%    ' num2str(100*dat(3),cc)  '%     '])
    disp([' '])
    disp(['Production taxes in SIMPLES / Production taxes'])
    var=Model_TaxYSimp_TaxY_j;cc=3;disp(['             Model' '        ' num2str(100*var(1),cc) '%    ' num2str(100*var(2),cc) '%    ' num2str(100*var(3),cc)  '%     '])
    dat=Data_TaxYSimp_TaxY_j ;cc=2;disp(['             Data ' '        ' num2str(100*dat(1),cc) '%    ' num2str(100*dat(2),cc) '%    ' num2str(100*dat(3),cc)  '%     '])
    disp([' '])
   end
    disp(['______________________________________________________________'])

    disp([' '])
    disp([' '])
    disp([' '])
    
    disp(['______________________________________________________________'])
    disp(['                      Number of firms                         '])
    disp(['______________________________________________________________'])
    disp([' '])
    disp(['Sector     2.Manuf   3.Serv                  '])
    disp([' '])
    cc=6;
    disp(['Formal:'     '       ' num2str(Totals_j{2}.NumEntryH,cc) '      ' num2str(Totals_j{3}.NumEntryH,cc) '     '])
    disp(['Informal:'   '   ' num2str(Totals_j{2}.NumEntryI,cc) '    ' num2str(Totals_j{3}.NumEntryI,cc) '     '])
   if SO.FirmOptions > 2
    disp(['SIMPLES:'    '   ' num2str(Totals_j{2}.NumEntryS,cc) '    ' num2str(Totals_j{3}.NumEntryS,cc) '     '])
   end
    disp([' '])
    disp(['______________________________________________________________'])
    
