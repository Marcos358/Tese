%%% Show results

%%%%%
% Set baseline scenario and industry
ib = 1;     % [1] before tax reform, [2] after tax reform
jb = 2;     % j = {1,2,3}
%%%%%

% Copy results
RT.name = 'Results table';
for is = SO.SimScenarios
    gov = Results.Gov{is};
    totals = Results.Totals{is};
    %
    RT.GDP(is) = gov.GDP;
    RT.TaxYGDP(is) = gov.TaxYGDP;
    RT.Maintax(is) = Results.TauY{is}(jb);
    if is == 2
        RT.Maintax(is) = SO.AfterRefMainTax;
    end
    %
    RT.Taxrate_j(is,:) = Results.TauY{is};
    RT.TaxYVA_j(is,:) = gov.TaxYVA_j;
    if SO.FirmOptions > 1           % add informal firms
        RT.InfshGDP(is) = gov.InfGDP;
        RT.InfVA_j(is,:) = gov.InfGDP_j;
        if SO.FirmOptions > 2       % add SIMPLES firms
            RT.SimpshGDP(is) = gov.SimpGDP;
            RT.SimpVA_j(is,:) = gov.SimpGDP_j;  
            RT.TaxrateSimp_j(is,:) = Results.TauYSimp{is};
            RT.TaxYSimp_TaxY_j(is,:) = Results.TaxYSimp_TaxY_j{is};
        end
    end
    for j=1:J
        RT.VAGDP_j(is,j) = totals{j}.VA / gov.GDP;
        RT.VA_j(is,j) = totals{j}.VA;
        RT.L_j(is,j) = totals{j}.L;
        RT.VAj_Lj(is,j) = totals{j}.VA/totals{j}.L;
        RT.Yh_j(is,j) = totals{j}.Yh;
        if j > 1
            RT.NumEntryH_j(is,j) = totals{j}.NumEntryH;
            if SO.FirmOptions > 1   % add informal firms
                RT.Yi_j(is,j) = totals{j}.Yi;
                RT.NumEntryI_j(is,j) = totals{j}.NumEntryI;
            end
            if SO.FirmOptions > 2   % add SIMPLES firms
                RT.Ys_j(is,j) = totals{j}.Ys;
                RT.NumEntryS_j(is,j) = totals{j}.NumEntryS;
            end
        end
    end
end


% Normalization
GDP0 = RT.GDP(ib);
VA0_j = RT.VA_j(ib,:);

%%% Real value added
if SO.FirmOptions == 1
    Y_j    = RT.Yh_j;
    shY_j  = RT.VAGDP_j;
elseif SO.FirmOptions == 2      % Informal
    Y_j    = [RT.Yh_j, RT.Yi_j];
    shY_j  = [(1-RT.InfVA_j).*RT.VAGDP_j, RT.InfVA_j.*RT.VAGDP_j];
    RT.Y_is_j = zeros(max(SO.SimScenarios),J);
    RT.Y_is_j(is,1) = 1;
    for j=2:J
        % by j, aggregate formal and informal
        Yj      = [RT.Yh_j(:,j), RT.Yi_j(:,j)];
        shYj    = [1-RT.InfVA_j(:,j),RT.InfVA_j(:,j)];
        Nk      = size(Yj,2);
        Ybj     = Yj(ib,:);
        shYbj   = shYj(ib,:);
        for is = SO.SimScenarios
            Y_k = zeros(2,Nk);
            Y_k(1,:) = Ybj;
            Y_k(2,:) = Yj(is,:);
            shY_k = zeros(2,Nk);
            shY_k(1,:) = shYbj;
            shY_k(2,:) = shYj(is,:);
            Y_agg = TR_Tornqvist(Y_k,shY_k);
            RT.Y_is_j(is,j) = Y_agg(2);
        end
    end
elseif SO.FirmOptions == 3    % SIMPLES
    Y_j    = [RT.Yh_j, RT.Yi_j, RT.Ys_j];
    shY_j  = [(1-RT.InfVA_j-RT.SimpVA_j).*RT.VAGDP_j, RT.InfVA_j.*RT.VAGDP_j, RT.SimpVA_j.*RT.VAGDP_j];
%shY_j(2,end-1)=1e-10 %GAMBIARRA>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.
    RT.Y_is_j = zeros(max(SO.SimScenarios),J);
    RT.Y_is_j(is,1) = 1;
    for j=2:J
        % by j, aggregate formal, informal, and simples
        Yj      = [RT.Yh_j(:,j), RT.Yi_j(:,j), RT.Ys_j(:,j)];
        shYj    = [1-RT.InfVA_j(:,j)-RT.SimpVA_j(:,j),RT.InfVA_j(:,j),RT.SimpVA_j(:,j)];
        Nk      = size(Yj,2);
        Ybj     = Yj(ib,:);
        shYbj   = shYj(ib,:);
        for is = SO.SimScenarios
            Y_k = zeros(2,Nk);
            Y_k(1,:) = Ybj;
            Y_k(2,:) = Yj(is,:);
            shY_k = zeros(2,Nk);
            shY_k(1,:) = shYbj;
            shY_k(2,:) = shYj(is,:);
            Y_agg = TR_Tornqvist(Y_k,shY_k);
            RT.Y_is_j(is,j) = Y_agg(2);
        end
    end

end
% Economy aggregate
Yb_j   = Y_j(ib,:);
shYb_j = shY_j(ib,:);
Yb_j = Yb_j(shYb_j>0);
shYb_j  = shYb_j(shYb_j>0);
Nk = size(Yb_j,2);
for is = SO.SimScenarios
    Y_k = zeros(2,Nk);
    Y_k(1,:) = Yb_j;
    shY_k = zeros(2,Nk);
    shY_k(1,:) = shYb_j;
    shYo_j = shY_j(is,:);
    shY_k(2,:) = shYo_j(shYo_j>0)
    Yo_j = Y_j(is,:);
    Y_k(2,:) = Yo_j(shYo_j>0)
    
    Y_agg = TR_Tornqvist(Y_k,shY_k);
    RT.Y_is(is) = Y_agg(2);
end

% %clc
% if SO.SetEndoParams==1    % Parametrization
%     %%% ARRUMAR ISSO AQUI
%     % Table A + C
%     %disp(['______________________________________________________________'])
%     %disp(['        Formal, informal and SIMPLES firms, by tax rate       '])
%     %disp(['______________________________________________________________'])
%     %disp([' '])
%     %is = 1;
%     %disp(['                             Before  tax reform               '])
%     %disp(['______________________________________________________________'])
%     %disp([' '])
%     %disp(['Sector                     1       2       3      Total               '])
%     %var=RT.Taxrate_j;cc=3;disp(['Tax ratesj           ' '     ' num2str(100*var(is,1),cc) '%    ' num2str(100*var(is,2),cc)   '%   ' num2str(100*var(is,3),cc)   '%'])
%     %var=RT.VAGDP_j  ;cc=2;disp(['VAj / GDP            ' '     ' num2str(100*var(is,1),cc) '%    ' num2str(100*var(is,2),cc+1) '%   ' num2str(100*var(is,3),cc+1) '%    ' num2str(100)  '%'])
%     %disp([' '])
%     
%     disp([' '])
%     disp([' '])
%     disp([' '])
%     
%     % Table E
%     is = 1;
%     disp(['______________________________________________________________'])
%     disp(['            Moments: model (before tax reform) x data         '])
%     disp(['______________________________________________________________'])
%     disp([' '])
%     disp(['Sector                      1       2       3      Total               '])
%     disp([' '])
%     disp(['Tax ratesj'])
%     var=RT.Taxrate_j    ;cc=3;disp(['             Model' '         ' num2str(100*var(is,1),cc) '%   ' num2str(100*var(is,2),cc) '%   ' num2str(100*var(is,3),cc)  '%    '  ])
%     disp([' '])
%     disp(['Simples tax ratesj'])
%     var=RT.TaxrateSimp_j;cc=3;disp(['             Model' '           '                         '     ' num2str(100*var(is,2),cc) '%   ' num2str(100*var(is,3),cc)  '%    '  ])
%     disp([' '])
%     disp(['Informal VAj / VAj'])
%     var=RT.InfVA_j   ;cc=2;disp(['             Model' '         '  '        ' num2str(100*var(is,2),cc) '%     ' num2str(100*var(is,3),cc)  '%      '  num2str(100*RT.InfshGDP(is),cc) '%'])
%     dat=Data.InfVA_j ;cc=2;disp(['             Data ' '         '  '        ' num2str(100*dat(2),cc)    '%     '   num2str(100*dat(3),cc)   '%      '  num2str(100*Data.InfGDP(is),cc) '%'])
%     disp([' '])
%     disp(['Simples VAj / VAj'])
%     var=RT.SimpVA_j  ;cc=2;disp(['             Model' '         '  '        ' num2str(100*var(is,2),cc) '%     ' num2str(100*var(is,3),cc)  '%      '  num2str(100*RT.InfshGDP(is),cc) '%'])
%     %dat=Data.SimpVA_j ;cc=2;disp(['             Data ' '         '  '        ' num2str(100*dat(2),cc)    '%     '   num2str(100*dat(3),cc)   '%      '  num2str(100*Data.InfGDP(is),cc) '%'])
%     disp([' '])
%     disp(['Indirect taxesj / VAj'])
%     var=RT.TaxYVA_j  ;cc=2;disp(['             Model' '        ' num2str(100*var(is,1),cc) '%     ' num2str(100*var(is,2),cc) '%     ' num2str(100*var(is,3),cc)  '%     '  num2str(100*RT.TaxYGDP(is),cc) '%'])
%     dat=Data.TaxYVA_j;cc=2;disp(['             Data ' '        ' num2str(100*dat(1),cc)    '%     ' num2str(100*dat(2),cc)    '%     ' num2str(100*dat(3),cc)     '%     '  num2str(100*Data.TaxGDP(is),cc) '%'])
%     disp([' '])
%     disp(['VAj / Lj'])
%     var0=(RT.GDP(is)/sum(RT.L_j(is,:)));dat0=Data.GDP_L;
%     var=RT.VAj_Lj    ;cc=3;disp(['             Model' '         ' num2str(100*var(is,1)/var0,cc) '    ' num2str(100*var(is,2)/var0,cc) '     ' num2str(100*var(is,3)/var0,cc) '     ' num2str(100*var0/var0,cc)    ' '])
%     dat=Data.VAj_Lj  ;cc=3;disp(['             Data ' '         ' num2str(100*dat(1)/dat0,cc-1)  '      '  num2str(100*dat(2)/dat0,cc)   '     ' num2str(100*dat(3)/dat0,cc)     '     ' num2str(100*dat0/dat0,cc+1)  ' '])
%     disp([' '])
%     disp(['VAj / GDP'])
%     var=RT.VAGDP_j  ;cc=2;disp(['             Model' '         ' num2str(100*var(is,1),cc-1) '%      ' num2str(100*var(is,2),cc) '%     ' num2str(100*var(is,3),cc)  '%'])
%     dat=Data.VAj_GDP;cc=2;disp(['             Data ' '         ' num2str(100*dat(1),cc-1)    '%      ' num2str(100*dat(2),cc)    '%     '   num2str(100*dat(3),cc)   '%'])
%     disp([' '])
%     disp(['______________________________________________________________'])
%     return
%     kapratio = FixPar.kapratio
%     %zeta_j = FixPar.zeta_j 
%     %EForm = FixPar.EForm_j
% 
% else
% 
% %if SO.FirmOptions == 2
% if SO.FirmOptions >1
% %is = tax scenario - [1] before tax reform, [2] single tax rate, [3] after tax reform
% 
% % Table A
% 
% disp(['______________________________________________________________'])
% disp(['                    Tax Reform in Brazil                      '])
% disp(['                  Formal and informal firms                   '])
% disp(['______________________________________________________________'])
% disp(['                           Before       After           '])
% disp(['Variable                 tax reform   tax reform        '])
% disp([' '])
% disp(['Main tax rate              ' num2str(100*RT.Maintax(1),3)  '%        ' num2str(100*RT.Maintax(2),3)  '%'])
% disp(['Indirect taxes / GDP       ' num2str(100*RT.TaxYGDP(1),3)  '%        ' num2str(100*RT.TaxYGDP(2),3)  '%'])
% disp(['Informal GDP / GDP         ' num2str(100*RT.InfshGDP(1),3) '%        ' num2str(100*RT.InfshGDP(2),3) '%'])
% disp(['GDP                        ' num2str(100*RT.Y_is(1),3)     '          ' num2str(100*RT.Y_is(2),4)    ' '])
% disp(['______________________________________________________________'])
% %num2str(100*RT.GDP(3)/GDP0,4)
% 
% disp([' '])
% disp([' '])
% disp([' '])
% 
% % Table C 
% disp(['______________________________________________________________'])
% disp(['            Formal and informal firms, by tax rate            '])
% disp(['______________________________________________________________'])
% disp([' '])
% disp(['Industry                   1.Agro  2.Manuf  3.Serv            '])
% disp(['______________________________________________________________'])
% disp([' '])
% disp(['                            Tax rates               '])
% var=RT.Taxrate_j;
% disp([' '])
% is=1;cc=3;disp(['Before tax reform    ' '       ' num2str(100*var(is,1),cc)   '%    ' num2str(100*var(is,2),cc) '%   ' num2str(100*var(is,3),cc) '%'])
% is=2;cc=3;disp(['After tax reform     ' '       ' num2str(100*var(is,1),cc-1) '%    ' num2str(100*var(is,2),cc) '%   ' num2str(100*var(is,3),cc) '%'])
% disp([' '])
% disp(['                            Indirect taxesj / VAj               '])
% var=RT.TaxYVA_j;
% disp([' '])
% is=1;cc=3;disp(['Before tax reform    ' '       ' num2str(100*var(is,1),cc-1) '%    ' num2str(100*var(is,2),cc) '%    ' num2str(100*var(is,3),cc) '%'])
% is=2;cc=3;disp(['After tax reform     ' '       ' num2str(100*var(is,1),cc-1) '%    ' num2str(100*var(is,2),cc) '%   ' num2str(100*var(is,3),cc) '%'])
% disp([' '])
% disp(['                            Informal VAj / GDPj               '])
% var=RT.InfVA_j;
% disp([' '])
% is=1;cc=3;disp(['Before tax reform    ' '       ' num2str(100*var(is,1),cc) '%    ' num2str(100*var(is,2),cc) '%   ' num2str(100*var(is,3),cc) '%'])
% is=2;cc=3;disp(['After tax reform     ' '       ' num2str(100*var(is,1),cc) '%    ' num2str(100*var(is,2),cc-1) '%   ' num2str(100*var(is,3),cc) '%'])
% disp([' '])
% if SO.FirmOptions > 2
% disp(['                            SIMPLES VAj / GDPj               '])
% var=RT.SimpVA_j;
% disp([' '])
% is=1;cc=3;disp(['Before tax reform    ' '       ' num2str(100*var(is,1),cc) '%    ' num2str(100*var(is,2),cc) '%   ' num2str(100*var(is,3),cc) '%'])
% is=2;cc=3;disp(['After tax reform     ' '       ' num2str(100*var(is,1),cc) '%    ' num2str(100*var(is,2),cc-1) '%   ' num2str(100*var(is,3),cc) '%'])
% disp([' '])
% disp(['                            SIMPLES taxesj / taxesj               '])
% var=RT.TaxYSimp_TaxY_j;
% disp([' '])
% is=1;cc=3;disp(['Before tax reform    ' '       ' num2str(100*var(is,1),cc) '%    ' num2str(100*var(is,2),cc) '%   ' num2str(100*var(is,3),cc) '%'])
% is=2;cc=3;disp(['After tax reform     ' '       ' num2str(100*var(is,1),cc) '%    ' num2str(100*var(is,2),cc-1) '%   ' num2str(100*var(is,3),cc) '%'])
% disp([' '])
% end
% disp(['                            VAj / GDP               '])
% var=RT.VAGDP_j;
% disp([' '])
% is=1;cc=3;disp(['Before tax reform    ' '       ' num2str(100*var(is,1),cc-1) '%     ' num2str(100*var(is,2),cc) '%   ' num2str(100*var(is,3),cc) '%'])
% is=2;cc=3;disp(['After tax reform     ' '       ' num2str(100*var(is,1),cc-1) '%   ' num2str(100*var(is,2),cc) '%   ' num2str(100*var(is,3),cc) '%'])
% disp([' '])
% disp(['______________________________________________________________'])
% 
% disp([' '])
% disp([' '])
% disp([' '])
% 
% % Table D
% disp(['______________________________________________________________'])
% disp(['            Formal and informal firms, GDP change             '])
% disp(['______________________________________________________________'])
% disp(['                             Value added by j            '])
% disp([' '])
% disp(['Industry                   1.Agro  2.Manuf  3.Serv            '])
% disp([' '])
% var=RT.Y_is_j;cc=3;
% is = 1;
% disp(['Before tax reform     ' '      ' num2str(100*var(is,1),cc) '    ' num2str(100*var(is,2),cc) '    ' num2str(100*var(is,3),cc) ' '])
% is = 2;
% disp(['After tax reform      ' '      ' num2str(100*var(is,1),cc) '    ' num2str(100*var(is,2),cc) '    ' num2str(100*var(is,3),cc) ' '])
% disp(['______________________________________________________________'])
% 
% disp([' '])
% disp([' '])
% disp([' '])
% 
%     % Table E
%     is = 1;
%     disp(['______________________________________________________________'])
%     disp(['            Moments: model (before tax reform) x data         '])
%     disp(['______________________________________________________________'])
%     disp([' '])
%     disp(['Industry                1.Agro  2.Manuf  3.Serv   Total               '])
%     disp([' '])
%     disp(['Informal VAj / VAj'])
%     var=RT.InfVA_j   ;cc=2;disp(['             Model' '        ' num2str(100*var(is,1),cc) '%     ' num2str(100*var(is,2),cc) '%     ' num2str(100*var(is,3),cc)  '%     '  num2str(100*RT.InfshGDP(is),cc) '%'])
%     dat=Data.InfVA_j ;cc=2;disp(['             Data ' '        ' num2str(100*dat(1),cc)    '%     ' num2str(100*dat(2),cc)    '%     '   num2str(100*dat(3),cc)   '%     '  num2str(100*Data.InfGDP(is),cc) '%'])
%     disp([' '])
%     disp(['Indirect taxesj / VAj'])
%     var=RT.TaxYVA_j  ;cc=2;disp(['             Model' '        ' num2str(100*var(is,1),cc) '%    ' num2str(100*var(is,2),cc) '%    ' num2str(100*var(is,3),cc)  '%     '  num2str(100*RT.TaxYGDP(is),cc) '%'])
%     dat=Data.TaxYVA_j;cc=2;disp(['             Data ' '        ' num2str(100*dat(1),cc)    '%    ' num2str(100*dat(2),cc)    '%    '   num2str(100*dat(3),cc)   '%     '  num2str(100*Data.TaxGDP(is),cc) '%'])
%     disp([' '])
%     disp(['VAj / Lj'])
%     var0=(RT.GDP(is)/sum(RT.L_j(is,:)));dat0=Data.GDP_L;
%     var=RT.VAj_Lj    ;cc=3;disp(['             Model' '        ' num2str(100*var(is,1)/var0,cc) '     ' num2str(100*var(is,2)/var0,cc) '      ' num2str(100*var(is,3)/var0,cc-1) '     ' num2str(100*var0/var0,cc)    ' '])
%     dat=Data.VAj_Lj  ;cc=3;disp(['             Data ' '         ' num2str(100*dat(1)/dat0,cc-1) '     '  num2str(100*dat(2)/dat0,cc)       '     ' num2str(100*dat(3)/dat0,cc)      '     ' num2str(100*dat0/dat0,cc+1)  ' '])
%     disp([' '])
%     disp(['______________________________________________________________'])
% 
%         disp([' '])
%     disp([' '])
%     disp([' '])
%     
%     disp(['______________________________________________________________'])
%     disp(['                      Number of firms                         '])
%     disp(['______________________________________________________________'])
%     disp([' '])
%     disp(['Sector                     2.Manuf   3.Serv                  '])
%     disp([' '])
%     cc=6;
%     disp(['Formal:'     ])
%     var=RT.NumEntryH_j;
%     is = 1;
%     disp(['Before tax reform     ' '         ' num2str(var(is,2),cc) '    ' num2str(var(is,3),cc) ' '])
%     is = 2;
%     disp(['After tax reform      ' '        ' num2str(var(is,2),cc) '    ' num2str(var(is,3),cc) ' '])
%     disp(['Informal:'   ])
%     var=RT.NumEntryI_j;
%     is = 1;
%     disp(['Before tax reform     ' '      ' num2str(var(is,2),cc) '    ' num2str(var(is,3),cc) ' '])
%     is = 2;
%     disp(['After tax reform      ' '      ' num2str(var(is,2),cc) '    ' num2str(var(is,3),cc) ' '])
%    if SO.FirmOptions > 2
%     disp(['SIMPLES:'    ])
%     var=RT.NumEntryS_j;
%     is = 1;
%     disp(['Before tax reform     ' '        ' num2str(var(is,2),cc) '    ' num2str(var(is,3),cc) ' '])
%     is = 2;
%     disp(['After tax reform      ' '        ' num2str(var(is,2),cc) '    ' num2str(var(is,3),cc) ' '])
%    end
%     disp([' '])
%     disp(['______________________________________________________________'])
%     
% Results.Totals{is} = Totals_j;
% %elseif SO.FirmOptions == 1
% else
%     
% % Table B
% disp(['______________________________________________________________'])
% disp(['                    Tax Reform in Brazil                      '])
% disp(['                     Only formal firms                        '])
% disp(['______________________________________________________________'])
% disp(['                           Before       After           '])
% disp(['Variable                 tax reform   tax reform        '])
% disp([' '])
% disp(['Main tax rate              ' num2str(100*RT.Maintax(1),3)  '%        ' num2str(100*RT.Maintax(2),3)  '%'])
% disp(['Indirect taxes / GDP       ' num2str(100*RT.TaxYGDP(1),3)  '%        ' num2str(100*RT.TaxYGDP(2),3)  '%'])
% disp(['GDP                        ' num2str(100*RT.Y_is(1),3)     '          ' num2str(100*RT.Y_is(2),4)    ' '])
% disp(['______________________________________________________________'])
% 
% disp([' '])
% disp([' '])
% disp([' '])
% 
% 
% % Table C 
% disp(['______________________________________________________________'])
% disp(['            Formal and informal firms, by tax rate            '])
% disp(['______________________________________________________________'])
% disp([' '])
% disp(['Industry                   1.Agro  2.Manuf  3.Serv            '])
% disp(['______________________________________________________________'])
% disp([' '])
% disp(['                            Tax rates               '])
% var=RT.Taxrate_j;
% disp([' '])
% is=1;cc=3;disp(['Before tax reform    ' '       ' num2str(100*var(is,1),cc)   '%    ' num2str(100*var(is,2),cc) '%   ' num2str(100*var(is,3),cc) '%'])
% is=2;cc=3;disp(['After tax reform     ' '       ' num2str(100*var(is,1),cc-1) '%    ' num2str(100*var(is,2),cc) '%   ' num2str(100*var(is,3),cc) '%'])
% disp([' '])
% disp(['                            Indirect taxesj / VAj               '])
% var=RT.TaxYVA_j;
% disp([' '])
% is=1;cc=3;disp(['Before tax reform    ' '       ' num2str(100*var(is,1),cc-1) '%    ' num2str(100*var(is,2),cc) '%    ' num2str(100*var(is,3),cc) '%'])
% is=2;cc=3;disp(['After tax reform     ' '       ' num2str(100*var(is,1),cc-1) '%    ' num2str(100*var(is,2),cc) '%   ' num2str(100*var(is,3),cc) '%'])
% disp([' '])
% disp(['                            VAj / GDP               '])
% var=RT.VAGDP_j;
% disp([' '])
% is=1;cc=3;disp(['Before tax reform    ' '       ' num2str(100*var(is,1),cc-1) '%     ' num2str(100*var(is,2),cc) '%   ' num2str(100*var(is,3),cc) '%'])
% is=2;cc=3;disp(['After tax reform     ' '       ' num2str(100*var(is,1),cc-1) '%   ' num2str(100*var(is,2),cc) '%   ' num2str(100*var(is,3),cc) '%'])
% disp([' '])
% disp(['______________________________________________________________'])
% 
% disp([' '])
% disp([' '])
% disp([' '])
% 
% 
% 
% 
% end
% end
% 
% 
% 
% 
% 
% return
