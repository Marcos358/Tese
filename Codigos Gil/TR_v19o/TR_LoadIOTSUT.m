%//////////////////////////////////////////////////////////////////////////
% Indirect Tax Reform
%//////////////////////////////////////////////////////////////////////////
% TR_LoadIOT
% Load Brazilian Input-output tables and compute shares
%//////////////////////////////////////////////////////////////////////////

clearvars -except FN PEND SO

%% Load Excel data on Input Output Tables
% Uses tables, basic prices, current values in R$ 1,000,000

%%%%%%%%%%%%%%%%%%%%
% Set options
%%%%%%%%%%%%%%%%%%%%
usesimp = 0;        % include imports in uses table [1] or home production only [0].
year0 = 2015;       % base year
%%%%%%%%%%%%%%%%%%%%
loadexcel = SO.LoadData;      % load Excel data [1] or use previously loaded values [0]

dataiotsut = fullfile(FN.codeoutputs,"iotdata.mat");        % file with IOT data in MATLAB format
if loadexcel == 1
    nameio = "Matriz_de_Insumo_Produto_"  + year0 + "_Nivel_67.xls";
    filetoread = fullfile(FN.dataIOT,nameio);
    q_sec  = xlsread(filetoread,'01',"H134:BV134");         % National production vector by sector
    Resources = xlsread(filetoread,'01',"H6:BW132");        % Resources of goods and services (products x sectors) 
    Supply = Resources(:,1:end-1)';
    q_prod = Resources(:,end);                              % National production vector by product
    Uses_Home = xlsread(filetoread,'03',"D6:BR132");        % Supply and demand of national production, basic prices (products x sectors) 
    Uses_Imported = xlsread(filetoread,'04',"D6:BR132");    % Supply and demand of imported products, basic prices (products x sectors)
    Imports = xlsread(filetoread,'04',"C6:C132");           % Imported products (products x sectors)
    Investment_Home = xlsread(filetoread,'03',"BX6:BX132");      % Investment, national production, basic prices (products x sectors) 
    Investment_Imported = xlsread(filetoread,'04',"BX6:BX132");  % Investment, imported products, basic prices (products x sectors)
    Exports_Home = xlsread(filetoread,'03',"BT6:BT132");         % Exports, national production, basic prices (products x sectors) 
    Exports_Imported = xlsread(filetoread,'04',"BT6:BT132");     % Exports, imported products, basic prices (products x sectors) - non zero in "airplanes, ships and transportation equipments"
    FinalDemand_Home = xlsread(filetoread,'03',"BZ6:BZ132");     % Final demand, national production, basic prices (products x sectors) 
    FinalDemand_Imported = xlsread(filetoread,'04',"BZ6:BZ132"); % Final demand, imported products, basic prices (products x sectors)
    IndTax_Home = xlsread(filetoread,'05',"C6:C132");            % Indirect taxes vector, national production (by product) 
    IndTax_Imported = xlsread(filetoread,'06',"C6:C132");        % Indirect taxes vector, imported products (by product)
    %Investment = xlsread(filetoread,'02',"BW134");              % Gross fixed capital formation
    %
    Bn = xlsread(filetoread,'11',"C6:BQ132");    % Technical coefficients matrix for national inputs (products x sectors) 
    Bm = xlsread(filetoread,'12',"C6:BQ132");    % Technical coefficients matrix for imported inputs (products x sectors)
    D  = xlsread(filetoread,'13',"C6:DY72");     % Market share matrix D (sectors x products) - sector shares in production of each product (column sums = 1)
    assert(max(abs(sum(D)-1))<1e-10)
    %
    namesut2 = "68_tab2_" + year0 + ".xls";
    filesut2 = fullfile(FN.dataIOT,namesut2);
    VA0 = xlsread(filesut2,'VA',"B6:BR19");                     % Value added, 68 industries
    VA = [VA0(:,1:40),sum(VA0(:,41:42),2),VA0(:,43:end-1)];     % Value added, 67 industries 
    VAtot = VA0(:,end);
    VA_sec = VA(1,:);
    L_sec = VA(end,:);
    %
    %nameprods = "Products_v3d.xlsx";
    nameprods = "Products_v3e.xlsx";
    if SO.NumberofSectors == 2
        nameprods = "Products_v2d.xlsx";
    end
    fileprods = fullfile(FN.dataIOT,nameprods);
    %
    Products = xlsread(fileprods,'Product',"J6:L132");
    group = Products(:,1);
    aftertaxfactor_prod = Products(:,2);
    %aftertax_prod = Products(:,3);
    %
    Industries = xlsread(fileprods,'DD',"EC6:EC72");
    indgroup = Industries;
    if SO.NumberofSectors == 2
        indgroup = indgroup - 1;
        indgroup(indgroup < 0) = 0;
    end    
    %
    labelInf0 = xlsread(fileprods,'VA',"B2:BQ2");                           % 68 industries
    labelInf = [labelInf0(1:40),labelInf0(41),labelInf0(:,43:end)];         % 67 industries
    Informal = xlsread(fileprods,'Informal',"B6:J18");
    shInfVA = 1-mean(Informal(:,2:3),2)/100;    % 1 - formal share on VA, mean 2014-2016 
    shInfL  = 1-mean(Informal(:,6:7),2)/100;    % 1 - formal share on employment, mean 2014-2016
    %
    labelExit0 = xlsread(fileprods,'VA',"B3:BQ3");                          % 68 industries
    labelExit = [labelExit0(1:40),labelExit0(41),labelExit0(:,43:end)];     % 67 industries
    Exit = xlsread(fileprods,'Exit',"C6:D25");                              % 
    ExitRate = Exit(:,1)/100;
    %
    assert(abs(VAtot(end-3)-79101)<1e-12)       % 79101 from "tab18.xls", 'other taxes on production', year 2015
    OtherPayrollTaxes = 38595;                  % "tab18.xls", 'Impostos sobre a folha de pagamento', year 2015
    Wages = VAtot(3);
    OtherTaxWw = OtherPayrollTaxes/Wages;
    %
    save(dataiotsut) 
else
    load(dataiotsut)
end

%% Check B and D tables

IndTax = IndTax_Home; % CHECAR, ignora imposto sobre produtos importados
ITBI_2015 = 10280;      % 'tab18.xls', 'Impostos sobre transmissão de bens imóveis inter-vivos'
i_Buildings = 90;       % product 'Buildings' 
IndTax(i_Buildings) = IndTax(i_Buildings) - ITBI_2015;

% Matrix B: products use by unit of sector production (Only for checking now. We compute B again after sector aggregation)
if usesimp == 1
    Exports = Exports_Home + Exports_Imported;
    Investment = Investment_Home + Investment_Imported; 
    FinalDemand = FinalDemand_Home + FinalDemand_Imported;
    Consumption = FinalDemand - Investment - Exports;
    Uses = Uses_Home + Uses_Imported;   % National + imported supply
    %IndTax = IndTax_Home + IndTax_Imported;
    B = Bn + Bm;
else
    Investment = Investment_Home; 
    FinalDemand = FinalDemand_Home;
    Consumption = FinalDemand - Investment;
    Uses = Uses_Home;                   % National supply 
    %IndTax = IndTax_Home;
    B = Bn;         
end
Bcheck = Uses*diag(1./q_sec);
assert(max(max(abs(Bcheck-B)))<1e-12)

% Market share matrix D (Only for checking now. We compute D again after sector aggregation)
Dcheck = Supply*diag(1./q_prod);
assert(max(max(abs(Dcheck-D)))<1e-12)

%% Informal production
Njj = length(labelInf);
shInfVA_j = zeros(1,Njj);
shInfL_j = zeros(1,Njj);
for j=1:Njj
    shInfVA_j(j) = shInfVA(labelInf(j));
    shInfL_j(j) = shInfL(labelInf(j));
end

% Inputed rents correction
%shImputedRents = 0.73;  % Share of imputed rents in real state production
%jRealState = 53;
%x = (shInfVA_j(jRealState) - shImputedRents)/(1-shImputedRents)

% Informal value added and employment by industry
VAInf_j = shInfVA_j.*VA_sec;
LInf_j = shInfL_j.*L_sec;

%% Exit rates
Nje = length(labelExit); 
ExitRate_j = zeros(1,Nje);
for j=1:Nje
    ExitRate_j(j) = ExitRate(labelExit(j));
end

%% Sector aggregation
load(FN.modelinputs);
indnum = SO.NumberofSectors;        % Number of industries
yearnum = 1;                        % Number of years

% selected products
allprod = 1:length(q_prod);
keepprod = allprod(group>0);

% selected industries
indcell = cell(indnum,1);
excludeind = length(indgroup(indgroup==0));
keepind = length(indgroup(indgroup>0));

    ie = 0;
    ik = 0;
    for jn=1:length(indgroup)
        ip = indgroup(jn);
        if ip > 0
            ik = ik + 1;
            keepind(ik) = jn;
            indcell{ip} = [indcell{ip}, jn];
        else
            ie = ie + 1;
            excludeind(ie) = jn;
        end
    end

Nprod = length(keepprod);
agg_q_prod = q_prod(keepprod);
agg_IndTax_prod  = IndTax(keepprod);
agg_Consumption_prod = Consumption(keepprod);
agg_Investment_prod = Investment(keepprod);
agg_FinalDemand_prod = FinalDemand(keepprod);
%
agg_Supply = zeros(indnum,Nprod);
agg_Uses = zeros(Nprod,indnum);
agg_VA = zeros(size(VA,1),indnum);
indzeros = zeros(1,indnum);
agg_q_j = indzeros;
agg_VAInf_j = indzeros;
agg_LInf_j = indzeros;
agg_ExitRate_j = indzeros;
agg_shInfVA_j = indzeros;
agg_shInfL_j = indzeros;
for j=1:indnum
    agg_Supply(j,:)     = sum(Supply(indcell{j},keepprod));
    agg_Uses(:,j)       = sum(Uses(keepprod,indcell{j}),2);
    agg_VA(:,j)         = sum(VA(:,indcell{j}),2);
    agg_q_j(1,j)        = sum(q_sec(indcell{j}));
    agg_VAInf_j(1,j)    = sum(VAInf_j(indcell{j}));
    agg_LInf_j(1,j)     = sum(LInf_j(indcell{j}));
    agg_ExitRate_j(1,j) = sum(q_sec(indcell{j}).*ExitRate_j(indcell{j}))/agg_q_j(1,j); % Mean exit rates, weighted by gross output
    agg_shInfVA_j(1,j)  = agg_VAInf_j(1,j)/sum(VA_sec(indcell{j}));
    agg_shInfL_j(1,j)   = agg_LInf_j(1,j)/sum(L_sec(indcell{j}));
end
agg_B  = agg_Uses*diag(1./agg_q_j);         % Technical coefficients matrix
agg_D0 = agg_Supply*diag(1./agg_q_prod);    
agg_D  = agg_D0*diag(1./sum(agg_D0));       % Market share matrix: by product, share produced in each industry
assert(max(abs(sum(agg_D)-1))<1e-12)
agg_DD = agg_Supply'*diag(1./agg_q_j);      % Product share matrix: by industry, share of each product in gross output
assert(max(abs(sum(agg_DD)-1))<1e-12)
%
agg_C_j = agg_D*agg_Consumption_prod;
agg_Inv_j = agg_D*agg_Investment_prod;
agg_FinalDemand_j = agg_D*agg_FinalDemand_prod;
%
agg_VA_j = agg_VA(1,:);             % value added by industry
agg_GO_j = agg_VA(end-1,:);         % gross output by industry
agg_L_j = agg_VA(end,:);            % ocupations ("employment") by industry
agg_totwages_j  = agg_VA(2,:);      % labor compensation (wages + payroll taxes + employer private pension contributions) by industry
agg_wages_j = agg_VA(3,:);          % wages by industry
agg_taxwages_j = agg_VA(4,:);       % payroll taxes + employer private pension contributions by industry
agg_GOS_j = agg_VA(end-4,:);        % gross operational surplus ("profits") by industry
%
agg_WVA = sum(agg_totwages_j)./sum(agg_totwages_j+agg_GOS_j);
agg_WVA_j = agg_totwages_j./(agg_totwages_j+agg_GOS_j);
agg_wVA = sum(agg_wages_j)./sum(agg_wages_j+agg_GOS_j);
agg_wVA_j = agg_wages_j./(agg_wages_j+agg_GOS_j);
agg_taxWw0 = sum(agg_taxwages_j)./sum(agg_wages_j);
agg_taxWw0_j = agg_taxwages_j./agg_wages_j;
agg_VAjLj = 1000*agg_VA_j./agg_L_j;
agg_GDPL = 1000*sum(agg_VA_j)./sum(agg_L_j);
agg_FormVAjLj = 1000*(agg_VA_j-agg_VAInf_j)./(agg_L_j-agg_LInf_j);
agg_InfVAjLj = 1000*agg_VAInf_j./agg_LInf_j;
agg_shInfGDP = sum(agg_VAInf_j)/sum(agg_VA_j);
agg_shInfL = sum(agg_LInf_j)/sum(agg_L_j);

%% Compute shares
% Input-output table, D*B
IOT = agg_D*agg_B;
% Value added share of gross output
vash = 1-sum(IOT)';
% Intermediate inputs shares (sum = 1 in the columns)
ioiish = IOT/diag(sum(IOT));
% Consumption expenditure shares
findemsh = agg_C_j/sum(agg_C_j);
% Consumption share on final demand
consh_j = agg_C_j./agg_FinalDemand_j;

% Check totals
assert(sum(sum(ioiish)-1) < 1e-12)
assert(sum(findemsh)-1 < 1e-12)

%% Tax rates

% After reform
%aftertax_prod0 = SO.AfterRefMainTax * aftertaxfactor_prod;
%aftertax_prod = aftertax_prod0(keepprod);
aftertax_prod = aftertaxfactor_prod(keepprod);
aftertax = agg_DD'*aftertax_prod;

% Before reform
agg_taxva_j = (agg_D*agg_IndTax_prod)./agg_VA_j';
agg_taxGDP = sum(agg_IndTax_prod)/sum(agg_VA_j);

%% Data moments
Data.name = 'Data moments';
Data.TaxGDP = agg_taxGDP;
Data.TaxWw = agg_taxWw0 + OtherTaxWw;
Data.TaxWw0 = agg_taxWw0;
Data.OtherTaxWw = OtherTaxWw;
%Data.GDP = sum(agg_VA_j);
Data.L = sum(agg_L_j);
Data.GDP_L = agg_GDPL;  %NEW
Data.W_GDP = agg_WVA;
Data.w_GDP = agg_wVA;
Data.InfGDP = agg_shInfGDP;
Data.InfL = agg_shInfL;
Data.TaxYVA_j = agg_taxva_j';
Data.TaxWw0_j = agg_taxWw0_j;
Data.L_j = agg_L_j;
Data.VAj_Lj = agg_VAjLj;
Data.FormVAjLj = agg_FormVAjLj;
Data.InfVAjLj = agg_InfVAjLj;
Data.WVA_j = agg_WVA_j;
Data.wVA_j = agg_wVA_j;
Data.InfVA_j = agg_shInfVA_j;
Data.InfL_j = agg_shInfL_j;
Data.VAj_GDP = agg_VA_j/sum(agg_VA_j);
Data.ExitRate_j = agg_ExitRate_j;
Data.Investment_GOS = sum(Investment) / sum(agg_GOS_j);
jE = 2;
Data.Consumption_FinalDemand_jE = consh_j(jE);


