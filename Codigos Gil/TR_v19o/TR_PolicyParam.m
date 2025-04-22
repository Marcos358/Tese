%%% Policy parameters
% Payroll taxes

if is == 1     % Brazil before 2023 tax reform
    %PolicyPar.tauY_j = (1+ Data.InfVA_j./(1-Data.InfVA_j)).* Data.TaxYVA_j;
    PolicyPar.tauY_j = TaxY_j;
elseif is == 2   % Brazil after 2023 tax reform
    PolicyPar.tauY_j = SO.AfterRefMainTax*aftertax';                           % Value added tax rates, EC 132, code 'TR_LoadIOTSUT.m'  
    PolicyPar.tauY_j(1) = 0.5*PolicyPar.tauY_j(1);
    if SO.lowercompliance == 1
        %FixPar.EForm_j = FixPar.ESimp_j;
        FixPar.EForm_j = FixPar.EForm_j/2;
    end
end


%% SIMPLES parameters

%%% WAGE
% From IBGE [VER CITAÇÂO], mean monthly wage in 2015 is R$ 2480
% Informality rate in 2015 [VER CITAÇÃO] = 45%
% Monthly to yearly wage: Gomes, Iachan, and Santos (2020) multiply
% informal wages by 12 and formal wages by 13.33 (additional thirteenth
% monthly salary every year plus one-third of a monthly salary as vacation allowance)
% then: w = 2480*(12*0.45 + 13.33*0.55) = 31574

%%% Revenue limit in SIMPLES in 2015: 3600 thousands of Reais
% then: R/w = 3600 / 31.574 = 114

%%% SIMPLES tax rates
%From Matsumoto (2021):
% Figure 2.3, average tax rate in 2012-2016,
 % firms with revenues below SIMPLES revenue limit
 % 19 to 17.5 in formal firms outside SIMPLES
 % 11 to 12.5 in SIMPLES firms
 % then: tauYSimp/tauY = mean([11/19,12.5/17.5]) = 0.65

%%% Compliance costs
%From Matsumoto (2021):
% Figure 2.4, average compliance cost rate in 2012-2016, share of Gross Revenue
 % firms with revenues below SIMPLES revenue limit
 % 4.5 to 3 in formal firms outside SIMPLES
 % 2 to 1.5 in SIMPLES firms
 % then: ESimp/EForm = mean([2/4.5,1.5/3]) = 0.4722, circa 0.5
 
 