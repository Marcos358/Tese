%//////////////////////////////////////////////////////////////////////////
% Indirect Tax Reform
%//////////////////////////////////////////////////////////////////////////
% Main code
% compute all the results in the paper
% before runing this code, chose simulation options in 'aSimOptions.m'
%//////////////////////////////////////////////////////////////////////////

clear;
close all;
%clc;
tic;

% Assign directories 
FN.name = "Folder names for saving and loading files";
FN.dataIOT = "..\RawData\IOTSUT";                               % Input-output tables data
FN.codeoutputs = "..\CodeOutputs";                              % Code generated outputs 
FN.modelinputs = fullfile(FN.codeoutputs,"modelinputs.mat");    % File with model inputs
FN.firmsimuls  = fullfile(FN.codeoutputs,"Estimates.mat");      % File with firm simulations
%FN.firmsimuls  = fullfile(FN.codeoutputs,"Estimates2.mat");    % File with firm simulations

% Chose simulation options
aSimOptions

% Firms random simulations
SimulateFirms = 0;
if SimulateFirms == 1
    RandomMSM
    save(FN.firmsimuls,'X1', 'X2', 'N') 
end

% Load data
TR_LoadIOTSUT

% Load micro data
TR_LoadMicrodata
%SO.SetEndoParams=1;

% Saving data
Data
save(FN.modelinputs,"FN","SO","vash","ioiish","findemsh","aftertax","Data") 
clearvars -except FN 
modelpar = fullfile(FN.codeoutputs,"modelparameters.mat");

% Parametrization 
load(FN.modelinputs)
if SO.SetEndoParams > 0         % Model parametrization
    %%% Prices and guesses
    %tY_j = (1./(1-Data.InfVA_j)).* Data.TaxYVA_j;     % Standard tax rate
    %tY_j(1) = Data.TaxYVA_j(1);

    % Prices
    w = 1;                              %[FIXO] Initial steady state wage 
    p_j = [1,1,1];%                     %[FIXO] Initial steady state prices [p_j(1) = 1]
    
    %%% %%% %%%
    % Single parameters
    A = 950;%871;%2283;%937;                         % Total factor productivity in agriculture
    tauS_tauY2 = 0.4118;%0.4; % SIMPLES Y tax rate / Formal standard tax rate in sector j=2;
    % All firms
    EForm_j = [1,9,6];%[1,7,7];%0.985[1,4.2,3.97];%[1,4.18,3.86];%[1,1,1]*4;        % Entry cost, formal firms standard taxation
    xi_j = [1,4.8,4.33941];%[1,4.6659,4.3941];%[1,4,3];%[1,4.58,2.47];%[1,4.46,2.62];%[1,1,1]*4;                % Pre-entry productivity, pareto distribution shape parameter  {TAMANHO DO SIMPLES -> reduz xi, reduz SIMPLES (deve aumentar size das firmas)]
    sigmaz_j = [0.245,0.17,0.5334];%[0.245,0.238,0.245];%0.245*[1,1,1];%[1,0.15,0.5];%[1,0.2,0.5];%[1,0.45,0.45];%[1,0.4056,0.4671];%           % Post-entry productivity, standard deviation of the unantecipated shock
    %%% %%% %%%
    % Informal firms
    kappaInfForm_j = [1,2,1];%[1,1,1];%[1,4.5,2.8];%[1,1.9,1.9];%[1,1.5,1.37];%[1,2.8,1.45];%         % Discount rates, informal to standard formal ratio, by sector. Ulyssea (2018), informal death rate circa 3x formal death rate
    EInfForm_j = [1,1,1]*0.47;%[1,0.477,0.448];%[1,0.477,0.448];%[1,1,1]*0.47;%[1,0.254,0.38];%[0.31,0.2613,0.3851];%        % Entry cost, informal to standard formal ratio, by sector. Einf/Eform = 0.47 in Ulyssea (2018)
    % Taxes
    %TaxY_j = (1./(1-Data.InfVA_j)).* Data.TaxYVA_j;    % Standard tax rate
    TaxY_j = [1,0.4872,0.1435];%0.1435;%[1,0.4935,0.1173];
    TaxY_j(1) = Data.TaxYVA_j(1);    
    % SIMPLES firms
    kappaSimpForm_j = [1,1,1];%[1,1,1];%[1,2.5,1.65];%[1,1.1765,1.0813];% [FIXED] Discount rates, SIMPLES to standard formal ratio, by sector. Ulyssea (2018), informal death rate circa 3x formal death rate
    ESimpForm_j = [1,1,1]*0.5;%[1,0.263,0.417];%[1,0.2835,0.2172];%       % Entry cost, SIMPLES to standard formal ratio, by sector, discussion in 'TR_PolicyParam'
    RmaxSimp_j = [1,114,114];%110.8102];%114*w*[1,1,1];               % Maximum revenue allowed for SIMPLES tax regime [discussion in PolicyPar]
    % Tax rates
    %TaxY_j = [tY_j(1),0.5152,0.1140];%[tY_j(1),0.4362,0.1139];%[tY_j(1),0.4290,0.1719];%[6.2,25.9,25]/100;       % Stantard value added tax rate [fixed for j=1, initial guess = tY_j]
    %tauYSimpForm_j = [1,0.4,0.7336];%[1,1,1]*0.65;%% Value added tax rates [discussion below]
    
    %%%
    
    
    % Remaining parameters (except policy parameters)
    TR_SetParam
    
else                            % Policy scenarios
    load(modelpar);             % Load parameters 
    %FixPar.zeta_j = Data.VAj_GDP;
    
    %TR_Priceguess               % Initial guess for prices
    Prices.w = w;
    Prices.p_j = p_j;
    Prices.pE = Prices.p_j(jE);
    Prices = TR_Piota(Prices,FixPar);                       % piota_j, intermediate inputs price index

end

SO.checkfig = 0;
    
    if SO.SetEndoParams > 0
        is = 1;
        TR_PolicyParam
        if SO.SetEndoParams == 1                % COLOCAR AQUI O ALGORITMO DE SMMM
            TR_SMM_Main
        end
        TR_PartialEq
        %
        clear SO
        save(modelpar)
    
    else
        for is = SO.SimScenarios
            IS = is
            TR_PolicyParam
            % Compute equilibrium
            %%%
            shnew = 2/3;%1/4;%2/3;%3/4;
            toler = 0.001;
            jmm = 2;        % index for the manufacturing industry 
            %%%
            Prices0 = Prices;
            if is > 2
                %Prices0.w = 1.5*Prices0.w;
                %Prices0.p_j(3) = 1.5*Prices0.p_j(3);
                %Prices0.p_j(1) = 1.3*Prices0.p_j(1);
            end
            [Prices,Totals_j,Msim_j,Gov,Hous,walras] = TR_FindEquilibrium(SO,FixPar,PolicyPar,EndoPar,Prices0,shnew,toler,jmm,is);    
    
        % Keep results
        Results.TauY{is} = PolicyPar.tauY_j;
        Results.Prices{is} = Prices;
        Results.Totals{is} = Totals_j;
        Results.Msim{is} = Msim_j;
        Results.Hous{is} = Hous;
        Results.Gov{is} = Gov;
        if SO.FirmOptions > 2       % add SIMPLES firms
            Results.TauYSimp{is} = PolicyPar.tauYSimp_j;
            Results.TaxYSimp_TaxY_j{is} = Gov.TaxYSimp_TaxY_j;
        end        
        end
        
        % Tax reform results
        TR_TabResults
        
    end

toc

