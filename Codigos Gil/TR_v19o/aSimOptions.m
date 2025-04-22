%//////////////////////////////////////////////////////////////////////////
% Indirect Tax Reform
%//////////////////////////////////////////////////////////////////////////
% aSimOptions
% chose simulation options (SO)
%//////////////////////////////////////////////////////////////////////////

SO.name = "Simulation options in 'aSimOptions.m'";

%% Calibration or policy scenarios
% [2] Alternative model parametrization {na mão} 
% [1] SMM estimation
% [0] Run policy scenarios

SO.SetEndoParams = 0

%% After reform: main tax rate
SO.AfterRefMainTax = 0.315;%0.315;

%% After reform: reduction in compliance costs

% [0] Do not reduce compliance costs
% [1] Reduce compliance costs
SO.lowercompliance = 1;

%% Firm status, alternative models
% [1] Only formal firms
% [2] Formal and informal firms
% [3] Formal, informal, and SIMPLES firms

% Chose firm status options (scalar)
SO.FirmOptions = 1; 

%if SO.FirmOptions == 1
%    SO.AfterRefMainTax = 0.32;%0.28;
%elseif SO.FirmOptions == 2
%    SO.AfterRefMainTax = 0.32;
%end    

%% Load data
% [0] Use previously loaded data
% [1] Load data

% Load data or not (scalar)
if SO.SetEndoParams == 0
    SO.LoadData = 0; 
else
    SO.LoadData = 1;
end
    
%% Tax Reform Scenarios
% [1] Before the tax reform - baseline for parametrization
% [2] 2023 tax reform

% List of simulated scenarios (vector)
if SO.SetEndoParams > 0
    SO.SimScenarios = 1;
else
    SO.SimScenarios = [1,2];
end
%% Preliminary simulations
% [0] Use previously simulated firm simulations
% [1] Simulate firm random shocks

% Run preliminary simulations or not (scalar)
SO.PreSimul = 0; 

%% Number of products

% Chose number of sectors (scalar)
SO.NumberofSectors = 3; 

%% Figures for checking model pamaretrization
% [0] do not plot check figures
% [1] plot check figures

% Plot check figures or not (scalar)
SO.plotcheckfig = 0; 

%% Initial guess for estimation, versions

% [0] Standard guess
% [1] Trying to match some moments TABELA LATEX
SO.initguess = 0; 

%% Productivity known before entry

% [0] Productivity unknown before entry
% [1] Productivity known before entry
SO.prodknown = 0; 

%% Save options
SO
save(FN.modelinputs);

