%% Load microdata

%%%
nu0e = 6;                       % fit Pareto Type I distribution for firms with 'nu0e' or more employees
%%%

J = SO.NumberofSectors;         % number of sectors
loadexcelmicro = SO.LoadData;   % load Excel microdata [1] or use previously loaded values [0]
FN.codeoutputs = "..\CodeOutputs";                              % Code generated outputs 
FN.microdata = "..\RawData\Microdata";                                      % Raw microdata 
datamicro = fullfile(FN.codeoutputs,"microdata.mat");                       % file microdata in MATLAB format

if loadexcelmicro == 1
    % Formal
    Data.n_empForm_j  = cell(J,1);
    Data.n_firmForm_j = cell(J,1);
    %nameform = "distrib_n_workers_formal.xlsx";
    nameform = "distrib_n_workers_formal3.xlsx";
    filetoread = fullfile(FN.microdata,nameform);
    %Data.n_empForm_j{2}  = xlsread(filetoread,'Industry',"A2:A1751");       % Number of employees, sector 2
    %Data.n_empForm_j{3}  = xlsread(filetoread,'Services',"A2:A1777");       % Number of employees, sector 3
    %Data.n_firmForm_j{2} = xlsread(filetoread,'Industry',"B2:B1751");       % Number of firms by size (employees), sector 2
    %Data.n_firmForm_j{3} = xlsread(filetoread,'Services',"B2:B1777");       % Number of firms by size (employees), sector 3
    Data.n_empForm_j{2}  = xlsread(filetoread,'Industry',"A2:A1254");       % Number of employees, sector 2
    Data.n_firmForm_j{2} = xlsread(filetoread,'Industry',"B2:B1254");       % Number of firms by size (employees), sector 2
    Data.n_empForm_j{3}  = xlsread(filetoread,'Services',"A2:A1790");       % Number of employees, sector 3
    Data.n_firmForm_j{3} = xlsread(filetoread,'Services',"B2:B1790");       % Number of firms by size (employees), sector 3
    % Informal
    %Data.n_empInf_j  = cell(J,1);
    %Data.sh_firmInf_j = cell(J,1);
    %nameinf = "distrib_n_workers_informal.xlsx";
    %filetoread = fullfile(FN.microdata,nameinf);
    %Data.n_empInf_j{2}  = (1:1:5)';                                         % Number of employees, sector 2
    %Data.n_empInf_j{3}  = (1:1:5)';                                         % Number of employees, sector 2
    %Data.sh_firmInf_j{2} = xlsread(filetoread,'Industry',"C2:C6");          % Share of firms by size (employees), sector 2
    %Data.sh_firmInf_j{3} = xlsread(filetoread,'Services',"C2:C6");          % Share of firms by size (employees), sector 3
    %
    % ADICIONAR: SIMPLES, Agricultura
    
    %
    save(datamicro,'Data')
else
    load(datamicro)
end


%% Compute employment distributions

% Number of formal employees by sector
empForm_j = zeros(1,J);
for j=2:J
    empForm_j(j) = sum(Data.n_empForm_j{j}.*Data.n_firmForm_j{j});
end
Data.empForm_j = empForm_j;

% Empirical employment grids
Data.emp_max_j = [0,Data.n_empForm_j{2}(end),Data.n_empForm_j{3}(end)];
emp_max = max(Data.emp_max_j);
emp_max_grid = floor(emp_max/1e3 +1)*1e3;
emp_grid_Form = (1:1:emp_max_grid)';      % employment grid
Ne = length(emp_grid_Form);

% Frequencies
%for each firm size (number of employees) in=1:Ne, 'empForm_num_j' counts the number of formal firms with the correspondig size in the data
empForm_num_j = zeros(Ne,J); % number of employees (1 to Ne) x sectors
for j=2:J
    for in = 1:length(Data.n_empForm_j{j})
        ne = Data.n_empForm_j{j}(in);
        empForm_num_j(ne,j) = Data.n_firmForm_j{j}(in);
    end
end
Data.Nf_j = sum(empForm_num_j);         % Number of firms by industry
Data.Nf_j(1) = 1;
empForm_freq_j = empForm_num_j./Data.Nf_j;
assert(abs(sum(sum(empForm_freq_j))-2) < 1e-10)
%
Data.empForm_num_j = empForm_num_j;     % Number of firms by firm size (1 to Ne employees) and sector j={2,3}
Data.empForm_freq_j = empForm_freq_j;   % Frequency of firms by firm size (1 to Ne employees) and sector j={2,3}

% cdf
empForm_cdf_j = zeros(size(empForm_freq_j));
for j=2:J
    empForm_cdf_j(1,j) = empForm_freq_j(1,j);
    for ie = 2:size(empForm_freq_j,1)
        empForm_cdf_j(ie,j) = empForm_cdf_j(ie-1,j) + empForm_freq_j(ie,j);
    end
end
Data.empForm_cdf_j = empForm_cdf_j;     % CDF of firm size distribution

% Formal firms population 
emp_f_j = cell(1,J);
for j=2:J
    i0=1;
    emp_f_j{j} = zeros(Data.Nf_j(j),1);
    for ie=1:Ne
        nfirms = empForm_num_j(ie,j);
        if nfirms > 0
            emp_f_j{j}(i0:i0-1+nfirms) = ie*ones(nfirms,1);
            i0 = i0+nfirms;
        end
    end
end
Data.empForm_f_j = emp_f_j;     % Formal firms population, RAIS dataset

% Formal employment empirical distribution
 % Fit a Generalized Pareto Distribution ('gp') in the data, for number of employees >= 6.
 % 'gp' needs 3 parameters (k,sigma,theta), we use a Pareto distribution proxy with 2 parameters,
 % which is a special case of 'gp' for k > 0 and sigma = k*theta. 

xi0e_j = zeros(1,J);
%
theta_gp=nu0e-0.0001;   % 'fitdist(...,'gp') requires observed values strictly greater than 'theta'
X=nu0e:Ne;
for j=2:J
    % Maximum likelihood estimator, Pareto Type I distribution    
    empForm_f = emp_f_j{j}(emp_f_j{j}>theta_gp);
    N = length(empForm_f);
    grid = min(empForm_f):1:max(empForm_f);
    xi0e_j(j) = N/sum(log(empForm_f/nu0e));
end

% Keep data
Data.Ne = Ne;                   % Maximum employment grid value
Data.nu0e = nu0e;
Data.xi0e_j = xi0e_j;           % Formal employment firm size pareto distribution parameter by sector j, number of employees >= nu0e

