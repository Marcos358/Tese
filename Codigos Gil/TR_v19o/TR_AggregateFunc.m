%//////////////////////////////////////////////////////////////////////////
% Indirect Tax Reform
%//////////////////////////////////////////////////////////////////////////
% TR_AggregateFunc
% Agregation of firms' decisions
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function Totals = TR_AggregateFunc(SO,Firms,lz_f,Nf,EForm,kappaForm,EInf,kappaInf,ESimp,kappaSimp,j)

% _f = by simulated firm
% h = formal firms, domestic production
% i = informal firms
% s = SIMPLES firms

if j == 1
    LH  = Firms.Eh_f;           % Labor demand, formal
    VAh = Firms.VAh_f;          % Nominal value added
    Yh  = Firms.Yh_f;           % Real value added
    Ph_f = Firms.Ph_f;          % Profits
    Qh  = Firms.qh_f;           % Gross output
    Iotah_k = Firms.iotah_kf;   % Intermediate inputs
    %
    Totals.name = "Firm choices and aggregated number of firms and employment";
    %Totals.MuH = MuH;
    %Totals.NumEntryH = NumEntryH;
    Totals.LH = LH;
    %Totals.Lh = Lh;
    %Totals.IdFormH_f = IdFormH_f;
    Totals.VAh = VAh;
    Totals.Yh = Yh;
    Totals.VAForm = VAh;
    Totals.Profits_Form = Ph_f;
    %Totals.Ecost_Form = Ecost_Form;
    %   
    Totals.L = LH;
    Totals.Q = Qh;
    Totals.VA = Totals.VAForm;
    Totals.Iota_k = Iotah_k;

else

z_f = exp(lz_f);                % Productivity

Vf0_f = Firms.Vf0_f;            % Before entry value function, formal
Eh_f  = Firms.Eh_f;             % Labor demand, formal
VAh_f = Firms.VAh_f;            % Nominal value added
Yh_f  = Firms.Yh_f;             % Real value added
Vh_f  = Firms.Vh_f;             % After entry value function
qh_f  = Firms.qh_f;             % Gross output
iotah_kf = Firms.iotah_kf;      % Intermediate inputs

if SO.FirmOptions > 1
    % Informal firms
    Vi0_f = Firms.Vi0_f;
    Ei_f  = Firms.Ei_f;
    VAi_f = Firms.VAi_f;
    Yi_f  = Firms.Yi_f;
    Vi_f  = Firms.Vi_f;
    qi_f  = Firms.qi_f;
    iotai_kf = Firms.iotai_kf;
end

if SO.FirmOptions > 2
    % SIMPLES firms
    Vs0_f = Firms.Vs0_f;
    Es_f  = Firms.Es_f;
    VAs_f = Firms.VAs_f;
    Ys_f  = Firms.Ys_f;
    Vs_f  = Firms.Vs_f;
    qs_f  = Firms.qs_f;
    iotas_kf = Firms.iotas_kf;
end
%figure;scatter(z_f,log(Vf0_f));hold on;scatter(z_f,log(Vi0_f));hold on;scatter(z_f,log(Vs0_f))

%--------------------------------------------------------------------------
% The index functions:
%--------------------------------------------------------------------------

zeros_f = zeros(Nf,1);
ones_f = ones(Nf,1);
yf0 = Vf0_f - EForm*ones_f;
if SO.FirmOptions > 1
    yi0 = Vi0_f - EInf*ones_f;
    if SO.FirmOptions > 2
        ys0 = Vs0_f - ESimp*ones_f;
    end
end
%figure;scatter(z_f,log(yf0+5));hold on;scatter(z_f,log(yi0+5));hold on;scatter(z_f,log(ys0+5))
%figure;scatter(z_f,log(diff_fi+200));hold on;scatter(z_f,log(diff_hs+200))

yy = zeros_f;
if SO.FirmOptions > 1
    yy = [yy,yi0];
    if SO.FirmOptions > 2
        yy = [yy,ys0];
    end
end
yy = [yy,yf0];

[TProfit,Status] = max(yy,[],2);
Status = Status - 1;

IdFormH_f = zeros_f;
if SO.FirmOptions == 1
% Firm Status: [0] no entry, [1] formal, 
    IdFormH_f(Status==1) = 1; 
elseif SO.FirmOptions == 2
% Firm Status: [0] no entry, [1] informal, [2] formal, 
    IdFormH_f(Status==2) = 1; 
    IdInf_f = zeros_f;
    IdInf_f(Status==1) = 1;    
elseif SO.FirmOptions == 3
% Firm Status: [0] no entry, [1] informal, [2] SIMPLES, [3] formal standard, 
    IdFormH_f(Status==3) = 1; 
    IdInf_f = zeros_f;
    IdInf_f(Status==1) = 1;
    IdSimp_f = zeros_f;
    IdSimp_f(Status==2) = 1;
end

%%% ARRUMAR ou limpar aqui
% Formal wants to export:
%ev = exp(diff_xh/lambda);
%EX = ev./(ones_f + ev);
%EX(isnan(EX)) = 1;
%clear ev

% Actual exporter:
%EntryH_f = EntryHX_f.*(ones_f - EX);
%EntryX_f = EntryHX_f.*EX;


% Number of entrants in each sector:
NumEntryH = sum(IdFormH_f);
if SO.FirmOptions > 1       % add informal firms
    NumEntryI = sum(IdInf_f);
    if SO.FirmOptions > 2   % add SIMPLES firms
        NumEntryS = sum(IdSimp_f);
    end
end

% Mass of firms in each sector:
MuH = NumEntryH/kappaForm;
if SO.FirmOptions > 1   % add informal firms
    MuI = NumEntryI/kappaInf;
    if SO.FirmOptions > 2   % add SIMPLES firms
        MuS = NumEntryS/kappaSimp;
    end
end


%--------------------------------------------------------------------------
% Getting total formal employment and other aggregates
%--------------------------------------------------------------------------

%%% Formal general case

x = MuH/NumEntryH;   % expansion factor
x(isnan(x)) = 1;

[Lh,LH,VAh,Yh,Qh,Iotah_k] = TR_Sectoragg(x,IdFormH_f,Eh_f,VAh_f,Yh_f,qh_f,iotah_kf);

Ecost_Form = EForm * NumEntryH * x; % Entry costs
Profits_Form  = x* sum(IdFormH_f.*Vh_f)* kappaForm  - Ecost_Form; % Total profits, liquid of entry costs


if SO.FirmOptions > 1       % add informal firms
    %--------------------------------------------------------------------------
    % Getting total informal employment and other aggregates
    %--------------------------------------------------------------------------
    
    xi = MuI/NumEntryI;     % expansion factor
    xi(isnan(xi)) = 1;

    [Li,LI,VAi,Yi,Qi,Iotai_k] = TR_Sectoragg(xi,IdInf_f,Ei_f,VAi_f,Yi_f,qi_f,iotai_kf);

    Ecost_Inf = EInf * NumEntryI * xi; % Entry costs
    Profits_Inf  = xi* sum(IdInf_f.* Vi_f)* kappaInf - Ecost_Inf; % Total profits, liquid of entry costs
end


%%% SIMPLES
if SO.FirmOptions > 2       % add SIMPLES firms
    %--------------------------------------------------------------------------
    % Getting total SIMPLES employment and other aggregates
    %--------------------------------------------------------------------------
    
    xs = MuS/NumEntryS;     % expansion factor
    xs(isnan(xs)) = 1;

    [Ls,LS,VAs,Ys,Qs,Iotas_k] = TR_Sectoragg(xs,IdSimp_f,Es_f,VAs_f,Ys_f,qs_f,iotas_kf);

    Ecost_Simp = ESimp * NumEntryS * xs; % Entry costs
    Profits_Simp  = xs* sum(IdSimp_f.* Vs_f)* kappaSimp - Ecost_Simp; % Total profits, liquid of entry costs
end


%%% Totals

Totals.name = "Firm choices and aggregated number of firms and employment";
Totals.MuH = MuH;
Totals.NumEntryH = NumEntryH;
Totals.LH = LH;
Totals.Lh = Lh;
Totals.IdFormH_f = IdFormH_f;
Totals.VAh = VAh;
Totals.Yh = Yh;
Totals.VAForm = VAh;
Totals.Profits_Form = Profits_Form;
Totals.Ecost_Form = Ecost_Form;
%
Totals.L = LH;
Totals.Q = Qh;
Totals.VA = Totals.VAForm;
Totals.Iota_k = Iotah_k;

if SO.FirmOptions > 1               % add informal firms
    Totals.MuI = MuI;
    Totals.NumEntryI = NumEntryI;
    Totals.LI = LI;
    Totals.Li = Li;
    Totals.IdInf_f = IdInf_f;
    Totals.VAi = VAi;
    Totals.Yi = Yi;
    Totals.Profits_Inf = Profits_Inf;
    Totals.Ecost_Inf = Ecost_Inf;
    %
    Totals.L = Totals.L + LI;
    Totals.Q = Totals.Q + Qi;
    Totals.VA = Totals.VA + VAi;
    Totals.Iota_k = Totals.Iota_k + Iotai_k;
    
    if SO.FirmOptions > 2       % add SIMPLES firms
        Totals.MuS = MuS;
        Totals.NumEntryS = NumEntryS;
        Totals.LS = LS;
        Totals.Ls = Ls;
        Totals.IdSimp_f = IdSimp_f;
        Totals.VAs = VAs;
        Totals.Ys = Ys;
        Totals.Profits_Simp = Profits_Simp;
        Totals.Ecost_Simp = Ecost_Simp;
        %
        Totals.L = Totals.L + LS;
        Totals.Q = Totals.Q + Qs;
        Totals.VA = Totals.VA + VAs;
        Totals.Iota_k = Totals.Iota_k + Iotas_k;
    end
end
%Totals.MuF = MuH + MuX + MuS;
%Totals.IdForm_f = IdFormH_f + IdFormX_f + IdFormS_f;
%Totals.VAForm = VAh + VAs + VAx;

% TESTES
figforcheck = 0
if figforcheck == 1
    jj=j
    figure;scatter(z_f,log(Vf0_f));hold on;scatter(z_f,log(Vi0_f));%hold on;scatter(z_f,log(Vs0_f))
    figure;scatter(log(z_f),Status)
    parada=1;
    close all
end


end

end


%%% Firm status
% Entrants: [1] entrant firm, [0] do not enter
%Entrant = zeros_f;
%Entrant(yf0>0) = 1;
%if SO.FirmOptions > 1
%    Entrant(yi0>0) = 1;
%    if SO.FirmOptions > 2
%        Entrant(ys0>0) = 1;
%    end
%end
% Can enter simples
%Status(yi0>yf0) = 1;
%Status(yi0>ys0) = 1;

% StatusFormal: [1] Formal standard > SIMPLES, [0] Formal standard <= SIMPLES
%if SO.FirmOptions > 2
%    StatusFormal = ones_f;
%    StatusFormal(ys0>yf0) = 0;
%end

% Status: [0] no entry, [1] standard formal, [2] informal, [3] SIMPLES
%Status = zeros_f;
%if SO.FirmOptions == 1
%    Status(Entrant>0) = 1;
%else
%    Status(yi0>0) = 2;
%    Status(yf0>yi0) = 1;
%    if SO.FirmOptions > 2
%        Status(ys0>yf0) = 3;
%        if ys0>yf0
%            Status(yi0>ys0) = 3;
%        end
%    end
%    Status(Entrant<1) = 0;
%end

