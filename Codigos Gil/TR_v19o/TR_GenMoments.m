%//////////////////////////////////////////////////////////////////////////
% Indirect Tax Reform
%//////////////////////////////////////////////////////////////////////////
% This function computes the vector of moments
%//////////////////////////////////////////////////////////////////////////

% OBS: ACRESCENTAR MOMENTOS REFERENTES A SIMPLES E CARGA TRIBUTÁRIA

function Msim = TR_GenMoments(SO,Totals,jo)

MuH = Totals.MuH;
%MuX = Totals.MuX;
%MuS = Totals.MuS;
%MuF = MuH + MuX + MuS;
MuF = MuH;
LH = Totals.LH;
%LX1 = Totals.LX1;
%LS1 = Totals.LS1;
%LF1 = Totals.LF1;
Lh = Totals.Lh;
%Lx1 = Totals.Lx1;
%Ls1 = Totals.Ls1;
%Lf1 = Totals.Lf1;
IdFormH = Totals.IdFormH_f;
%IdFormX = Totals.IdFormX_f;
%IdFormS = Totals.IdFormS_f;
%IdForm = IdFormH + IdFormX + IdFormS;
IdForm = IdFormH;
if SO.FirmOptions > 1
    % Informal firms
    MuI = Totals.MuI;
    LI  = Totals.LI;
    Li = Totals.Li;
    IdInf = Totals.IdInf_f;
end


%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Share of informal workers. 
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


LH(imag(LH)~=0) = nan;
%LX(imag(LX)~=0) = nan;
%LS(imag(LS)~=0) = nan;

LH(isnan(LH))=0;
%LX(isnan(LX))=0;
%LS(isnan(LS))=0;

if SO.FirmOptions > 1
    LI(imag(LI)~=0) = nan;
    LI(isnan(LI))=0;
end

LF = LH;
%LF1 = LH1 + LX1 + LS1;
if SO.FirmOptions == 1
    L = LF;
    infwork = 0;
elseif SO.FirmOptions == 2
    L = LI + LF;
    infwork = LI / L;
end


%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Share of informal firms:
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% Size variable:

%%INF%% Li(isnan(Li))=0;
Lh(isnan(Lh))=0;
%Lx1(isnan(Lx1))=0;
%Ls1(isnan(Ls1))=0;

%Lf = Lh1 + Lh2 + Lx1 + Lx2 + Ls1 + Ls2;
Lf = Lh;
Lf = round(Lf);

if SO.FirmOptions == 2  % Formal + informal
    Li = round(Li);
    ShareFirmsInf_all = MuI/(MuI + MuF);

    %-----------------------------
    % Firms with <=2 employees
    %-----------------------------
    
    % Share of informal
    y1=(Li<=2 & IdInf==1);
    x0=sum(y1)/sum(IdInf);
    
    % Share of formal 
    y2=(Lf<=2 & IdForm==1);
    x1=sum(y2)/sum(IdForm); 
    
    ShareFirmsInf_0to2emp = (x0*MuI)/(x0*MuI + x1*MuF);
    
    clear y1 y2 x0 x1
    %-----------------------------
    % Firms with 3-4 employees
    %-----------------------------
    
    % Share of informal 
    y1=(Li>=3 & Li<=4 & IdInf==1);
    x0=sum(y1)/sum(IdInf); 
    
    % Share of formal 
    y2=(Lf>=3 & Lf<=4 & IdForm==1);
    x1=sum(y2)/sum(IdForm); 
    
    ShareFirmsInf_3to4emp = (x0*MuI)/(x0*MuI + x1*MuF);
    clear y1 y2 x0 x1
    
    %-----------------------------
    % Firms with 5-10 employees
    %-----------------------------
    
    % Share of informal 
    y1=(Li>=5 & Li<=10 & IdInf==1);
    x0=sum(y1)/sum(IdInf); 
    
    % Share of formal 
    y2=(Lf>=5 & Lf<=10 & IdForm==1);
    x1=sum(y2)/sum(IdForm); 
    
    ShareFirmsInf_5to10emp = (x0*MuI)/(x0*MuI + x1*MuF);
    clear y1 y2 x0 x1

    %-----------------------------
    % Firms with 10+ employees
    %-----------------------------
    
    % Share of informal 
    y1=(Li>10 & IdInf==1);
    x0=sum(y1)/sum(IdInf); 
    
    % Share of formal 
    y2=(Lf>10 & IdForm==1);
    x1=sum(y2)/sum(IdForm); 
    
    ShareFirmsInf_10moreemp = (x0*MuI)/(x0*MuI + x1*MuF);
    clear y1 y2 x0 x1
    
    
    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    % Firm size distribution: Informal sector
    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    y=(Li<=2 & IdInf==1);
    LiDist_0to2emp = sum(y)/sum(IdInf);
    clear y
    
    y=(Lf>2 & Lf<=5  & IdInf==1);
    LiDist_0to5emp = sum(y)/sum(IdInf);
    clear y
end

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Firm size distribution: formal sector
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

Lf(isnan(Lf))=0;

y=(Lf<=5 & IdForm==1);
LfDist_0to5emp = sum(y)/sum(IdForm);
clear y

y=(Lf>5 & Lf<=10 & IdForm==1);
LfDist_5to10emp = sum(y)/sum(IdForm);
clear y

y=(Lf>10 & Lf<=20 & IdForm==1);
LfDist_10to20emp = sum(y)/sum(IdForm);
clear y

y=(Lf>20 & Lf<=50 & IdForm==1);
LfDist_20to50emp = sum(y)/sum(IdForm);
clear y

y=(Lf>50 & IdForm==1);
LfDist_50moreemp = sum(y)/sum(IdForm);
clear y

Msim.name = "Simulated moments";
Msim.LfDist_0to5emp = LfDist_0to5emp;
Msim.LfDist_5to10emp = LfDist_5to10emp;
Msim.LfDist_10to20emp = LfDist_10to20emp;
Msim.LfDist_20to50emp = LfDist_20to50emp;
Msim.LfDist_50moreemp = LfDist_50moreemp;
if SO.FirmOptions == 2  % Formal + informal
    Msim.infwork = infwork;
    Msim.ShareFirmsInf_all = ShareFirmsInf_all;
    Msim.ShareFirmsInf_0to2emp = ShareFirmsInf_0to2emp;
    Msim.ShareFirmsInf_3to4emp = ShareFirmsInf_3to4emp; 
    Msim.ShareFirmsInf_5to10emp = ShareFirmsInf_5to10emp;
    Msim.ShareFirmsInf_10moreemp = ShareFirmsInf_10moreemp;
    Msim.LiDist_0to2emp = LiDist_0to2emp;
    Msim.LiDist_3to5emp = LiDist_0to5emp;
end

%% Figures for moments check
if SO.plotcheckfig == 1
    
    % RAIS, share of formal employment by firm size and industry
    DData.LfDist_0tolessthan3emp    = [0.321792387	0.46580732	0.441614746];
    DData.LfDist_3to5emp            = [0.212460797	0.205325414	0.247013911];
    DData.LfDist_6to10emp           = [0.166725033	0.13137197	0.155498379];
    DData.LfDist_11to20emp          = [0.130343427	0.091124816	0.086376368];
    DData.LfDist_21to50emp          = [0.09633538	0.06473481	0.046855152];
    DData.LfDist_morethan50         = [0.072342975	0.041635671	0.022641443];
    J=3
    j0=2
    DData.cdf = zeros(5,2,J-j0+1);
    for j = j0:J
        %
        k = 0;
        k=k+1; DData.cdf(k,1,j) = 2; 
        k=k+1; DData.cdf(k,1,j) = 5; 
        k=k+1; DData.cdf(k,1,j) = 10; 
        k=k+1; DData.cdf(k,1,j) = 20; 
        k=k+1; DData.cdf(k,1,j) = 50; 
        %
        k = 0;
        k=k+1; DData.cdf(k,2,j) = DData.LfDist_0tolessthan3emp(j); 
        k=k+1; DData.cdf(k,2,j) = DData.LfDist_3to5emp(j) + DData.cdf(k-1,2,j);
        k=k+1; DData.cdf(k,2,j) = DData.LfDist_6to10emp(j) + DData.cdf(k-1,2,j);
        k=k+1; DData.cdf(k,2,j) = DData.LfDist_11to20emp(j) + DData.cdf(k-1,2,j);
        k=k+1; DData.cdf(k,2,j) = DData.LfDist_21to50emp(j) + DData.cdf(k-1,2,j);
        %
        assert(abs(DData.cdf(k,2,j) + DData.LfDist_morethan50(j) - 1)<1e-5)
    end
    
        




    
    y0 = IdForm.*Lf;
    empo_f = y0(y0>0);
    emp_f = sort(empo_f);
    N_f = length(emp_f);
    num_f = (1:N_f)'/N_f;

    sample1 = log(empo_f);
    
%    % cdf
%    pd0 = fitdist(log(empo_f), 'Normal');
%    sig0 = pd0.sigma;
%    mu0 = pd0.mu;
%    cdf0 = normcdf(log(emp_f),mu0,sig0);
%    figure;
%    [fo,x1] = ecdf(sample1); plot(x1,fo); %CDF
%    hold on
%    plot(log(emp_f),cdf0)
%    hold on    
%    j=jo;
%    scatter(log(DData.cdf(:,1,j)),DData.cdf(:,2,j))
    %grid
%    xlabel('Log of number of workers')
%    ylabel('CDF')
%    legend({'Model simulation','Normal approximation','Data'},'Location','southeast')
%    hold off
    
   piecesplot = 0;
   if piecesplot == 1
    % cdf PEDAÇOS
    pd0 = fitdist(log(empo_f), 'Normal');
    sig0 = pd0.sigma;
    mu0 = pd0.mu;
    cdf0 = normcdf(log(emp_f),mu0,sig0);
    figure;
    [fo,x1] = ecdf(sample1); plot(x1,fo); %CDF
    xlabel('Log of number of workers')
    ylabel('CDF')
    %legend({'Model simulation','Normal approximation','Data'},'Location','southeast')
    legend({'Model simulation'},'Location','southeast')
    hold off
    
    % cdf
    pd0 = fitdist(log(empo_f), 'Normal');
    sig0 = pd0.sigma;
    mu0 = pd0.mu;
    cdf0 = normcdf(log(emp_f),mu0,sig0);
    figure;
    [fo,x1] = ecdf(sample1); plot(x1,fo); %CDF
    hold on
    plot(log(emp_f),cdf0)
    hold on    
    j=jo;
    %scatter(log(DData.cdf(:,1,j)),DData.cdf(:,2,j))
    %grid
    xlabel('Log of number of workers')
    ylabel('CDF')
    legend({'Model simulation','Normal approximation'},'Location','southeast')
    hold off
   end
   
    jumpdf=1
    if jumpdf==0
        
        % pdf
        %[p,x0] = hist(sample1); figure; plot(x0,p/sum(p)); %PDF
        figure;
        emp_pdf = histogram(sample1);
    
        pdf0 = normpdf(log(emp_f),mu0,sig0);
        figure
        plot(log(emp_f),pdf0/sum(pdf0))
        hold on
        
        n0 = emp_pdf.BinEdges(1:end-1);
        n1 = emp_pdf.BinEdges(2:end);
        nn=mean([n0;n1]); 
        pdo = emp_pdf.Values/sum(emp_pdf.Values);
        plot(nn,pdo)
        
        figure
    
    
        f = figure;
%       in = in+1;subplot(nl,nc,in,'Parent',figu)
        scatter(log(emp_f),num_f)
        xlabel('Log of number of workers');ylabel('Cumulative distribution')
        title('Formal firms')
    
        hold on
        pd0 = fitdist(log(empo_f), 'Normal');
        sig0 = pd0.sigma;
        mu0 = pd0.mu;
        cdf0 = normcdf(log(emp_f),mu0,sig0);
        plot(log(emp_f),cdf0)
        
        xo = (1:N_f)';
        yo = pdf(pd0,xo);
        figure; plot(xo, yo)
        
        [p,x0] = hist(sample1); figure; plot(x0,p/sum(p)); %PDF
        [f,x1] = ecdf(sample1); figure; plot(x1,f); %CDF
    
        % empo_f = y0(y0>0);
        % lpd = fitdist(log(empo_f), 'Normal')
        
        loo=1
        
        f = figure;
%        in = in+1;subplot(nl,nc,in,'Parent',figu)
        scatter(log(emp_f),num_f*N_f)
        xlabel('Log of number of workers');ylabel('Cumulative distribution')
        title('Formal firms')
        
    
    
        hold on
        sig = 1.1
        muu = 3.7
        cdf0 = normcdf(log(emp_f),muu,sig)
        scatter(log(emp_f),cdf0)
        
        a=0
    end
    
    
end
end

