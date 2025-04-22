%//////////////////////////////////////////////////////////////////////////
% Indirect Tax Reform
%//////////////////////////////////////////////////////////////////////////
% TR_SimEconomy
% Finds Firm Value Functions 
% (and updates labor policy for SIMPLES firms, using the revenue cap)
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [Firms,LaborPolicy] = TR_ValueFunc(SO,wage,p,piota,p_k,lambda_k,alpha,theta,...
            tauY,tauW,kappaForm,Probzz,lz_grid,lz0_f,lz_f,LaborPolicy,j,kappaInf,tauYSimp,tauWSimp,kappaSimp,RmaxSimp)

if j == 1
    
    z_f = exp(lz_grid);
    taxY = tauY;
    taxW = tauW;
    Eh_f = LaborPolicy.ellFormH;
    Emp = Eh_f;
    [VAh_f,Yh_f,Ph_f,Rh_f,qh_f,iotah_kf] = TR_Profits(Emp,z_f,wage,p,piota,p_k,lambda_k,alpha,theta,taxY,taxW);
    
    Vh_f = Ph_f/kappaForm;
    Vf0_f = Vh_f;     
    
else
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Expected before entry profit functions
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

z = exp(lz_grid);
z0_f = exp(lz0_f);

% Formal firms
taxY = tauY;
taxW = tauW;
Emp = LaborPolicy.ellFormH;
%Emp = round(LaborPolicy.ellFormH);
[~,~,ProfitFormH,~,~,~] = TR_Profits(Emp,z,wage,p,piota,p_k,lambda_k,alpha,theta,taxY,taxW);
VForm0=(Probzz*ProfitFormH)/kappaForm;
%Vf0_f = spline(z,VForm0,z0_f);
lVf0_f = spline(z,log(VForm0),z0_f);
Vf0_f = exp(lVf0_f);

if SO.FirmOptions > 1
    % Informal firms
    taxY = 0;
    taxW = 0;
    Emp = LaborPolicy.ellInf;
    %Emp = round(LaborPolicy.ellInf);
    [~,~,ProfitInf,~,~,~] = TR_Profits(Emp,z,wage,p,piota,p_k,lambda_k,alpha,theta,taxY,taxW);
    VInf0=(Probzz*ProfitInf)/kappaInf;
    %Vi0_f = spline(z,VInf0,z0_f);
    lVi0_f = spline(z,log(VInf0),z0_f);
    Vi0_f = exp(lVi0_f);
end
%figure;plot(z,log(LaborPolicy.ellFormH));hold on;plot(z,log(LaborPolicy.ellInf))
%figure;plot(z,log(VForm0));hold on;plot(z,log(VInf0))

if SO.FirmOptions > 2
    % SIMPLES firms
    taxY = tauYSimp;
    taxW = tauWSimp;
    % Revenue cap
    Emp = LaborPolicy.ellSimp;
    %
    [~,~,ProfitSimp,~,~,~,Emp] = TR_ProfitsSimp(Emp,z,wage,p,piota,p_k,lambda_k,alpha,theta,taxY,taxW,RmaxSimp);
    %
    %[~,~,~,~,~,~,Emp] = TR_ProfitsSimp(Emp,z,wage,p,piota,p_k,lambda_k,alpha,theta,taxY,taxW,RmaxSimp);
    LaborPolicy.ellSimp = Emp;
    % Profits
    %Emp = round(LaborPolicy.ellSimp);
    %[~,~,ProfitSimp,~,~,~,~] = TR_ProfitsSimp(Emp,z,wage,p,piota,p_k,lambda_k,alpha,theta,taxY,taxW,RmaxSimp);
    %
    VSimp0=(Probzz*ProfitSimp)/kappaSimp;
    %Vs0_f = spline(z,VSimp0,z0_f);
    lVs0_f = spline(z,log(VSimp0),z0_f);
    Vs0_f = exp(lVs0_f);
    %figure;scatter(z0_f,log(Vf0_f));hold on;scatter(z0_f,log(Vi0_f));hold on;scatter(z0_f,log(Vs0_f))
   
    %%%%%%
    % ARRUMAR, peguei aqui da v6c
    %Prof = [ProfitSimp,ProfitFormH];
    %y = (Prof==max(Prof,[],2));
    %ProfitForm = sum(y.*Prof,2);
    %  ProfitForm = ProfitFormH;
    
    %VForm0=(Probzz*ProfitForm)/kappaForm; 
    
    % Interpolating
    %Vi0_f = spline(z,VInf0,z0_f);
    %Vf0_f = spline(z,VForm0,z0_f);
    %%%%%%

end

%%
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Post-entry firm choices
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

z_f = exp(lz_f);

% Formal firms
taxY = tauY;
taxW = tauW;
Eh_f = exp(spline(lz_grid,log(LaborPolicy.ellFormH),lz_f));
Emp = Eh_f;
%Emp = round(Eh_f);
[VAh_f,Yh_f,Ph_f,Rh_f,qh_f,iotah_kf] = TR_Profits(Emp,z_f,wage,p,piota,p_k,lambda_k,alpha,theta,taxY,taxW);
Vh_f = Ph_f/kappaForm;
%pYh_f = VAh_f ./ Yh_f; % value added price index
%assert(max(pYh_f)-min(pYh_f)<1e-10)

if SO.FirmOptions > 1
    % Informal firms
    taxY = 0;
    taxW = 0;
    Ei_f = exp(spline(lz_grid,log(LaborPolicy.ellInf),lz_f));
    Emp = Ei_f;
    %Emp = round(Ei_f);
    [VAi_f,Yi_f,Pi_f,Ri_f,qi_f,iotai_kf] = TR_Profits(Emp,z_f,wage,p,piota,p_k,lambda_k,alpha,theta,taxY,taxW);
    Vi_f = Pi_f/kappaInf;
    %pYi_f = VAi_f ./ Yi_f; 
    %assert(max(pYi_f)-min(pYi_f)<1e-10)
end

if SO.FirmOptions > 2
    % SIMPLES firms
    taxY = tauYSimp;
    taxW = tauWSimp;
    Es_f = exp(spline(lz_grid,log(LaborPolicy.ellSimp),lz_f));
    Emp = Es_f;
    %Emp = round(Es_f);
    [VAs_f,Ys_f,Ps_f,Rs_f,qs_f,iotas_kf,~] = TR_ProfitsSimp(Emp,z_f,wage,p,piota,p_k,lambda_k,alpha,theta,taxY,taxW,RmaxSimp);
    Vs_f = Ps_f/kappaSimp;
    
    % Revenue cap for SIMPLES firms
    Rcheck = (Rs_f > RmaxSimp + 1e-3);
    assert(sum(Rcheck) < 1)
    
end

end


%% Keep results

% Formal firms
Firms.name = "Firm value functions, employment policies, value added and profits";
Firms.Vf0_f = Vf0_f;
Firms.Eh_f = Eh_f;
Firms.VAh_f = VAh_f;
Firms.Yh_f = Yh_f;
Firms.Ph_f = Ph_f;
Firms.Vh_f = Vh_f;
Firms.qh_f = qh_f;
Firms.iotah_kf = iotah_kf;

if j > 1
if SO.FirmOptions > 1
    % Informal firms
    Firms.Vi0_f = Vi0_f;
    Firms.Ei_f = Ei_f;
    Firms.VAi_f = VAi_f;
    Firms.Yi_f = Yi_f;
    Firms.Pi_f = Pi_f;
    Firms.Vi_f = Vi_f;
    Firms.qi_f = qi_f;
    Firms.iotai_kf = iotai_kf;
end

if SO.FirmOptions > 2
    % SIMPLES firms
    Firms.Vs0_f = Vs0_f; 
    Firms.Es_f = Es_f;
    Firms.VAs_f = VAs_f;
    Firms.Ys_f = Ys_f;
    Firms.Ps_f = Ps_f;
    Firms.Vs_f = Vs_f;
    Firms.qs_f = qs_f;
    Firms.iotas_kf = iotas_kf;
end
end

%% Figures for parametrization check
if j > 1
if SO.plotcheckfig == 1
    
    darkred = [0.6350 0.0780 0.1840];
    cgreen  = [0.2431 0.5882 0.3176];
    
    %% Labor demand
    
    x  = LaborPolicy.ellFormH;
    xi = LaborPolicy.ellInf;
    if SO.FirmOptions > 2
        xs = LaborPolicy.ellSimp;
    end
    
    % Subplot options
    nl = SO.FirmOptions; % number of lines
    nc = 2; % number of columns
    in = 0; % subplot index
    
    f = figure;
    figu = uipanel('Parent',f,'BorderType','none'); 
    figu.Title = ['Labor demand policy function, file "TR_LaborFunc", sector ' num2str(j)];
    figu.TitlePosition = 'centertop'; 
    figu.FontSize = 12;
    figu.FontWeight = 'bold';
    darkred = [0.6350 0.0780 0.1840];
    cgreen  = [0.2431 0.5882 0.3176];

    
    in = in+1;subplot(nl,nc,in,'Parent',figu)
    scatter(z,x)
    xlabel('Firm productivity (z)');ylabel('Number of workers')
    title('Formal firms')

    in = in+1;subplot(nl,nc,in,'Parent',figu)
    plot(log(z),log(x))
    xlabel('Log of firm productivity (log z)');ylabel('Log of number of workers')
    title('Formal firms')
    
    if SO.FirmOptions > 1 % informal firms
        in = in+1;subplot(nl,nc,in,'Parent',figu)
        scatter(z,xi,[],darkred)
        xlabel('Firm productivity (z)');ylabel('Number of workers')
        title('Informal firms')

        in = in+1;subplot(nl,nc,in,'Parent',figu)
        plot(log(z),log(xi),'color',darkred)
        xlabel('Log of firm productivity (log z)');ylabel('Log of number of workers')
        title('Informal firms')
    end    
    
    if SO.FirmOptions > 2 % SIMPLES firms
        in = in+1;subplot(nl,nc,in,'Parent',figu)
        scatter(z,xs,[],cgreen)
        xlabel('Firm productivity (z)');ylabel('Number of workers')
        title('SIMPLES firms')

        in = in+1;subplot(nl,nc,in,'Parent',figu)
        plot(log(z),log(xs),'color',cgreen)
        xlabel('Log of firm productivity (log z)');ylabel('Log of number of workers')
        title('SIMPLES firms')
    end    
    
    

    %% Before entry
    
    % add this constant to prevent negative values while computing log of value function
    if SO.FirmOptions == 1
        cte = abs(min(VForm0)) + 1e-3;
    elseif SO.FirmOptions == 2
        cte = abs(min([VForm0;VInf0])) + 1e-3;   
    elseif SO.FirmOptions == 3
        cte = abs(min([VForm0;VInf0;VSimp0])) + 1e-3; 
    end
    
    % Subplot options
    nl = 1; % number of lines
    nc = 1; % number of columns
    in = 0; % subplot index
    if SO.FirmOptions > 1
        nc = SO.FirmOptions + 1;
    end
    
    ff = figure;
    ffigu = uipanel('Parent',ff,'BorderType','none'); 
    ffigu.Title = ['Before entry, file "TR_ValueFunc", sector ' num2str(j)];
    ffigu.TitlePosition = 'centertop'; 
    ffigu.FontSize = 12;
    ffigu.FontWeight = 'bold';

    in = in+1;subplot(nl,nc,in,'Parent',ffigu)
    plot(lz_grid,log(VForm0+cte));hold on;scatter(log(z0_f),log(Vf0_f+cte),[],'blue')
    xlabel('Log of firm productivity (log z)');ylabel('Lof of value fn before entry')
    legend({'Value function','Simulated'},'Location','southeast')
    title('Formal firms')
    
    if SO.FirmOptions > 1
        % Informal firms
        in = in+1;subplot(nl,nc,in,'Parent',ffigu)
        plot(lz_grid,log(VInf0+cte),'color',darkred);hold on;scatter(log(z0_f),log(Vi0_f+cte),[],darkred)
        xlabel('Log of firm productivity (log z)');ylabel('Lof of value fn before entry')
        legend({'Value function','Simulated'},'Location','southeast')
        title('Informal firms')
    end
    
    if SO.FirmOptions > 2
        % SIMPLES firms
        in = in+1;subplot(nl,nc,in,'Parent',ffigu)
        plot(lz_grid,log(VSimp0+cte),'color',cgreen);hold on;scatter(log(z0_f),log(Vs0_f+cte),[],cgreen)
        xlabel('Log of firm productivity (log z)');ylabel('Lof of value fn before entry')
        legend({'Value function','Simulated'},'Location','southeast')
        title('SIMPLES firms')
    end
    
    if SO.FirmOptions > 1
        in = in+1;subplot(nl,nc,in,'Parent',ffigu)
        plot(lz_grid,log(VForm0+cte));%hold on;scatter(log(z0_f),log(Vf0_f+cte),[],'blue'); 
        hold on;plot(lz_grid,log(VInf0+cte),'color',darkred);%hold on;scatter(log(z0_f),log(Vi0_f+cte),[],darkred);
        if SO.FirmOptions == 3
            hold on;plot(lz_grid,log(VSimp0+cte),'color',cgreen);%hold on;scatter(log(z0_f),log(Vs0_f+cte),[],cgreen);
        end
        hold on;plot([0,max(lz_grid)],log(cte)*[1,1],':','color','black');
        %hold on;yline(log(cte),':'); % MATLAB r2018b
        xlabel('Log of firm productivity (log z)');ylabel('Lof of value fn before entry')
        if SO.FirmOptions == 2
            %legend({'Value fn, formal','Simulated, formal','Value fn, informal','Simulated, informal','No entry'},'Location','northwest')
            legend({'Value fn, formal','Value fn, informal','No entry'},'Location','northwest')
            title('Formal x informal firms')
        elseif SO.FirmOptions == 3
            %legend({'Value fn, formal','Simulated, formal','Value fn, informal','Simulated, informal','Value fn, SIMPLES','Simulated, SIMPLES','No entry'},'Location','northwest')
            legend({'Value fn, formal','Value fn, informal','Value fn, SIMPLES','No entry'},'Location','northwest')
            title('Formal x informal x SIMPLES firms')
        end
    end
    
    %%% After entry
    % Subplot options
    nl = 2; % number of lines
    nc = 1; % number of columns
    in = 0; % subplot index
    if SO.FirmOptions > 1
        nc = SO.FirmOptions + 1;
    end
    
    f = figure;
    figu = uipanel('Parent',f,'BorderType','none'); 
    figu.Title = ['After entry, file "TR_ValueFunc", sector ' num2str(j)]; 
    figu.TitlePosition = 'centertop'; 
    figu.FontSize = 12;
    figu.FontWeight = 'bold';
    
    %%% Labor Demand
    in = in+1;subplot(nl,nc,in,'Parent',figu)
    plot(lz_grid,log(LaborPolicy.ellFormH));hold on;scatter(lz_f,log(Eh_f),[],'blue')
    xlabel('Log of firm productivity (log z)');ylabel('Log of number of workers');
    legend({'Policy function','Simulated'},'Location','southeast')
    title('Labor demand, formal firms')
    
    if SO.FirmOptions > 1
        in = in+1;subplot(nl,nc,in,'Parent',figu)
        plot(lz_grid,log(LaborPolicy.ellInf),'color',darkred);hold on;scatter(lz_f,log(Ei_f),[],darkred)
        xlabel('Log of firm productivity (log z)');ylabel('Log of number of workers');
        legend({'Policy function','Simulated'},'Location','southeast')
        title('Labor demand, informal firms')

        if SO.FirmOptions == 2
            in = in+1;subplot(nl,nc,in,'Parent',figu)
            plot(lz_grid,log(LaborPolicy.ellFormH));%hold on;scatter(lz_f,log(Eh_f),[],'blue');
            hold on;plot(lz_grid,log(LaborPolicy.ellInf),'color',darkred);%hold on;scatter(lz_f,log(Ei_f),[],darkred);
            xlabel('Log of firm productivity (log z)');ylabel('Log of number of workers')
            %legend({'Policy fn, formal','Simulated, formal','Policy fn, informal','Simulated, informal'},'Location','southeast')
            legend({'Policy fn, formal','Policy fn, informal'},'Location','southeast')
            title('Labor demand, formal x informal firms')
      
        elseif SO.FirmOptions > 2
            in = in+1;subplot(nl,nc,in,'Parent',figu)
            plot(lz_grid,log(LaborPolicy.ellSimp),'color',cgreen);hold on;scatter(lz_f,log(Es_f),[],cgreen)
            xlabel('Log of firm productivity (log z)');ylabel('Log of number of workers');
            legend({'Policy function','Simulated'},'Location','southeast')
            title('Labor demand, SIMPLES firms')
            
            if SO.FirmOptions == 3
                in = in+1;subplot(nl,nc,in,'Parent',figu)
                plot(lz_grid,log(LaborPolicy.ellFormH));%hold on;scatter(lz_f,log(Eh_f),[],'blue');
                hold on;plot(lz_grid,log(LaborPolicy.ellInf),'color',darkred);%hold on;scatter(lz_f,log(Ei_f),[],darkred);
                hold on;plot(lz_grid,log(LaborPolicy.ellSimp),'color',cgreen);%hold on;scatter(lz_f,log(Es_f),[],cgreen);
                xlabel('Log of firm productivity (log z)');ylabel('Log of number of workers')
                legend({'Policy fn, formal','Policy fn, informal','Policy fn, SIMPLES'},'Location','northwest')
                title('Labor demand, formal x informal x SIMPLES firms')
                
            else
                % Exporters
            end
        end
    end
    
    
    %%% Profits
    % ARRUMAR DAQUI PARA BAIXO
    
    % add this constant to prevent negative values while computing log of profits 
    %if SO.FirmOptions == 1
    %    ctep = abs(min(Ph_f)) + 1e-3;
    %elseif SO.FirmOptions == 2
    %    ctep = abs(min([Ph_f;Pi_f])) + 1e-3;   
    %elseif SO.FirmOptions == 3
    %    ctep = abs(min([Ph_f;Pi_f;Ps_f])) + 1e-3; 
    %end

    %in = in+1;subplot(nl,nc,in,'Parent',figu)
    %scatter(lz_f,log(Ph_f +ctep),[],'blue');
    %hold on;plot([0,max(lz_f)],log(ctep)*[1,1],':','color','black');
    
    in = in+1;subplot(nl,nc,in,'Parent',figu)
    scatter(exp(lz_f),Ph_f,[],'blue')
    xlabel('Firm productivity (z)');ylabel('Profits');
    %legend({'Policy function','Simulated'},'Location','northwest')
    title('Profits, formal firms')
    
    if SO.FirmOptions > 1
        in = in+1;subplot(nl,nc,in,'Parent',figu)
        scatter(exp(lz_f),Pi_f,[],darkred)
        xlabel('Firm productivity (z)');ylabel('Profits');
        %legend({'Policy function','Simulated'},'Location','northwest')
        title('Profits, informal firms')
        
        if SO.FirmOptions == 2
            ctep = abs(min([Ph_f;Pi_f])) + 1e-3;  % add this constant to prevent negative values while computing log of profits 
            in = in+1;subplot(nl,nc,in,'Parent',figu)
            scatter(lz_f,log(Ph_f +ctep),[],'blue');hold on;
            scatter(lz_f,log(Pi_f +ctep),[],darkred);
            hold on;plot([0,max(lz_f)],log(ctep)*[1,1],':','color','black');
            %hold on;yline(log(ctep),':'); % MATLAB r2018b
            xlabel('Log of firm productivity (log z)');ylabel('Log of profits')
            legend({'Profits, formal','Profits, informal','Zero profits'},'Location','northwest')
            title('Profits, formal x informal firms')
        
        elseif SO.FirmOptions > 2
            in = in+1;subplot(nl,nc,in,'Parent',figu)
            scatter(exp(lz_f),Ps_f,[],cgreen)
            xlabel('Firm productivity (z)');ylabel('Profits');
            %legend({'Policy function','Simulated'},'Location','northwest')
            title('Profits, SIMPLES firms')
            
            if SO.FirmOptions == 3
                ctep = abs(min([Ph_f;Pi_f;Ps_f])) + 1e-3;  % add this constant to prevent negative values while computing log of profits 
                in = in+1;subplot(nl,nc,in,'Parent',figu)
                scatter(lz_f,log(Ph_f +ctep),[],'blue');hold on;
                scatter(lz_f,log(Pi_f +ctep),[],darkred);hold on;
                scatter(lz_f,log(Ps_f +ctep),[],cgreen);
                hold on;plot([0,max(lz_f)],log(ctep)*[1,1],':','color','black');
                %hold on;yline(log(ctep),':'); % MATLAB r2018b
                xlabel('Log of firm productivity (log z)');ylabel('Log of profits')
                legend({'Profits, formal','Profits, informal','Profits, SIMPLES','Zero profits'},'Location','northwest')
                title('Profits, formal x informal x SIMPLES firms')

            else
                % Exporters
            end
        end
    
    
    end
    
    
end
end

end


