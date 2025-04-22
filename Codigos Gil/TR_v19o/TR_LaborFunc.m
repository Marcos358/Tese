%//////////////////////////////////////////////////////////////////////////
% Indirect Tax Reform
%//////////////////////////////////////////////////////////////////////////
% TR_LaborFunc
% Finds the aggregate labor demand policy function
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


function [LaborPolicy]=TR_LaborFunc(SO,alpha,theta,tauY,tauW,lz_grid,wage,p,piota,j,LmaxInf,tauYSimp,tauWSimp)

LaborPolicy.name = "Policy functions for firms after entry";

z_grid = exp(lz_grid);

%--------------------------------------------------------------------------
% Formal firms: 
%--------------------------------------------------------------------------

% Labor demand
z = z_grid;
ell_fn = @(taxW,taxY,wage)  ( (theta / piota)^theta * (alpha * (1-taxY) / ((1+taxW)*wage))^(1-theta) * p* z ).^(1 / (1 - alpha - theta));

taxY = tauY;
taxW = tauW;

%--------------------------------------------------------------------------------------------------------
% Formal firms, main tax code (H):
%--------------------------------------------------------------------------------------------------------

x = ell_fn(taxW,taxY,wage);
x(x<0)=0;

LaborPolicy.ellFormH = x;
%LaborPolicy.ellFormH = round(x);

if j > 1
    
    if SO.FirmOptions > 2
        %--------------------------------------------------------------------------------------------------------
        % Formal SIMPLES firms:
        %--------------------------------------------------------------------------------------------------------

        % Firm maximization before revenue cap
        taxY = tauYSimp;
        taxW = tauWSimp;
        xs = ell_fn(taxW,taxY,wage);
        xs(xs<0)=0;
        
        LaborPolicy.ellSimp = xs;
        %LaborPolicy.ellSimp = round(xs);
    end
    
    if SO.FirmOptions > 1
        %--------------------------------------------------------------------------
        % Informal firms
        %--------------------------------------------------------------------------
        
        taxY = 0;
        taxW = 0;
        xi = ell_fn(taxW,taxY,wage);
        xi(xi<0)=0;
        xi(xi>LmaxInf)=LmaxInf;
    
        LaborPolicy.ellInf = xi;
        %LaborPolicy.ellInf = round(xi);
    end
    
end

end
