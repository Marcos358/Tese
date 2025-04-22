%//////////////////////////////////////////////////////////////////////////
% Indirect Tax Reform
%//////////////////////////////////////////////////////////////////////////
% TR_Profits
% Profits function
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [VA,Y,Profit,R,q,iota_k,Emp] = TR_ProfitsSimp(Emp,z,wage,p,piota,p_k,lambda_k,alpha,theta,taxY,taxW,RmaxSimp)

    % Before revenue cap:
    iota = (p / piota * theta * z .* Emp.^alpha).^(1/(1 - theta));
    q =  z .* Emp.^alpha .* iota.^theta;
    %q = ((theta * p / piota)^theta * z .* Emp.^alpha).^(1/(1 - theta));
    R = p * q;
    
    % Revenue cap
    Emp(R>RmaxSimp) = (RmaxSimp ./ (p * z(R>RmaxSimp) .* (theta*(1+taxW)*wage./(alpha*(1-taxY)*piota)).^theta) ).^(1./(alpha+theta));
    iota(R>RmaxSimp) = (theta*(1+taxW)*wage./(alpha*(1-taxY)*piota))* Emp(R>RmaxSimp);    
    R(R>RmaxSimp) = RmaxSimp;
    q = R / p;
    
    % Value added and profits:
    Y = (z .* Emp.^alpha).^(1./(1-theta)); % Real value added
    VA = R - piota * iota; % Nominal value added
    Profit = (1 - taxY)* VA - (1 + taxW)* wage* Emp;
    iota_k = (lambda_k ./ p_k) * piota .* iota;
end


