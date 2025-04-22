%//////////////////////////////////////////////////////////////////////////
% Indirect Tax Reform
%//////////////////////////////////////////////////////////////////////////
% TR_Profits
% Profits function
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function [VA,Y,Profit,R,q,iota_k] = TR_Profits(Emp,z,wage,p,piota,p_k,lambda_k,alpha,theta,taxY,taxW)

    % Value added and profits:
    iota = (p / piota * theta * z .* Emp.^alpha).^(1/(1 - theta));
    q =  z .* Emp.^alpha .* iota.^theta;
    %q = ((theta * p / piota)^theta * z .* Emp.^alpha).^(1/(1 - theta));
    R = p * q;
    Y = (z .* Emp.^alpha).^(1./(1-theta));  % Real value added
    VA = (1-theta)*R;                       % Nominal value added
    %VA = R - piota * iota; 
    Profit = (1 - taxY)* VA - (1 + taxW)* wage* Emp;
    iota_k = (lambda_k ./ p_k) * piota .* iota;
    %iota_k = theta * p * (lambda_k ./ p_k) .* q;
        
end