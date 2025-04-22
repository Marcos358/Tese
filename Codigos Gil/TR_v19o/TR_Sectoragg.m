%//////////////////////////////////////////////////////////////////////////
% Indirect Tax Reform
%//////////////////////////////////////////////////////////////////////////

function [Le,LE,VA,Y,Q,Iota_k] = TR_Sectoragg(x,Id_f,E_f,VA_f,Y_f,q_f,iota_kf)

Id_f(isnan(Id_f)) = 0;

Le = Id_f.*E_f;
Le(isnan(Le)) = 0;
LE = sum(Le)*x;
%LS = round(LS);

va = Id_f.*VA_f;
va(isnan(va)) = 0;
VA = sum(va)*x;

y = Id_f.*Y_f;
y(isnan(y)) = 0;
Y = sum(y)*x;

q = Id_f.*q_f;
q(isnan(q)) = 0;
Q = sum(q)*x;

J = size(iota_kf,2);
Iota_k = zeros(1,J);
for k=1:J
    iota_k = Id_f.*iota_kf(:,k);
    iota_k(isnan(iota_k)) = 0;
    Iota_k(k) = sum(iota_k)*x;
end

end
