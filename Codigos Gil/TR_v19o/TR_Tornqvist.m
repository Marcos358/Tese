%%% Tornqvist index number

function Y_agg = TR_Tornqvist(Y_k,shY_k)

Nperiods = size(Y_k,1) - 1
secnum = size(Y_k,2)
%
logY_agg = zeros(Nperiods,1)
Y_agg = zeros(Nperiods+1,1)
Y_agg(1) = 1

for t=1:Nperiods
    for j=1:secnum
        logY_agg(t) = logY_agg(t)+(log(Y_k(t+1,j))-log(Y_k(t,j)))*0.5*(shY_k(t,j)+shY_k(t+1,j))        
    end
    if max(shY_k(t+1,:))>0.999
        Y_agg(t+1) = 1;
    else
        Y_agg(t+1) = Y_agg(t)*exp(logY_agg(t))    
    end
end

end