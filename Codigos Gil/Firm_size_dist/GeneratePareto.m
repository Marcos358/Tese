function sample = GeneratePareto(X1,Xm,xi,data,formal)

N = length(X1);
nu = Xm*(ones(N,1) - X1).^(-1/xi);

sample = log(nu);



histogram(sample, 'Normalization','probability'); hold on
histogram(data, 'Normalization','probability')
end