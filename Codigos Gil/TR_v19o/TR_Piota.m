% Intermediate inputs price index, by sector

function Prices = TR_Piota(Prices,FixPar)

p_k = Prices.p_j;
J = length(p_k);
piota_j = zeros(1,J);
for j=1:J
    lambda_k = FixPar.lambda_jk(j,:);
    Piota = (p_k ./ lambda_k).^lambda_k;
    piota_j(j) = prod(Piota);
end
Prices.piota_j = piota_j;

end