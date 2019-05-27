function B = constructVariablesPolynomial(base, d, k)
%CONSTRUCTVARIABLESPOLYNOMIAL   Generates monomial evaluations and derivatives.
%   B = constructVariablesPolynomial(base, d, k)

numberOfEquations = length(base);
numberOfCoefficients = length([0-k:d-k]);

B = repmat(base,1,numberOfCoefficients).^repmat([0-k:d-k],numberOfEquations,1);
for i = [1:k]
    factors = zeros(1,d+1);
    factors(1,i+1:end) = [1:(d+1)-i].';
    B = B.*repmat(factors,numberOfEquations,1);
end
    
end %function