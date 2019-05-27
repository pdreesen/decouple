function Fn = estimateFunctionFromNoisyData(F,Yn,U)
%ESTIMATEFUNCTIONFROMNOISYDATA   Estimate 'noisy' version of a function. 
%   Fn = estimateFunctionFromNoisyData(F,Yn,U) returns a symbolic representation
%   of F that has been obtained by taking noisy measurements of the outputs
%   of F and then performing parameter estimation to find the coefficients.
%
%   The degree of the function is assumed to remain the same.

N = length(U);
m = size(symvar(F).',1);
n = size(F,1);

% Determine max degree
d = 0;
for i = [1:size(F,1)]
   d = max(feval(symengine, 'degree', F(i)),d);
end

% Create all possible monomials with degree up to d
X_sym = sym('x', [1 m]).';
monomial_powers = generate_mons_full(m, d);
monomials = sym(zeros(size(monomial_powers,1), 1));
for i = [1: size(monomial_powers,1)]
    monomials(i,1) = prod((X_sym.').^monomial_powers(i,:));
end

% Make sure that our system is not underdetermined
if N < length(monomial_powers) * n
   error('Trying to solve underdetermined system, add more measurements');
end

% Calculate the coefficients
F_mf = matlabFunction(monomials, 'vars', X_sym);
K = zeros(n * N, size(monomials,1) * n);
for i = [1:N]
    for j = [1:n]        
        Utmp = num2cell(U(:,i));
        K(j + (i-1)*n, size(monomials,1)*(j-1)+1:size(monomials,1) * j) = F_mf(Utmp{:});
    end
end

C = pinv(K) * Yn(:);

% Create new (polynomial) function using the calculated coefficients
Fn = sym(zeros(n,1));
for i = [1:n]
    Fn(i,:) = C(size(C, 1)/n * (i-1) + 1:size(C, 1)/n * (i-1) + size(C, 1)/n)' * monomials;
end

end
