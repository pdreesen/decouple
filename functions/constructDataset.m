function [X,Y] = constructDataset(F,N) 
%CONSTRUCTDATASET   Generate data from given multivariate vector function.
%   [X,Y] = constructDataset(F,N) generates N inputs X and outputs Y from 
%   a given multivariate vector function F. The inputs X are drawn from a 
%   uniform distribution on the interval [0,1]. 

syms = sort(symvar(F)).';
m = size(syms,1);
n = size(F,1);

X = rand(m,N);
Y = zeros(n,N);
for i = [1:N]
    Y(:,i) = subs(F, syms, X(:,i));
end

end %function