function [F,V,W,G] = generateFunction(m,n,d,r)
%GENERATEFUNCTION   Generates a coupled function and its decoupled representation.
%   [F,V,W,G] = generateFunction(n,m,d,r) generates randomly a function F,
%   having factors V, W and internal functions G.
%
%       m = number of inputs
%       n = number of outputs
%       d = degree
%       r = number of branches

U_sym = sym('u',[1 r]).';

V = randn(m,r);
W = randn(n,r);

G = sym(zeros(r,1));
for i = 1:r
    G(i) = sum((U_sym(i).^[0:d]).*randn(1,d+1));
end

X_sym = sym('x', [1 m]).';
F = expand(W * subs(G,U_sym,V' * X_sym));

end