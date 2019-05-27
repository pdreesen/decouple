function Ystar = evaluateDecoupledFunction(W,V,G,Xstar)
%EVALUATEDECOUPLEDFUNCTION   Evaluate function from factors and inputs. 

N = size(Xstar,2);
n = size(W,1);
m = size(V,1);

X_sym = sym('x', [1 m]).';

Ystar = zeros(n,N);

F = expand(W*subs(G,sort(symvar(G)).',V.' * X_sym));
F_mf = matlabFunction(F, 'vars', X_sym);

for k = [1:N]
    Xtmp = num2cell(Xstar(:,k));
    Ystar(:,k) = F_mf(Xtmp{:});
end

end

