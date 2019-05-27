function Ystar = evaluateFunction(F,Xstar)
%EVALUATEFUNCTION   Evaluate function from symbolic handle and inputs. 

N = size(Xstar,2);
n = size(F,1);

X_sym = sort(symvar(F)).';
F_mf = matlabFunction(F,'vars',X_sym);

Ystar = zeros(n,N);

for k = 1:N
    Xtmp = num2cell(Xstar(:,k));
    Ystar(:,k) = F_mf(Xtmp{:});
end
    
end %function
