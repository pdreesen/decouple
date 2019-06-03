function [W,V,G,sol,output] = decoupleSDF_structureless(F,X,r,d,Y,varargin)
%DECOUPLESDF  Decouple a given function using Structured Data Fusion. 
%   [W,V,G,sol,output] = decoupleSDF(F,X,r,d,Y) returns a decoupled
%   representation of a given multivariate polynomial vector function F by
%   combining zeroth (function evaluation), first, and second-order 
%   derivative information in the sampling points X, while imposing that 
%   the internal functions G(i) are polynomial. This leads to a Structured
%   Data Fusion problem and is solved by tensorlab. 
%         F = multivariate vector function
%         X = sampling points to take function evaluations and derivatives
%         r = number of internal functions G
%         d = polynomial degree
%         Y = function evaluations at sampling points X

F = formula(F);
X_sym = sort(symvar(F)).';
N = size(X,2);
n = size(F,1);
m = size(X_sym,1);

p = inputParser;
addOptional(p, 'W', randn(n,r));
addOptional(p, 'V', rand(m,r));
addOptional(p, 'G', {randn(N,r),randn(N,r),randn(N,r)});
addOptional(p, 'info', [0 1 2]);
addOptional(p, 'plot', false);
parse(p,varargin{:});   
%'check_model'

% calculate function evaluation
Ystar = Y; 
%ystar = zeros(n,N);
%for i = [1:N]
%   ystar(:,i) = subs(F, X_sym, X(:,i));
%end

% calculate first order information of the function
if any(p.Results.info==1)
    J = jacobian(F,X_sym);
    J_mf = matlabFunction(J, 'vars', X_sym); 
    
    Jstar = zeros(n,m,N);
    for k = 1:N
        % Jstar(:,:,k) = subs(J, X_sym, X(:,k));
        Xtmp = num2cell(X(:,k));
        Jstar(:,:,k) = J_mf(Xtmp{:});
    end
end

% calculate second order information of the function
if any(p.Results.info==2)
    for i = 1:n, 
        H(i,:,:) = hessian(F(i), X_sym);
    end
    H_mf = matlabFunction( H, 'vars', X_sym );
    
    Hstar = zeros(n,m,m,N);
    for k = 1:N
        % Hstar(:,:,:,k) = subs(H, X_sym, X(:,k));
        Xtmp = num2cell(X(:,k));
        Hstar(:,:,:,k) = H_mf(Xtmp{:});
    end
end

model  =  struct;

model.variables.w = p.Results.W;
model.variables.v = p.Results.V;
for i = p.Results.info, model.variables.(sprintf('g%d',i)) = p.Results.G{i+1}; end

model.factors.W = {'w'};
model.factors.V = {'v'};
for i = p.Results.info, model.factors.(sprintf('G%d',i)) = {sprintf('g%d',i)}; end 

if any(p.Results.info==0), model.factorizations.tensorF.data = Ystar; end
if any(p.Results.info==1), model.factorizations.tensorJ.data = Jstar; end
if any(p.Results.info==2), model.factorizations.tensorH.data = Hstar; end

if any(p.Results.info==0), model.factorizations.tensorF.cpd = {'W', 'G0'}; end
if any(p.Results.info==1), model.factorizations.tensorJ.cpd = {'W', 'V', 'G1'}; end
if any(p.Results.info==2), model.factorizations.tensorH.cpd = {'W', 'V', 'V', 'G2'}; end

%model.factorizations.tensorF.weight = 3;
%model.factorizations.tensorJ.weight = 2;
%model.factorizations.tensorH.weight = 1;    

%sdf_check(model, 'print');

% tolerances and solution
options.TolFun = eps^2; %relative function tolerance
options.TolX = eps;     %step size tolerance
options.maxiter = 300; %2000;
options.cgmaxiter = 500;
options.Display = 50; %10;  %default is zero
% eps = 2.2204e-16

[sol, output] = ccpd_nls(model, options);
%[sol, output] = sdf_nls(model, options);
%[sol, output] = sdf_minf(model, options);

% Output
W = sol{1};
V = sol{2};

Z = V'*X;

% construct block diagonal matrix W
block_diagonal_W = zeros(N*n,N*r);
for i = 0:N-1
    block_diagonal_W(1+i*n:i*n+n,1+i*r:i*r+r) = W;
end

% construct block vandermonde-like matrix
Xk = zeros(N*r,r*(d+1));
for k = 1:N
    for j = 1:r
        for i = 0:d
            Xk((k-1)*r+j,i+1+(j-1)*(d+1)) = Z(j,k)^i;
        end
    end
end

% calculate the coefficients of the decoupled functions
Z = reshape(Y,N*n,1); %Y(:)
C = pinv(block_diagonal_W*Xk)*Z;

%build the decoupled functions
G = sym(zeros(r,r*(d+1)));
symbolic_vars = sym('u',[1 r]); 
for j = 1:r
    for i = 0:d
        G(j,i+1+(j-1)*(d+1)) = symbolic_vars(j)^i;
    end
end
G = G * C;


% if (p.Results.plot)
%     for i = p.Results.info
%         figure
%         for j = 1:r
%             subplot(1,r,j);
%             eval(sprintf('plot((V(:,j)'' * X)'', sol.factors.G%i(:,j), ''*'')', i)); 
%             %set(gca,'YTick',[])
%             %set(gca,'XTick',[])
%         end
%     end
% end
    
end %function