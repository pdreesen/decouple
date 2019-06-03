function [W,V,G,sol,output] = decoupleSDF(F,X,r,d,Y,varargin)
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

%TODO: 
% * add parameter with regards to accuracy
% * passing ystar is an optional parameter
% * detect polynomial degree automatically

F = formula(F);
X_sym = sort(symvar(F)).';
N = size(X,2);
n = size(F,1);
m = size(X_sym,1);

p = inputParser;
addOptional(p, 'W', randn(n,r));
addOptional(p, 'V', randn(m,r));
addOptional(p, 'C', randn(d+1,r));
addOptional(p, 'info', [0 1 2]);
addOptional(p, 'plot', false);
addOptional(p, 'structure', @helperMethod);
parse(p,varargin{:});   

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
        Xtmp = num2cell(X(:,k));
        Jstar(:,:,k) = J_mf(Xtmp{:});
    end
end


% calculate second order information of the function
if any(p.Results.info==2)
    for i = 1:n 
        H(i,:,:) = hessian(F(i), X_sym);
    end
    H_mf = matlabFunction( H, 'vars', X_sym );
    
    Hstar = zeros(n,m,m,N);
    for k = 1:N
        Xtmp = num2cell(X(:,k));
        Hstar(:,:,:,k) = H_mf(Xtmp{:});
    end
end

model  =  struct;

for i = p.Results.info eval(sprintf('g%d = @(z,task) p.Results.structure(z,task,X,d,%d);',i,i)); end
select_v = @(z,task) struct_select(z,task,2);
select_c = @(z,task) struct_select(z,task,1);

model.variables.cv = {p.Results.C,p.Results.V}; %C,V
model.variables.w = p.Results.W;    

model.factors.W = {'w'};
model.factors.V = {'cv', select_v};
model.factors.C = {'cv', select_c};
for i = p.Results.info, eval(sprintf('model.factors.G%d = {''cv'', g%d};',i,i)); end

if any(p.Results.info==0), model.factorizations.tensorF.data = Ystar; end
if any(p.Results.info==1), model.factorizations.tensorJ.data = Jstar; end
if any(p.Results.info==2), model.factorizations.tensorH.data = Hstar; end

if any(p.Results.info==0), model.factorizations.tensorF.cpd = {'W', 'G0'}; end
if any(p.Results.info==1), model.factorizations.tensorJ.cpd = {'W', 'V', 'G1'}; end
if any(p.Results.info==2), model.factorizations.tensorH.cpd = {'W', 'V', 'V', 'G2'}; end

%model.factorizations.tensorF.weight = 1;
%model.factorizations.tensorJ.weight = 1;
%model.factorizations.tensorH.weight = 1;    

%sdf_check(model, 'print');

% tolerances and solution
options.TolFun = eps^2; %relative function tolerance
options.TolX = eps;     %step size tolerance
options.maxiter = 300; %2000;
options.cgmaxiter = 500;
options.Display = 50; %10;  %default is zero
% eps = 2.2204e-16

%[sol, output] = ccpd_nls(model, options);
[sol, output] = sdf_nls(model, options);
%[sol, output] = sdf_minf(model, options);

% Output
W = sol.factors.W;
V = sol.factors.V;
Ce = sol.factors.C; % pdreesen: why not 'C'?

%cpderr({W, V, G1, G2}, {W, V, G1e, G2e}) % zero error (if converged)

symbolic_vars = sym('u', [1 r]); 
G = sym(zeros(r, 1));
for i=1:r
    G(i) = symbolic_vars(i).^[0:d] * Ce(:,i);
end

if (p.Results.plot)
    for i = p.Results.info
        figure
        for j = 1:r
            subplot(1,r,j);
            eval(sprintf('plot((V(:,j)'' * X)'', sol.factors.G%i(:,j), ''*'')', i)); 
            %set(gca,'YTick',[])
            %set(gca,'XTick',[])
        end
    end
end
    
end %function
