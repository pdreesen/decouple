function [W,V,G,output] = decoupleJ(F,X,r,d,Y,varargin)
%DECOUPLEJ   Decouple a given function using Jacobian information. 
%   [W,V,G,output] = decoupleJ(F,X,r,d,Y) returns a decoupled
%   representation of a given multivariate polynomial vector function F by
%   collecting first-order derivative information in the sampling points X,
%   and performing a CPD on the obtained Jacobian tensor. 
%
%   Options:
%       - init: so either you enter gevd, or you give {[],[],[]}
%       - plot: plot the values contained in the H variables

% ---- 
% check options
% ----
p = inputParser;
addOptional(p,'init',{[],[],[]}); %W,V,H
addOptional(p,'plot',false);
parse(p,varargin{:});
options = p.Results;

% ---- 
% variable initialization
% ---- 
F = formula(F);
X_sym = sort(symvar(F)).';            % get all required symbolic variables
N = size(X,2);                  % number of points
n = size(F,1);                  % number of output functions
m = size(X_sym,1);              % number of inputs

cpdOptions.TolFun = eps^2;      %relative function tolerance
cpdOptions.TolX = eps;          %step size tolerance
cpdOptions.maxiter = 2000;
cpdOptions.Compression = false;

if iscell(options.init)
    % check if three cell arrays are provided {W,V,H}, throw error
    % otherwise
    if length(options.init) ~= 3 
        error('When initializing init option with cell, 3 cell arrays are required (W,V,H)'); 
    end
    
    if isempty(options.init{1})
        options.init{1} = randn(n,r);
    elseif ~isequal(size(options.init{1}), [n,r])
        error('dimensions of input for W are wrong');
    end
    
    if isempty(options.init{2})
        options.init{2} = randn(m,r);
    elseif ~isequal(size(options.init{2}), [m,r])
        error('dimensions of input for V are wrong');
    end
    
    if isempty(options.init{3})
        options.init{3} = randn(N,r);
    elseif ~isequal(size(options.init{3}), [N,r])
        error('dimensions of input for H are wrong');
    end
    options.gevd = false;
elseif (ischar(options.init) && strcmp(options.init, 'gevd'))
    options.gevd = true;
else 
    error('chosen setting for init is unknown');
end

% ---- 
% algorithm
% ---- 
% calculate first order information of the function
J = jacobian(F,X_sym);
J_mf = matlabFunction(J,'vars',X_sym); 

% create tensor with first order information of the samples
T = zeros(n,m,N);
for k = 1:N
    Xtmp = num2cell(X(:,k));
    T(:,:,k) = J_mf(Xtmp{:});
end

% decompose the tensor
if (options.gevd) 
    [CPD,output] = cpd(T,r,cpdOptions);
else
    [CPD,output] = cpd(T,options.init,cpdOptions);
end

W = CPD{1};
V = CPD{2};
H = CPD{3};

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

% plot the derivative of the internal functions
if (p.Results.plot)
    figure
    number_of_plots = size(G,1);
    for j = 1:number_of_plots
        subplot(1,number_of_plots,j);
        plot((V(:,j)'*X)',H(:,j),'*');
        %set(gca,'YTick',[])
        %set(gca,'XTick',[])
    end
end
    
end%function
