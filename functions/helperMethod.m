function [X, state] = helperMethod(z,task,U,d,k) %z = {C,V,W}
%HELPERMETHOD   structure function for SDF to impose polynomial internal functions G.
%   helperMethod(C,U,V,d,k)
%   helperMethod(z,task,U,d,k)
% 
%      - k determines the derivative = 0,1,2

C = z{1}; V = z{2};
r = size(V,2); %Dit kan ook bepaald worden op basis van V of V

if nargin < 2, task = []; end
right = ~isempty(task) && isfield(task, 'r') && ~isempty(task.r);
left  = ~isempty(task) && isfield(task, 'l') && ~isempty(task.l);

if ~isempty(task) && isfield(task, 'persistent')
   N = task.persistent.N;
   state = [];
else % not present, compute
   N = size(U,2);
   state.persistent.N = N;
end

VU = V' * U;

% function evaluation
if ~left && ~right
    X = zeros(size(U,2), r);
    for i = [1:r]
        B = constructVariablesPolynomial(VU(i,:)', d, k);
        X(:,i) = B * C(:,i);
    end

% right jacobian vector product
elseif right
    % TODO: optimize
    % derivative w.r.t. C
    X1 = [];
    for i = [1:r]
        B = constructVariablesPolynomial(VU(i,:)', d, k);
        X1 = [X1 (B * task.r{1}(:,i))];
    end

    % derivatvive w.r.t. V
    X2 = []; %zeros(size(V)); %dit zal hier ook niet meer kloppen
    for j = [1:r] %for each of the branches
        temp = [];
        for i = [1:size(V,1)] %for each of the input dimensions
            B = constructVariablesPolynomial(VU(j,:)', d, k+1);
            B = bsxfun(@times, B, U(i,:)');
            temp = [temp (B * C(:,j))];
        end
        X2 = [X2 (temp * task.r{2}(:,j))];
    end

    X = X1 + X2;

% left jacobian vector product
elseif left
    % derivative w.r.t. C
    % TODO: optimize
    cl = [];
    for i = [1:r]
        B = constructVariablesPolynomial(VU(i,:)', d, k);
        cl = [cl (B.' * task.l(:,i))];
    end

    % derivatvive w.r.t. V
    tl = zeros(size(V));
    for j = [1:r] %for each of the branches
        for i = [1:size(V,1)] %for each of the input dimensions
            B = constructVariablesPolynomial(VU(j,:)', d, k+1);
            B = bsxfun(@times, B, U(i,:)');
            tl(i,j) = conj(B * C(:,j)).' * task.l(:,j);
        end
    end

    X = {cl, tl};
end

end
