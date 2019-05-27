function [X,Yn,Ynd] = constructDataset2(F,N,SNR)
%CONSTRUCTDATASET   Generate data from given multivariate vector function.
%   [X,Y] = constructDataset(F,N) generates N inputs X and outputs Y from 
%   a given multivariate vector function F. The inputs X are drawn from a 
%   uniform distribution on the interval [0,1]. 

% X-values
% Yn : Y values met noise op
% Ynd : afgeleide 
% Yndd: 2de orde afgeleide

% TODO: add extra noise 

syms = sort(symvar(F)).';
m = size(syms,1);
n = size(F,1);
step = 0.000001;

% generate points
X = rand(m,N);          %[1 2 3; 4 5 6; 7 8 9;];
Y = zeros(n,N);
for i = 1:N, Y(:,i) = subs(F, syms, X(:,i)); end

% generate points for differential
X2 = zeros(m,m*N);
for i = 1:N         %elk punt in X
    for j = 1:m     %elke variable
        point = X(:,i);
        point(j) = point(j) + step;
        X2(:,(i-1)*m+j) = point;
    end
end
Y2 = zeros(n,m*N);
for i = 1:m*N, Y2(:,i) = subs(F, syms, X2(:,i)); end

% add noise
Yn12 = [Y Y2];
for i = 1:n, Yn12(i,:) = noisy(Yn12(i,:),SNR); end
Yn = Yn12(:,1:N);
Yn2 = Yn12(:,N+1:end);
clear Yn12;

% calculate derivatives
Ynd = zeros(n,m,N);
for i = 1:N         % for each point
    for j = 1:n     % for each variable
        Ynd(:,j,i) = (Yn2(:,(i-1)*m+j) - Yn(:,i)) ./ step;
    end
end

end %function