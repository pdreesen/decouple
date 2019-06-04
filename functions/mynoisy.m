function [y, n] = mynoisy(x, SNR, dist),
% Add noise with given SNR to a vector. 
% 
% y = mynoisy(x, SNR, dist);
% 
% such that 20*log10(norm(x)/norm(y-x)) = SNR
% 
% x is a given vector
% SNR is desired signal-to-noise ratio [20]
% dist is desired noise distribution [randn]
% y is noisy vector y = x + n
% n is the added noise itself
% 

if nargin<2 || isempty(SNR), SNR = 20; end
if nargin<3 || isempty(dist), dist = @randn; end

if ~isvector(x), warning('first argument should be vector!'); end
if ~isreal(x), error('first argument should be real'); end

n = dist(size(x));
scale = norm(x)/norm(n)*10^(-SNR/20); 
n = scale*n; 
y=x+n;

end % function
