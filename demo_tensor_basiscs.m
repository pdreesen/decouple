%% DEMO file for MATLAB decoupling toolbox
% 
% Basic tensor operations

echo off
clear variables
close all
clc

echo on

%% Preliminaries
% check that decoupling toolbox is present
if exist('checkfortensorlab') ~= 2,
    warning('Decoupling toolbox not found. Add directory "functions" to the working path'); 

    % manually add path to the MATLAB path
    tlpath = uigetdir(pwd,'Select "functions" path');

    if exist('tlpath'), addpath(tlpath); end
end

if exist('checkfortensorlab')==2, disp('decoupling toolbox found!'); end

% check that tensorlab toolbox is present
if checkfortensorlab,
    disp('tensorlab found!'); 
end



%% Matrix SVD

% Let us construct a Matrix A 
% as the product of three known matrices U, S and V'

U = [1, 1; 3, 1]; % first row of U is composed of 1s
S = diag([3 1]); % diagonal matix
V = [1, 1; -1, 2]; % first row of V is composed of 1s
A = U * S * V';

% The SVD of the matrix A can be computed in the following way

[U_svd, S_svd, V_svd] = svd(A);

% We can then reconstruct A as the product of  U_svd, S_svd and V_svd'

A_svd = U_svd * S_svd * V_svd';

% Let us convince ourselves that A and A_svd are the same

A, A_svd

if (norm(A - A_svd) < 1e-8)
    disp('A and A_svd are the same.');
else
    disp('A and A_svd are different.')
end


% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The original factors and the ones obtained from the SVD 
% are different though.

U, U_svd

% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S, S_svd

% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V, V_svd

% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Tensor CPD

% Let us now move on to the tensor case.

% We will reuse the matrices U and V and will define a third matrix W. 

W = [1, 1; 0, 3]; % first row of W is composed of 1s

% We construct a tensor T from the matrices V, W, and W

T = cpdgen({U, V, W});

% and compute its rank with the rankest function. The rank should be 2.

R = rankest(T)

% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We will now decompose the tensor T
% and try to reconstruct T from these factors.

res = cpd(T, 2);

U_ten = res{1};
V_ten = res{2};
W_ten = res{3};

T_ten = cpdgen({U_ten, V_ten, W_ten});

% If this number is very small, then the have a good reconstruction:

norm(T(:) - T_ten(:))

% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Let us now compare the obtained factors with the original ones

U, U_ten

% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V, V_ten

% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W, W_ten

% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% At first sight, the factors look different from the original ones.
% However, if we normalize them, so that the first rows are composed of 1s,
% as is the case for the original factors, 
% then we see that the new factors are in fact the same as the original ones, 
% possibly up to permutation of the columns


U, U_ten * diag(1./U_ten(1,:))

% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V, V_ten * diag(1./V_ten(1,:))

% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W, W_ten * diag(1./W_ten(1,:))

% Great!

% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Can we apply the same normalization trick on SVD-factos?
% --> We can, but the factors are still not the same

U, U_svd * diag(1./U_svd(1,:))
V, V_svd * diag(1./V_svd(1,:))

% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Are you convinced that tensors are the way to go?

% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If yes, great! Go to the next demo.
% If no, no problem. Change the matrices U, V and W and try this demo again.

% Tip: If you keep the first rows of the original matrices equal to [1 1], 
% it will be easier to compare the results.

% For advanced users: 
% you can also change the dimensions of the matrices U, V, and W.
% Can you modify this script, so that it works for a more general case
% (not only for 2 x 2 matrices, but n x r, for example)?

echo off
