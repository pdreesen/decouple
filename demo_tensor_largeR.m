echo off
clear variables
close all
clc

echo on
%% CPD wih r (rank) > n (dimension)
% This is a more challenging case than the case r <= n.

n = 4;
r = 5;

% Let us define three n X r matrices U, V, and W.

U = [1, 1, 1, 1, 1; ... % first row of U is composed of ones
    2, 3, 4, 5, 6; ...
    1, -1, 3, 1, 4; ...
    -2, -1, -7, 2, 5];

V = [1, 2, 1, 1, 1; ... % first row of V is composed of ones
    2, 1, 2, 5, 3; ...
    -1, -1, 5, 1, 3; ...
    2, -2, -1, -2, 2];

W = [1, 1, 1, 1, 1; ... % first row of W is composed of ones
    1, -2, 3, 7, 6; ...
    -1, -1, 2, -1, 4; ...
    -1, 1, 7, -2, 4];

% We construct a tensor T from the matrices V, W, and W:

T = cpdgen({U, V, W});

% and compute its rank with the rankest function. The rank should be 5.

R = rankest(T)

% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We will now decompose the tensor T and reconstruct T from its factors.

T_ten = zeros(size(T));
i = 0;

while  (norm(T(:) - T_ten(:)) > 1e-8 && i < 20)
    % We have a good reconstruction if norm(T(:) - T_ten(:)) is very small.
    % Otherwise we have probably reached a local minimum.
    
    res = cpd(T, 5);
    
    U_ten = res{1};
    V_ten = res{2};
    W_ten = res{3};
    
    T_ten = cpdgen({U_ten, V_ten, W_ten});
    i = i+1;
end

if i < 20
    % It took i initialization(s) to compute the decomposition
    i
else
    i
    % The problem seems to be too difficult.
    return
end

% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Let us now compare the original factors
% with the normalized versions of the new factors


U, U_ten * diag(1./U_ten(1,:))

% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V, V_ten * diag(1./V_ten(1,:))

% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W, W_ten * diag(1./W_ten(1,:))

% Are the factors the same, up to permutation of the columns?

% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix decompositions cannot do this!

% Are you convinced that tensors are the way to go?

% Press any key to continue.
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If yes, great! Go to the next demo.
% If no, no problem. Change the matrices U, V and W and try this demo again.

% Tip: If you keep the first rows of the original matrices as [1 ... 1],
% it will be easier to compare the results.

% For advanced users:
% you can also change the dimensions of the matrices U, V, and W.
% Can you modify this script, so that it works for a more general case
% (not only for 4 x 5 matrices, but n x r, for example)?

% For a fixed n, how large can r be?

echo off
