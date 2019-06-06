%% DEMO file for MATLAB decoupling toolbox
% 
% Multivariate nonlinear vector function decoupling 
%
% version 0.2
% June 2019
%
% Philippe Dreesen, Jeroen De Geeter and Mariya Ishteva
% Vrije Universiteit Brussel

clear all; close all; clc; echo off;

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


%% A simple example

echo on % this echos also the commented lines to the Command Window

% Consider a simple coupled function F that we will decouple as

syms x1 x2;
syms u1 u2;

F = [ - 108*x1^3 + 108*x1^2*x2 - 18*x1^2 - 36*x1*x2^2 + 12*x1*x2 - 12*x1 + 1540*x2^3 + 1342*x2^2 - 68*x2 + 6; ...
        - 378*x1^3 + 378*x1^2*x2 - 63*x1^2 - 126*x1*x2^2 + 42*x1*x2 - 42*x1 - 2034*x2^3 - 1799*x2^2 + 110*x2 - 8];

% F has two inputs (x1 and x2) and two outputs (F(1) and F(2)) and admits a
% rank-two representation. The total degree d occurring of F is three.

m = 2; 
n = 2; 
r = 2; 

d = 3; 

% Press any key to continue...
pause

% It can be checked that the F has the following decoupled representation,
%
%       F(x) = Wo Go(Vo'*x), where x = [x1; x2], with
 
Vo = [-3  0 ;1 -8];
Wo = [ -2  3 ; -7 -4];
Go = [- 2*u1^3 +   u1^2 - 2*u1 ; - u2^3 + 7*u2^2 + 3*u2 + 2];

% Sanity check: putting everything together returns F as above:

Fo = expand(Wo*subs(Go,symvar(Go).',Vo'*sym('x', [1 m]).'))

% (For development, it is possible, using "generateFunction" to generate a 
% random coupled function starting from a decoupled representation.)

% Now let us try to recover from F the decoupled form. 
% Press any key to continue...
pause 




% Let us now try to unravel from the coupled function F, its decoupled 
% representation involving W, V and G. 

% First evaluate the coupled function F in a number N of sampling points.
% The function "constructDataset" evaluates F in N sampling points drawn 
% from a uniform distribution. 

N = 500;
[X,Y] = constructDataset(F,N);

% Decouple the function by performing CPD on the Jacobian tensor containing 
% the first-order information (evaluated at the N sampling points):

[W_J, V_J, G_J, output_J] = decoupleJ(F,X,r,d,Y);

% Although W_J and V_J look different from the original factors Wo and Jo,
% it can be verified that they are correct modulo scaling and permutation 
% (cf. 'essential uniqueness'). The cpderr should be (very) small.

cpderrorJ = cpderr({Wo,Vo},{W_J,V_J})

if ( norm(cpderrorJ) < 0.0001),
    disp('Good result!');
else
    warning('There seems to be a problem; local minimum? Try to re-run decoupleJ.'); 
end

% Finally, it is possible to reconstruct a function F_J from the decoupled
% representation, and 'visually' check the coefficients. 

F_J = vpa(expand(W_J*subs(G_J,symvar(G_J).',V_J'*sym('x', [1 m]).')),4)

% Press any key to continue...
pause 



% Decouple the same function using a SDF-based method.
% By default this method uses function evaluation, first- and second order
% information. Again we look at 'cpderr' to evaluate the error. 

[W_SDF012, V_SDF012, G_SDF012, sol_SDF012, output_SDF012] = decoupleSDF(F,X,r,d,Y);

cpderror012= cpderr({Wo,Vo},{W_SDF012,V_SDF012})

if ( norm(cpderror012) < 0.0001),
    disp('Good result!');
else
    warning('There seems to be a problem; local minimum? Try to re-run.'); 
end

F_SDF012 = vpa(expand(W_SDF012*subs(G_SDF012,symvar(G_SDF012).',V_SDF012'*sym('x', [1 m]).')),4)

% Press any key to continue...
pause 



% Decouple the same function using function evaluations (zeroth order) and
% first-order derivatives. Also show plots of the recovered internal
% functions G(i). Again we look at 'cpderr' to evaluate the error. 

[W_SDF01, V_SDF01, G_SDF01, sol_SDF01, output_SDF01] = decoupleSDF(F,X,r,d,Y, 'info', [0 1], 'plot', true);

cpderror01= cpderr({Wo,Vo},{W_SDF01,V_SDF01})

if ( norm(cpderror01) < 0.0001),
    disp('Good result!');
else
    warning('There seems to be a problem; local minimum? Try to re-run decoupleJ.'); 
end

F_SDF012 = vpa(expand(W_SDF01*subs(G_SDF01,symvar(G_SDF01).',V_SDF01'*sym('x', [1 m]).')),4)


%% print internal functions
Uo=Vo.'*X;
usymvars=sym('u', [1 r]);

figure;
for i=1:r,
    Go_mf = matlabFunction(Go(i),'vars', usymvars(i));
    subplot(2,r,i);
    plot(Uo(i,:), Go_mf(Uo(i,:)),'.')
    title('true') 
end
Ue=V_SDF01.'*X;
Ge=sol_SDF01.factors.G0;
for i=1:r,
    subplot(2,r,i+r);
    plot(Ue(i,:), Ge(:,i),'.')
    title('reconstructed') 
end



echo off


%% Questions
% - Fix m and n and experiment with incresing r (r = 1, 2, 3, ...).
% For each r, initialize a few times to get (close to) zero CPD-errors (cpderr). This is to avoid local minima.
% For which r are the factors correctly recovered (up to scaling and permutation)? 
% Fow which r is the decomposition unique?
% - Should we always run the most "powerful" algorithms? 
% Are the simpler algorithms faster in solving simpler problems?

