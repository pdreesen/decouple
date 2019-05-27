%% DEMO file for decoupling methods toolbox
% 
% Authors: Jeroen De Geeter and Philippe Dreesen
% (c) Vrije Universiteit Brussel (VUB), June 2018. 
% 
% version 0.1: For internal use only. 


%% Clear workspace, close windows and clear the Command Window
clear all; close all; clc;


%% Preliminaries: path to this toolbox and to tensorlab 

% add path to tensorlab folder to working path

%check what is already present on the working path
if exist('sdf_check') ~= 2, % this is a function from tensorlab
    warning('This code requires tensorlab. Please download tensorlab from "http://www.tensorlab.net/" (click <a href="http://www.tensorlab.net">here</a>).');

    % to manually add tensorlab path to the MATLAB path
    tlpath = uigetdir(pwd,'Select tensorlab path');
    if exist('tlpath'), addpath(tlpath); end
end

% add decoupling functions to path
addpath(pwd); %this is a dangereous assumption
addpath('./functions')
% TODO: automatically figure out where
% the demo.m file is loacated, and add the directory to MATLAB path
% test = mfilename('fullpath'); 




%% A simple example

syms x1 x2;
syms u1 u2;

echo on % this echos also the commented lines to the Command Window

% Consider a simple coupled function F that we will decouple as

syms x1 x2;

F = [ - 108*x1^3 + 108*x1^2*x2 - 18*x1^2 - 36*x1*x2^2 + 12*x1*x2 - 12*x1 + 1540*x2^3 + 1342*x2^2 - 68*x2 + 6; ...
        - 378*x1^3 + 378*x1^2*x2 - 63*x1^2 - 126*x1*x2^2 + 42*x1*x2 - 42*x1 - 2034*x2^3 - 1799*x2^2 + 110*x2 - 8]

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
 
Vo = [-3  0 ;1 -8],
Wo = [ -2  3 ; -7 -4], 
Go = [- 2*u1^3 +   u1^2 - 2*u1 ; - u2^3 + 7*u2^2 + 3*u2 + 2],

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
[U,Y] = constructDataset(F,N);

% Decouple the function by performing CPD on the Jacobian tensor containing 
% the first-order information (evaluated at the N sampling points):

[W_J,V_J,G_J,output_J] = decoupleJ(F,U,r,d,Y);

% Although W_J and V_J look different from the original factors Wo and Jo,
% it can be verified that they are correct modulo scaling and permutation 
% (cf. 'essential uniqueness'). The cpderr should be (very) small.

cpderrorJ = cpderr({Wo,Vo},{W_J,V_J})
if ( norm(cpderrorJ) > 0.0001), 
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
[W_SDF012,V_SDF012,G_SDF012,output_SDF012] = decoupleSDF(F,U,r,d,Y);
cpderror012= cpderr({Wo,Vo},{W_SDF012,V_SDF012})
if ( norm(cpderror012) > 0.0001), 
    warning('There seems to be a problem; local minimum? Try to re-run decoupleSDF.'); 
end

% Press any key to continue...
pause 



% Decouple the same function using function evaluations (zeroth order) and
% first-order derivatives. Also show plots of the recovered internal
% functions G(i). Again we look at 'cpderr' to evaluate the error. 
[W_SDF01, V_SDF01, G_SDF01, output_SDF01] = decoupleSDF(F,U,r,d,Y, 'info', [0 1], 'plot', true);
cpderror01= cpderr({Wo,Vo},{W_SDF01,V_SDF01})
if ( norm(cpderror01) > 0.0001), 
    warning('There seems to be a problem; local minimum? Try to re-run decoupleSDF.'); 
end




% End of DEMO
echo off




