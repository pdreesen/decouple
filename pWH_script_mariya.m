%% initialize
clear all; close all;
addpath('/home/noot/Dropbox/MATLAB/tensorlab3')
s=rng
% load testss
% rng(s);

mall = 5
mg1 = mall;
mg2 = mall;
mh1 = mall;
mh2 = mall;

lg1 = mg1+1; % paper definition of memory: mg1 = lg1-1 
lg2 = mg2+1;
lh1 = mh1+1;
lh2 = mh2+1;

Nt = 30;
Ne = 3000;
Nv = 1000;

noiselevelrel=0.05; % relative noiselevel
plotstuff=0;

Nruns=50; % to plot results of 10 runs on same dataset

test_with_good_initial_approximation=0; deviation_from_exact_init = 1e-1;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
%% Define true pWH system
% FIR filters
ag11 = -0.8 + 1.6*rand;
ag12 = -0.8 + 1.6*rand;
ag13 = -0.8 + 1.6*rand;
g1 = ag11.^[0:lg1-1]+ag12.^[0:lg1-1]+ag13.^[0:lg1-1]; g1=g1(:)./g1(1); 
ah11 = -0.8 + 1.6*rand;
h1 = ah11.^[0:lh1-1]; h1 = h1(:)./h1(1);

ag21 = -0.8 + 1.6*rand;
g2 = ag21.^[0:lg2-1]; g2 = g2(:)./g2(1);

ah21 = -0.8 + 1.6*rand;
h2 = ah21.^[0:lh2-1]; h2 = h2(:)./h2(1);

% static NL poly coeffs
c11 = randn; % coef ^2, branch 1
c12 = randn; % coef ^3, branch 1
c21 = randn;
c22 = randn;

% %%% covar matrix experiment
% mem=max([lg1+lh1 lg2+lh2])-1;
% nbmons2 = nb_mons_partial(mem,2); 
% nbmons3 = nb_mons_partial(mem,3); 
% Nexp = 250;
% Hmatrix = zeros(nbmons2+nbmons3,Nexp);
% for kvar = 1:Nexp
% %%% /covar matrix experiment

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
% Generate pWH i/o data
N  = Ne+Nv;

u = 0.7*randn(2*Nt+N,1); % scaled to avoid very high amplitudes; one transient at beginning, one between est and val data

yg1 = filter(g1,1,u);
uh1 = c11*yg1.^2 + c12*yg1.^3;
yh1 = filter(h1,1,uh1); 

yg2 = filter(g2,1,u);
uh2 = c21*yg2.^2 + c22*yg2.^3;
yh2 = filter(h2,1,uh2);

y0 = yh1+yh2;

ny = noiselevelrel*std(y0)*randn(size(y0));
y=y0+ny;

% SNR in dB
SNRdb = db(rms(y)/rms(ny))

% estimation data set
ue =   u(1:Nt+Ne);
y0e = y0(1:Nt+Ne);
ye =   y(1:Nt+Ne);

% validation data set
uv =   u(2*Nt+Ne+1:end);
y0v = y0(2*Nt+Ne+1:end);
yv =   y(2*Nt+Ne+1:end);


% % plot data
% figure;
%     plotrange=[100:300];
%     subplot(2,1,1);
%     hold all;
%     plot(plotrange,u(plotrange),'x:');
%     legend('u');
%     subplot(2,1,2);
%     hold all;
%     plot(plotrange,y(plotrange),'*:');
%     plot(plotrange,y0(plotrange),'o-');
%     axis([min(plotrange) max(plotrange) -3 3])
%     legend('y=y0+ny', 'y0');
%  
    

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
%% Building data matrix containing products of shifted input variables 

% figuring out positions of indices of Volterra coeffs vector to tensor 
%   (this is useful to know which shifted inputs are to be used). 

% third order Volterra kernel 
mem = max([lg1+lh1 lg2+lh2])-1; % memory length
if lg1+lh1-lg2-lh2 ~=0,
    error('filter lengths incorrect: code requires lg1+lh1=lg2+lh2');
end

% unique elements (wrt symmetry)
K2uniq = zeros(mem,mem);
for i = 1:mem,
    for j = i:mem,
        K2uniq(j,i) = 1;
    end
end

K3uniq = zeros(mem, mem, mem);
for i = 1:mem
    for j = i:mem
        for k = j:mem
            K3uniq(k, j, i) = 1;
        end
    end
end

idxd2 = find(K2uniq > 0); % indices of unique elements
idxd3 = find(K3uniq > 0);

% % build mask matrix/tensor
% K2mask = zeros(mem, mem);
% K3mask = zeros(mem, mem, mem);

nbmons2 = nb_mons_partial(mem,2); 
nbmons3 = nb_mons_partial(mem,3); 

% K2mask(idxd2) = 1:nbmons2;
% K3mask(idxd3) = 1:nbmons3;

% build data matrix A containing powers of shifted inputs
% top: most recent (high index) data
% vector b contains noisy outputs y
A2 = zeros(Ne-mem+1,nbmons2);
A3 = zeros(Ne-mem+1,nbmons3);

b = zeros(Ne-mem+1,1);

% for rowi = 1:Ne-mem+1,
%     nn=Ne+Nt-rowi+1; % added Nt for bottom rows of A needing samples from >N
%     for i = 1:nbmons2,
%         [s21, s22] = find(K2mask == i);
%         A2(rowi,i) = u(nn-s21+1)*u(nn-s22+1);
%     end
%     for i = 1:nbmons3,
%         [s31,s32,s33] = ind2sub([mem mem mem], find(K3mask==i));
%         A3(rowi,i) = u(nn-s31+1)*u(nn-s32+1)*u(nn-s33+1); 
%     end
%     b(rowi) = y(nn); 
% end
for rowi = 1:Ne-mem+1,
    nn=Ne+Nt-rowi+1; % added Nt for bottom rows of A needing samples from >N
    u_temp = flipud(u(nn-mem+1:nn));

    tmp2 = kron(u_temp, u_temp);
    A2(rowi,:) = tmp2(idxd2);

    tmp3 = kron(tmp2, u_temp);
    A3(rowi,:) = tmp3(idxd3);

    b(rowi) = y(nn);
end


hest = [A2 A3]\b;
norm([A2 A3]*hest-b)

hest2 = hest(1:nbmons2);
hest3 = hest(nbmons2+1:end);

% %%% covar matrix exp
% Hmatrix(:,kvar) = hest;
% end
% Hnormmatrix = Hmatrix-repmat(mean(Hmatrix,2),[1 Nexp]);
% figure;
%     surf(Hnormmatrix*Hnormmatrix')
% %%% /covar matrix exp
 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
%% reconstructing system parameters from estimated kernels
% putting estimated Volterra coeffs into 3rd order tensor
K3p1 = zeros(mem, mem, mem);
K3p1(idxd3) = hest3;

K3p2 = permute(K3p1, [1, 3, 2]);
K3p3 = permute(K3p1, [2, 1, 3]);
K3p4 = permute(K3p1, [2, 3, 1]);
K3p5 = permute(K3p1, [3, 1, 2]);
K3p6 = permute(K3p1, [3, 2, 1]);

K3 = K3p1 + K3p2 + K3p3 + K3p4 + K3p5 + K3p6;
K3 = K3/6;

K2p1 = zeros(mem, mem); 
K2p1(idxd2) = hest2;

K2p2 = permute(K2p1, [2 1]); 

K2 = K2p1 + K2p2; K2 = K2/2;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % for loop Nruns
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G1est = zeros(mall+1,Nruns);
G2est = zeros(mall+1,Nruns);
H1est = zeros(mall+1,Nruns);
H2est = zeros(mall+1,Nruns);
Cest = zeros(4,Nruns); % c11 c12 c21 c22
Y0est = zeros(N+2*Nt,Nruns);
E0val = zeros(1,Nruns);
Eval = zeros(1,Nruns);

for runi = 1:Nruns
% solving structured/simultaneous CPD using Structured Data Fusion (Tensorlab 3.0)
model  =  struct;

% parameters
if test_with_good_initial_approximation
    disp('goodapprox used')
    model.variables.a1 = g1(2:end) + deviation_from_exact_init * randn(lg1-1, 1);
    model.variables.a2 = g2(2:end) + deviation_from_exact_init * randn(lg2-1,1 );
    model.variables.b = [1, h1(2:end)' + deviation_from_exact_init * randn(1, lh1 -1 ), ...
                               1 h2(2:end)' + deviation_from_exact_init * randn(1, lh2 -1 )];
    model.variables.c11 = c11 + deviation_from_exact_init * randn;
    model.variables.c12 = c12 + deviation_from_exact_init * randn;
    model.variables.c21 = c21 + deviation_from_exact_init * randn;
    model.variables.c22 = c22 + deviation_from_exact_init * randn;
else
    model.variables.a1 = randn(lg1 -1, 1);
    model.variables.a2 = randn(lg2 -1, 1);
    model.variables.b = [1, randn(1, lh1-1), 1, randn(1, lh2-1)];
    model.variables.c11 = randn;
    model.variables.c12 = randn;
    model.variables.c21 = randn;
    model.variables.c22 = randn;
end

% structure of the factors
toepl1 = @(z,task) struct_toeplitz (z,task,[lg1+lh1-1 lh1],  ...
    [zeros(lh1-1,1); 1], zeros(lh1-1,1));

toepl2 = @(z,task) struct_toeplitz (z,task,[lg2+lh2-1 lh2],  ...
    [zeros(lh2-1,1); 1], zeros(lh2-1,1));

repscalar1 =  @(z,task)  struct_matvec(z,  task,  [],  ones(1, lh1));
repscalar2 =  @(z,task)  struct_matvec(z,  task,  [],  ones(1, lh2));
% MI: there should be a way to combine these two

mask = true(1, lh1+lh2);
mask(1) = false;
mask(lh1+1) = false;

% factors/factorization
model.factors.A  = {{'a1', toepl1},{'a2', toepl2}};
model.factors.B  =  {'b',  @(z,task) struct_const(z,task,mask)};
model.factors.C1 = {{'c12', repscalar1}, {'c22', repscalar1}};
model.factors.C2 = {{'c11', repscalar2}, {'c21', repscalar2}};

model.factorizations.tensor.data = K3;
model.factorizations.tensor.cpd = {'A', 'A', 'A', 'B', 'C1'};
%if test_with_second_order_kernel
    model.factorizations.matrix.data = K2;
    model.factorizations.matrix.cpd = {'A', 'A', 'B', 'C2'};
%end

sdf_check(model, 'print');

% tolerances and solution
options.Display = 25;
options.TolFun = 1e-20;
options.TolX =  1e-20;
options.MaxIter = 500;

tic
%[sol, output] = sdf_nls(model, options); 
[sol, output] = sdf_minf(model, options);
toc


% Output
g1est = [1 sol.variables.a1'];
g2est = [1 sol.variables.a2'];
h1est = sol.variables.b(1:lh1);% / sol.variables.b(1,1);
h2est = sol.variables.b(lh1+1:end);% / sol.variables.b(1,lh1+1);

c11est = sol.variables.c11; % * sol.variables.b(1,1);
c12est = sol.variables.c12; % * sol.variables.b(1,1);
c21est = sol.variables.c21; % * sol.variables.b(1,lh1+1);
c22est = sol.variables.c22; % * sol.variables.b(1,lh1+1);

if output.relerr > 1e-6
    warning('relerr=%f',output.relerr)
end

% simulate (validation) dataset with estimated system
yg1e = filter(g1est,1,u);
uh1e = c11est*yg1e.^2 + c12est*yg1e.^3;
yh1e = filter(h1est,1,uh1e); 

yg2e = filter(g2est,1,u);
uh2e = c21est*yg2e.^2 + c22est*yg2e.^3;
yh2e = filter(h2est,1,uh2e);

y0e = yh1e+yh2e;

if plotstuff,
figure;
hold all; 
plot(y,'x');
plot(y0,'o:');
plot(y0e,'*-')
legend('y','y0','y0e');
end%ifplot

% report rse on validation data 
rmse0val = rms(y0(2*Nt+Ne+1:end)-y0e(2*Nt+Ne+1:end))
rrmsyy0val = rms(y0(2*Nt+Ne+1:end)-y(2*Nt+Ne+1:end))
rmseval = rms(y(2*Nt+Ne+1:end)-y0e(2*Nt+Ne+1:end))

% % save data
G1est(:,runi) = g1est(:);
G2est(:,runi) = g2est(:);
H1est(:,runi) = h1est(:);
H2est(:,runi) = h2est(:);
Cest(:,runi)  = [c11est; c12est; c21est; c22est];
Y0est(:,runi) = y0e;
E0val(runi) = rmse0val; 
Eval(runi)  = rmseval; 

end % for expi

% plot all results
% look at good one
goodones = find(E0val < 0.2);
%goodones = 1:length(E0val);

figure
subplot(2,3,[1, 2])
hold on 
plot(1:lg1, g1,'ko-');
plot(1:lg2, g2,'ko-');
for runi=1:length(goodones),
plot(1:lg1, G1est(:,goodones(runi)),'r.-');
plot(1:lg2, G2est(:,goodones(runi)),'m.-');
end
ylim([-1 1])
title('estimates of g')
legend('true','estimated')
subplot(2,3,[4, 5])
hold on 
plot(1:lh1, h1, 'ko-');
plot(1:lh2, h2, 'ko-');
for runi=1:length(goodones),
plot(1:lh1, H1est(:,goodones(runi)),'r.-');
plot(1:lh2, H2est(:,goodones(runi)),'m.-');
end
ylim([-1 1])
title('estimates of h')
subplot(2,3,[3, 6])
hold on 
plot([2 3], [c11 c12], 'ko-');
plot([2 3], [c21 c22], 'ko-');
for runi=1:length(goodones),
plot([2 3 2 3], Cest(:,goodones(runi)), 'r.');
plot([2 3 2 3], Cest(:,goodones(runi)), 'm.');
end
title('polynomial coefficents')
xlabel('degree')
xlim([1 4]);




figure; 
hold all;
plot(Eval,'.');
plot(E0val,'.');
legend('|y-yhat|^2','|y0-yhat|^2')

datestring=datestr(clock,30)
save([datestring 'Volterra_pWH_experiment1.mat'])


% look at one particular solution
xi=goodones(1);
[g1 g2 G1est(:,xi) G2est(:,xi)]

