% LPV model
% Mass-Spring-Dampler Mechanical systems
clc,clear
lpv = {@(p)-0.1*p(1)^2,@(p)-0.02-0.67*p(2)^2;
    @(p)1,@(p)0
        };
% number of states (size of the A matrix)
n = 2;
% parameter dependencies:
% dep(i,j,k) is 1 if Sp{i,j} depends on p(k)
dep = zeros([size(lpv) 2]);
dep(1,1,:) = [1 0];
dep(1,2,:) = [1 1];
% sampling intervals for each parameter
% domain = [-45/180*pi 45/180*pi; -0.5 0.5];
a=1.5;b=1.5;
domain = [-a,a;-b,b];
% grid size: number of grid points for each parameter
% gridsize = [23 23];
gridsize = [400,400];
%% TP transformation, same as:
%   [S U] = tptrans(lpv, dep, domain, gridsize, 'close');

% sampling
lpvdata = sampling_lpv(lpv, dep, domain, gridsize);
% hosvd
[S U sv tol] = hosvd_lpv(lpvdata, dep, gridsize, 0.05);

% generating tight polytopic representation
hull = 'cno';
U = genhull(U, hull);
S = coretensor(U, lpvdata, dep);

% plot the results
;

% check model approximation error
[maxerr meanerr] = tperror(lpv, S, U, domain, 100);
disp('max and mean error:'); disp(maxerr); disp(meanerr);

save('tora_data', 'S', 'U', 'n', 'domain', 'gridsize');

