lpv = {...
    @(p)0         @(p)1    @(p)0;
    @(p)-0.5         @(p)p    @(p)-0.4;
    @(p)-1         @(p)1    @(p)0
}
% sampling intervals for each parameter
domain = [1.3881,1.6236]
% grid size: number of grid points for each parameter
gridsize = [137]
%% TP transformation, same as:
%   [S U] = tptrans(lpv, dep, domain, gridsize, 'close');

% number of states (size of the A matrix)
n = 2;
% parameter dependencies:
% dep(i,j,k) is 1 if Sp{i,j} depends on p(k)
dep = zeros([size(lpv) 1]);
dep(2,2,:) = [1];


% sampling
lpvdata = sampling_lpv(lpv, dep, domain, gridsize);

% hosvd
[S U sv tol] = hosvd_lpv(lpvdata, dep, gridsize, 0.001);

% generating tight polytopic representation
hull = 'cno';
U = genhull(U, hull);
S = coretensor(U, lpvdata, dep)