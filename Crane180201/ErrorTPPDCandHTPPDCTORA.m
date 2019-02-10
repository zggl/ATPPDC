%% TPPDC和HTPPDC的网格误差比较
% 2018-01-29
epsilon = 0.2;
% functions in the system matrix
F = @(p) epsilon * p(2) * sin(p(1));
G = @(p) epsilon * cos(p(1));
H = @(p) 1 - epsilon^2 * cos(p(1))^2;
% system matrix: lpv = [A(p) B(p)]
lpv = {...
    @(p)0         @(p)1    @(p)0    @(p)0               @(p)0; 
    @(p)-1/H(p)   @(p)0    @(p)0    @(p)F(p)/H(p)       @(p)-G(p)/H(p);
    @(p)0         @(p)0    @(p)0    @(p)1               @(p)0;
    @(p)G(p)/H(p) @(p)0    @(p)0    @(p)-F(p)*G(p)/H(p) @(p)1/H(p);
};
% number of states (size of the A matrix)
n = 4;
% parameter dependencies:
% dep(i,j,k) is 1 if Sp{i,j} depends on p(k)
dep = zeros([size(lpv) 2]);
dep(2,1,:) = [1 0];
dep(2,4,:) = [1 1];
dep(2,5,:) = [1 0];
dep(4,1,:) = [1 0];
dep(4,4,:) = [1 1];
dep(4,5,:) = [1 0];

% sampling intervals for each parameter
domain = [-45/180*pi 45/180*pi; -45/180*pi 45/180*pi];
% grid size: number of grid points for each parameter
ParaStruct={};
ModelingError=[];
errorind=0;
for ii=50:10:150
    errorind=errorind+1;
    gridsize1 = [ii,ii];    
    %% TP transformation
    U1 = genhull(U1, hull);
    S1 = coretensor(U1, lpvdata1, dep1);
    [S1 U1] = hosvd_lpv(lpvdata1, dep1, gridsize1,svtol, keep);
    % check model approximation error
    [maxerrt,meanerrt] = tperror(lpvux, S1, U1, domain1, 1000);
    %% Hammersley based TP transformation
    svtol=1E-3;
    keep=[3,2];   % reserved sigular value
    primebase=410;%  prime base
    pbaseindx=3;        %  prime base vector index
    hull = 'cno';
    [S1,U1,maxerrh meanerrh]=TPtransformaiton(gridsize1,lpvux,dep1,domain1,hull,svtol, keep,primebase,pbaseindx);
    ModelingError(errorind,:)=[ii,ii,maxerrt,meanerrt,maxerrh meanerrh];
end
save('ModelingErrorTORA.mat','ModelingErrorTORA')
load('ModelingErrorTORA.mat', 'ModelingErrorTORA')
ModelingErrorTORA
