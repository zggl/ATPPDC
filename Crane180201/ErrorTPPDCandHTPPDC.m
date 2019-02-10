%% TPPDC和HTPPDC的网格误差比较
% 2018-01-29
disflag =1; % display figures 1:display
distime=15; % display figure time
simflag=0;  %simflag=1，run the sim
B_eqx = 14;
B_eqtheta = 25;
B_palpha = 0.0015;
B_pphi = 1.5E-3;
eta_gx = 0.55;
eta_mx = 0.65;
eta_gtheta = 0.35;
eta_mtheta = 0.69;
g = 9.81;
J_b= 7.34E-5;
J_arm= 0.7;
J_mx = 7.32E-7;
J_cm= 9.44E-5;
K_gx = 76.64;
K_gtheta = 275;
K_tx= 0.032;
K_ttheta= 0.0195;
K_m=0.03;
R_mx = 25;
R_mtheta = 0.6;
M = 2.78;
m = 0.32;
r_mp = 0.0375;
L=0.6;
% system matrix: lpv = [A(p) B(p)]
D_1=@(p)((M+m)*J_cm+m*L^2*(M+m*sin(p(1))^2));
% p(t): x(3),x(4)
lpvux = {...
    @(p)0         @(p)1    @(p)0    @(p)0               @(p)0; 
    @(p)0 ...
    @(p)-(J_cm+m*L^2*B_eqx+eta_gx*K_gx^2 *eta_mx*K_tx*K_m*R_mx*r_mp^2) /D_1(p)...
    @(p)m^2*L^2*g*sinc(p(1)/pi)*cos(p(1))/D_1(p)...
    @(p)((m^2*L^3+L*m*J_cm)*sin(p(1))* p(2)+m*L*cos(p(1))*B_palpha)/D_1(p)...
    @(p)(-(m*L^2+J_cm)^2*eta_gx*K_gx*eta_mx*K_tx)/(D_1(p) *R_mx*r_mp^2);
   @(p)0         @(p)0    @(p)0    @(p)1               @(p)0;
   @(p)0 ...
   @(p)(m*L*cos(p(1))*B_eqx+eta_gx*K_gx^2 *eta_mx*K_tx*K_m*R_mx*r_mp^2)/D_1(p)...
   @(p)-(M+m)*m *g*L*sinc(p(1))/D_1(p)...
   @(p)-(M+m)*B_palpha-m^2*L^2*sin(p(1))*cos(p(1))* p(2)/D_1(p)...
   @(p)(-(m*L*cos(p(1)))^2*eta_gx*K_gx* eta_mx*K_tx)/(D_1(p)* R_mx*r_mp);
};
% number of states (size of the A matrix)
n1 = 4;
dep1 = zeros([size(lpvux) 2]);
dep1(2,2,:) = [1 0];
dep1(2,3,:) = [1 0]; 
dep1(2,4,:) = [1 1];
dep1(2,5,:) = [1 0];
dep1(4,2,:) = [1 0];
dep1(4,3,:) = [1 0];
dep1(4,4,:) = [1 1];
dep1(4,5,:) = [1 0];
% sampling intervals for each parameter
domain1 = [-5/180*pi 5/180*pi; -0.45 0.45];
% grid size: number of grid points for each parameter
ParaStruct={};
ModelingError=[];
errorind=0;
for ii=100:10:160
    errorind=errorind+1;
    gridsize1 = [ii,ii];
    %% TPmodel transformaiton
    svtol=1E-3;
    keep=[3,2];   % reserved sigular value
    primebase=410;%  prime base
    pbaseindx=3;        %  prime base vector index
    hull = 'cno';

    [S1,U1,maxerrt,meanerrt]=TPtransformaiton(gridsize1,lpvux,dep1,domain1,hull,svtol, keep,primebase,pbaseindx);

    %% Hammersley based TP transformation
    p = primes(1010);    %  prime base
    pbaseindx=3;        %  prime base vector index
    p=p(pbaseindx:end); %  p_1=5 5x2
    TP.siz=prod(gridsize1);
    TP.lpv=lpvux;
    TP.dep=dep1;
    TP.domain =domain1;
    TP.gridsize=gridsize1;
    Datanum=TP.siz;
    Xij= DIST_HAMMERSLEYDiy(domain1, Datanum,p);
    TP.X_scaled=Xij;
    lpvdata1 = sampling_lpvud(TP);
    [S1 U1 sv tol] = hosvd_lpv(lpvdata1, dep1, gridsize1,svtol, keep);
    % generating tight polytopic representation
    hull = 'cno';
    U1 = genhull(U1, hull);
    S1 = coretensor(U1, lpvdata1, dep1);
    % check model approximation error
    [maxerrh meanerrh] = tperror(lpvux, S1, U1, domain1, 1000);
    ModelingError(errorind,:)=[ii,ii,maxerrt,meanerrt,maxerrh meanerrh];
end
save('ModelingError.mat','ModelingError')
load('ModelingError.mat', 'ModelingError')
ModelingError
