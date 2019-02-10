% function ToraSim_NoSimulink()
% LPV model
% reference:
% 	R.T. Bupp, D.S. Bernstein, V.T. Coppola
%	A benchmark problem for nonlinear control design
%	International Journal of Robust and Nonlinear Control, 8:307-310 1998

% state vector:
%	x1: position of the cart
%	x2: velocity of the cart
%	x3: angular position of the proof body
%	x4: angular velocity of the proof body

% parameters of the LPV model:
%	p1: x3 (angular position)
%	p2: x4 (angular velocity)

% constant: coupling coefficient between translational and rotational motion
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
gridsize = [101,101];

%% TP transformation, same as:
%   [S U] = tptrans(lpv, dep, domain, gridsize, 'close');

% sampling
lpvdata = sampling_lpv(lpv, dep, domain, gridsize);

% hosvd
[S U sv tol] = hosvd_lpv(lpvdata, dep, gridsize, 0.001);

% generating tight polytopic representation
hull = 'cno';
U = genhull(U, hull);
S = coretensor(U, lpvdata, dep);

% plot the results
% plothull(U, domain);

% check model approximation error
[maxerr meanerr] = tperror(lpv, S, U, domain, 100);
disp('max and mean error:'); disp(maxerr); disp(meanerr);
% save('tora_data', 'S', 'U', 'n', 'domain', 'gridsize');

% State feedback TP controller design
lmi = lmistruct(S, n);
lmi = lmi_asym_decay(lmi, 0.05);  % |x(t)| < c exp(-0.05 t)
umax = 8;
phi = 1;
lmi = lmi_input(lmi, umax, phi);  % if |x| < phi then |u| < umax
K = lmi_solve(lmi);
x_0=[0.1,0,pi/3,0]';
opt = odeset('RelTol',1e-4,'AbsTol',1e-5);
Tspan=0.01;
SimTime=50;
IterationTimes=floor(SimTime/Tspan);
X=zeros(4,IterationTimes);time=X(1,:);
for i=1:IterationTimes
    time(i)=i*Tspan;
    X(:,i)=x_0;
     % model of TORA
    u = tora_CalU(x_0,K,U,domain);
    [t,x]=ode45(@tora_plant,[0,Tspan],x_0,opt,u,lpv);
    x_0=x(end,:)';
end
%% plot
% figure 1
pos_screen=get(0,'screensize');
fig_pos=linspace(20,pos_screen(3),4);
figure('position',[fig_pos(1),10,600,860])

subplot(2,1,1)
plot(time,X(1,:),'--');
xlabel('time (sec)');ylabel('Position');
% legend({'$q_1$','$q_{d1}$','$\hat{q}_1$'},'Interpreter','latex');

subplot(2,1,2)
plot(time,X(2,:),'--')
xlabel('time (sec)');ylabel('Velocity');
% legend({'$q_2$','$q_{d2}$','$\hat{q}_2$'},'Interpreter','latex');

% figure 2 
figure('position',[fig_pos(2),10,600,860])

subplot(2,1,1), 
plot(time,X(3,:),'-.')
xlabel('time (sec)');ylabel('Angular');
% legend({'$\dot{q}_1$','$\hat{\dot{q}}_1$'},'Interpreter','latex');

subplot(2,1,2), 
plot(time,X(4,:),'-.')
xlabel('time (sec)');ylabel('Angular Velocity');
% legend({'$\dot{q}_2$','$\hat{\dot{q}}_2$'},'Interpreter','latex');


