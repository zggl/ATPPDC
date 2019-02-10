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
% $B_{eqx} = 14$\si{[Ns/m]} &Equivalent viscous damping coeffcient as seen at the motor pinion\\
% $B_{eq\theta} = 25$\si{[Nms/rad]} &Equivalent viscous damping coeffcient as seen at the jib axis\\
% $B_p = 0.0015$\si{[Nm/s]} &viscous damping coeffcient as seen at the pendulum axis\\
% $B_p = 0.0015$\si{[Nm/s]} &viscous damping coeffcient as seen at the pendulum axis\\
% $\eta_{gx} = 0.55$ &Gearbox effciency\\
% $\eta_{mx} = 0.65$ &Motor effciency\\
% $\eta_{g\theta} = 0.35$ &Gearbox effciency\\
% $\eta_{m\theta} = 0.69$ &Motor effciency\\
% $g = 9.81$\si{[m/s^2]} &Gravitational constant of earth\\
% $J_p= 7.34 \times   10^{-5}$ \si{[kgm^2]} &Load moment of inertia\\
% $J_{arm }= 0.7$\si{[kgm^2]} &Arm moment of inertia\\
% $J_{mx} = 7.32\times  10^{-7}$\si{[kgm^2]} &Rotor moment of inertia\\
% $J_m= 9.44 \times   10^{-5}$\si{[kgm^2]} &Rotor moment of inertia\\
% $K_{gx} = 76.64$ &Planetary gearbox gear ratio\\
% $K_g = 275$ &Planetary gearbox gear ratio\\
% $K_{tx }= 0.032$ &Motor torque constant\\
% $K_t= 0.0195$\si{[Nm/A]} &Motor torque constant\\
% $R_{mx} = 25$\si{[\Omega]} &Motor armature resistance\\
% $R_{m\theta} = 0.6$\si{[\Omega]} &Motor armature resistance\\
% $M = 2.78$\si{[kg]} &Mass of the cart system, including the rotor inertia\\
% $m = 0.32$\si{[kg]} &Load mass\\
% $r_{mp} = 0.0375$\si{[m]}& Motor pinion radius\\
% $B_{eqx} = 14$\si{[Ns/m]} &Equivalent viscous damping coeffcient as seen at the motor pinion\\
% $B_{eq\theta} = 25$\si{[Nms/rad]} &Equivalent viscous damping coeffcient as seen at the jib axis\\
% $B_p = 0.0015$\si{[Nm/s]} &viscous damping coeffcient as seen at the pendulum axis\\
% $B_p = 0.0015$\si{[Nm/s]} &viscous damping coeffcient as seen at the pendulum axis\\
% $\eta_{gx} = 0.55$ &Gearbox effciency\\
% $\eta_{mx} = 0.65$ &Motor effciency\\
% $\eta_{g\theta} = 0.35$ &Gearbox effciency\\
% $\eta_{m\theta} = 0.69$ &Motor effciency\\
% $g = 9.81$\si{[m/s^2]} &Gravitational constant of earth\\
% $J_p= 7.34 \times   10^{-5}$ \si{[kgm^2]} &Load moment of inertia\\
% $J_{arm }= 0.7$\si{[kgm^2]} &Arm moment of inertia\\
% $J_{mx} = 7.32\times  10^{-7}$\si{[kgm^2]} &Rotor moment of inertia\\
% $J_m= 9.44 \times   10^{-5}$\si{[kgm^2]} &Rotor moment of inertia\\
% $K_{gx} = 76.64$ &Planetary gearbox gear ratio\\
% $K_g = 275$ &Planetary gearbox gear ratio\\
% $K_{tx }= 0.032$ &Motor torque constant\\
% $K_t= 0.0195$\si{[Nm/A]} &Motor torque constant\\
% $R_{mx} = 25$\si{[\Omega]} &Motor armature resistance\\
% $R_{m\theta} = 0.6$\si{[\Omega]} &Motor armature resistance\\
% $M = 2.78$\si{[kg]} &Mass of the cart system, including the rotor inertia\\
% $m = 0.32$\si{[kg]} &Load mass\\
% $r_{mp} = 0.0375$\si{[m]}& Motor pinion radius\\
% constant: coupling coefficient between translational and rotational motion
% epsilon = 0.2;
B_eqx = 14
B_eqtheta = 25
B_palpha = 0.0015
B_pphi = 0.0015
eta_gx = 0.55
eta_mx = 0.65
eta_gtheta = 0.35
eta_mtheta = 0.69
g = 9.81
J_b= 7.34E-5
J_arm= 0.7
J_mx = 7.32E-7
J_cm= 9.44E-5
K_gx = 76.64
K_gtheta = 275
K_tx= 0.032
K_ttheta= 0.0195
K_m=0.03
R_mx = 25
R_mtheta = 0.6
M = 2.78
m = 0.32
r_mp = 0.0375
L=0.6
% system matrix: lpv = [A(p) B(p)]
D_1=@(p)((M+m)*J_cm+m*L^2*(M+m*sin(p(1))^2));
D_2=@(q)((J_cm+m*L^2)*(J_b+(m+M)*q(1)^2+eta_gtheta*K_gtheta^2*J_mtheta)-(m*L*L*cos(q(2)))^2);
D_3=@(q)((J_cm+m*L^2)*(J_b+(m+M)*q(1)^2+eta_gtheta*K_gtheta^2*J_mtheta)-(m*L*q(1)*cos(q(2)))^2);
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
% q(t) x(1),x(7),x(8)
% lpvutheta = {...
%    @(q)0         @(q)1    @(q)0    @(q)0               @(q)0;
%    @(q)0 ...
%    @(q)-(J_cm+m*L^2)*(B_eqtheta+(eta_mtheta*eta_gtheta*K_ttheta*K_gtheta^2*K_m/R_mtheta))/D_2(q)...
%    @(q)m^2*g*L^2*q(1)*sinc(q(2)/pi)*cos(q(2))/D_2(q)...
%    @(q)-m*L* q(1)*(J_cm+m*L^2)*sin(q(2))* q(3)/D_2(q)...
%    @(q)((J_cm+m*L^2)*eta_mtheta*eta_gtheta*K_ttheta*K_gtheta)/(D_2(q)*R_mtheta);
%    @(q)0         @(q)0    @(q)0    @(q)1               @(q)0;
%    @(q)0 ...
%    @(q)-m*L*L*(B_eqtheta+(eta_mtheta*eta_gtheta*K_ttheta*K_gtheta^2*K_m)/(R_mtheta)*cos(q(2)))/D_3(q)...
%    @(q)((J_b+(m+M)*q(1)^2+eta_gtheta*K_gtheta^2*J_mtheta)*m*g*L*sinc(q(2)/pi))/D_3(q)...
%    @(q)-(m*L*q(1))^2*sin(q(2))*cos(q(2))* q(3)/D_3(q)...
%    @(q)m*L*L*eta_mtheta*eta_gtheta*K_ttheta*K_gtheta*cos(q(2))/(D_3(q)*R_mtheta)
}
% number of states (size of the A matrix)
n1 = 4;
% parameter dep1endencies:
% dep1(i,j,k) is 1 if Sp{i,j} dep1ends on p(k)
dep1 = zeros([size(lpvux) 2]);
dep1(2,2,:) = [1 0];
dep1(2,3,:) = [1 0]; 
dep1(2,4,:) = [1 1];
dep1(2,5,:) = [1 0];
dep1(4,2,:) = [1 0];
dep1(4,3,:) = [1 0];
dep1(4,4,:) = [1 1];
dep1(4,5,:) = [1 0];

% n2 = 4;
% dep2 = zeros([size(lpvutheta) 3]);
% dep2(2,2,:) = [1 1 0];
% dep2(2,3,:) = [1 1 0]; 
% dep2(2,4,:) = [1 1 1];
% dep2(2,5,:) = [1 1 0];
% dep2(4,2,:) = [1 1 0];
% dep2(4,3,:) = [1 1 0];
% dep2(4,4,:) = [1 1 1];
% dep2(4,5,:) = [1 1 0];

% sampling intervals for each parameter
domain1 = [-5/180*pi 5/180*pi; -0.45 0.45];
% domain2 = [0.24, 0.6;-5/180*pi 5/180*pi; -0.45 0.45];
% grid size: number of grid points for each parameter
gridsize1 = [101,101];
% gridsize2 = [101,101,101];
%% TP transformation, same as:
%   [S U] = tptrans(lpv, dep1, domain, gridsize, 'close');

% sampling
lpvdata1 = sampling_lpv(lpvux, dep1, domain1, gridsize1);
% lpvdata2 = sampling_lpv(lpvutheta, dep2, domain2, gridsize2);
% hosvd
[S1 U1 sv tol] = hosvd_lpv(lpvdata1, dep1, gridsize1, 0.001);
% [S2 U2 sv tol] = hosvd_lpv(lpvdata2, dep2, gridsize2, 0.01);
% generating tight polytopic representation
hull = 'cno';
U1 = genhull(U1, hull);
S1 = coretensor(U1, lpvdata1, dep1);
% U2 = genhull(U2, hull);
% S2 = coretensor(U2, lpvdata2, dep2);
% plot the results
% plothull(U, domain);

% check model approximation error
[maxerr meanerr] = tperror(lpvux, S1, U1, domain1, 100);
disp('max and mean error:'); disp(maxerr); disp(meanerr);
% [maxerr2 meanerr2] = tperror(lpvutheta, S2, U2, domain2, 100);
% disp('max and mean error:'); disp(maxerr2); disp(meanerr2);
% save('tora_data', 'S', 'U', 'n', 'domain', 'gridsize');

% State feedback TP controller design
lmi1 = lmistruct(S1, n1);
lmi1 = lmi_asym_decay(lmi1, 0.05);  % |x(t)| < c exp(-0.05 t)
umax = 8;
phi = 1;
lmi1 = lmi_input(lmi1, umax, phi);  % if |x| < phi then |u| < umax
K1 = lmi_solve(lmi1);

% lmi2 = lmistruct(S2,[ n2,4]);
% lmi2 = lmi_asym_decay(lmi2, 0.05);  % |x(t)| < c exp(-0.05 t)
% umax = 8;
% phi = 1;
% lmi2= lmi_input(lmi2, umax, phi);  % if |x| < phi then |u| < umax
% K2 = lmi_solve(lmi2);

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
    u = towercraneux_CalU(x_0,K1,U1,domain1);
    [t,x]=ode45(@towercraneUx_plant,[0,Tspan],x_0,opt,u,lpvux);
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


