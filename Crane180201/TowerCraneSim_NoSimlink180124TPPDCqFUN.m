function ParaStruct=TPPDCMain()
%% TPPDC ¸ú×Ù¿ØÖÆÆ÷
% 2018-01-26
clc
clear
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
% $B_{eq\theta} = 25$\si{[Nms/rad]} &Equivalent viscous damping coeffcient as seen at the jib axis\\
% $B_p = 0.0015$\si{[Nm/s]} &viscous damping coeffcient as seen at the pendulum axis\\
% $B_p = 0.0015$\si{[Nm/s]} &viscous damping coeffcient as seen at the pendulum axis\\
% $\eta_{g\theta} = 0.35$ &Gearbox effciency\\
% $\eta_{m\theta} = 0.69$ &Motor effciency\\
% $g = 9.81$\si{[m/s^2]} &Gravitational constant of earth\\
% $J_p= 7.34 \times   10^{-5}$ \si{[kgm^2]} &Load moment of inertia\\
% $J_{arm }= 0.7$\si{[kgm^2]} &Arm moment of inertia\\
% $J_m= 9.44 \times   10^{-5}$\si{[kgm^2]} &Rotor moment of inertia\\
% $K_g = 275$ &Planetary gearbox gear ratio\\
% $K_t= 0.0195$\si{[Nm/A]} &Motor torque constant\\
% $R_{m\theta} = 0.6$\si{[\Omega]} &Motor armature resistance\\
% $M = 2.78$\si{[kg]} &Mass of the cart system, including the rotor inertia\\
% $m = 0.32$\si{[kg]} &Load mass\\
% $r_{mp} = 0.0375$\si{[m]}& Motor pinion radius\\
% $B_{eq\theta} = 25$\si{[Nms/rad]} &Equivalent viscous damping coeffcient as seen at the jib axis\\
% $B_p = 0.0015$\si{[Nm/s]} &viscous damping coeffcient as seen at the pendulum axis\\
% $B_p = 0.0015$\si{[Nm/s]} &viscous damping coeffcient as seen at the pendulum axis\\
% $\eta_{g\theta} = 0.35$ &Gearbox effciency\\
% $\eta_{m\theta} = 0.69$ &Motor effciency\\
% $g = 9.81$\si{[m/s^2]} &Gravitational constant of earth\\
% $J_p= 7.34 \times   10^{-5}$ \si{[kgm^2]} &Load moment of inertia\\
% $J_{arm }= 0.7$\si{[kgm^2]} &Arm moment of inertia\\
% $J_m= 9.44 \times   10^{-5}$\si{[kgm^2]} &Rotor moment of inertia\\
% $K_g = 275$ &Planetary gearbox gear ratio\\
% $K_t= 0.0195$\si{[Nm/A]} &Motor torque constant\\
% $R_{m\theta} = 0.6$\si{[\Omega]} &Motor armature resistance\\
% $M = 2.78$\si{[kg]} &Mass of the cart system, including the rotor inertia\\
% $m = 0.32$\si{[kg]} &Load mass\\
% $r_{mp} = 0.0375$\si{[m]}& Motor pinion radius\\
% constant: coupling coefficient between translational and rotational motion
% epsilon = 0.2;

    disflag =1; % display figures 1:display
    distime=15; % display figure time
    simflag=0;  %simflag=1£¬run the sim

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

    % sampling intervals for each parameter
    domain1 = [-5/180*pi 5/180*pi; -0.45 0.45];
    % grid size: number of grid points for each parameter
    ParaStruct={};
    gridsize1 = [101,101];
    %% TPmodel transformaiton
    svtol=1E-3;
    keep=[3,2];   % reserved sigular value
    hull = 'cno';
    % State feedback TP controller design
    %lmi = lmistruct(S1, n1);
    %% LMI solve

    [S1,U1,maxerr,meanerr]=TPtransformaiton(gridsize1,lpvux,dep1,domain1,hull,svtol, keep);
    phi = 1;
    ind=0;
    for alpha=0.01:1E-3:1
        for umax = 0:20
            ind=ind+1;
            [K1,P1,Kim,diagnostic]=LMIsolve(S1,U1,n1,m,alpha,umax,phi);
            data=[umax,alpha,diagnostic.solveroutput.copt]; % save variables
            ParaStruct.feasible{ind,1}=data;
        end
    end
    % feasible value
    % diagnostic.solveroutput.copt  
    if simflag==1 %run simulation
        Kim=squeeze(Kim);
        x_0=[0.3,0,pi/3,0]';
        opt = odeset('RelTol',1e-4,'AbsTol',1e-5);
        Tspan=0.01;
        SimTime=50;
        IterationTimes=floor(SimTime/Tspan);
        X=zeros(4,IterationTimes);
        time=X(1,:);
        %%y_r
        y_r=0.4;
        q=(y_r-x_0(1))*Tspan;% extended signal q

        for i=1:IterationTimes
            time(i)=i*Tspan;
            X(:,i)=x_0;
             % model of cranes
            u = towercraneux_CalUq(x_0,K1,U1,domain1,n1,q);
            [t,x]=ode45(@towercraneUx_plantq,[0,Tspan],x_0,opt,u,lpvux);
            q=q+(y_r-x_0(1))*Tspan;
            x_0=x(end,:)';
        end
    end
    if (disflag ==1) &(simflag==1) 
       %% plot
        % figure 1
        pos_screen=get(0,'screensize');
        fig_pos=linspace(20,pos_screen(3),4);
        figure('position',[fig_pos(1),10,600,860])

        subplot(2,1,1)
        plot(time,X(1,:),'--');
        xlabel('time (sec)');ylabel('Trolley Position');
        % legend({'$q_1$','$q_{d1}$','$\hat{q}_1$'},'Interpreter','latex');

        subplot(2,1,2)
        plot(time,X(2,:),'--')
        xlabel('time (sec)');ylabel('Trolley Velocity');
        % legend({'$q_2$','$q_{d2}$','$\hat{q}_2$'},'Interpreter','latex');

        % figure 2 
        figure('position',[fig_pos(2),10,600,860])

        subplot(2,1,1), 
        plot(time,X(3,:),'-.')
        xlabel('time (sec)');ylabel('Pendulum Swing Angle');
        % legend({'$\dot{q}_1$','$\hat{\dot{q}}_1$'},'Interpreter','latex');

        subplot(2,1,2), 
        plot(time,X(4,:),'-.')
        xlabel('time (sec)');ylabel('Pendulum Swing Angular Velocity');
        % legend({'$\dot{q}_2$','$\hat{\dot{q}}_2$'},'Interpreter','latex');
        pause(distime)
        close all    
    end
    save('ParaStructData.mat','ParaStruct')
    load('ParaStructData.mat', 'ParaStruct')
    datapfea=reshape(cell2mat(ParaStruct.feasible),[],length(ParaStruct.feasible{1}))
    figure(1)
    for iii=1:size(datapfea,1)
        if datapfea(iii,3)<0
            plot(datapfea(iii,1),datapfea(iii,2),'*')
        end
        hold on
    end
    hold off
    xlabel('\phi');
    ylabel('\alpha')
    box on
    name='TPPDFeasibleRegion';
    epsname1=strcat(name,'.eps' );
    saveas(gcf,epsname1,'epsc2')
    pause(distime)
    close all   
    
end
%%%%%%%%%%%%%%%%%%%%%%
function [S1,U1,maxerr,meanerr]=TPtransformaiton(gridsize1,lpvux,dep1,domain1,hull,svtol, keep,primebase,pbaseindx)
%% TP transformation, same as:
    %   [S U] = tptrans(lpv, dep1, domain, gridsize, 'close');
    % sampling
  
    lpvdata1 =sampling_lpv(lpvux, dep1, domain1, gridsize1);
    %lpvdata1 = sampling_lpv(lpvux, dep1, domain1, gridsize1);
    % hosvd
    % [S U sv tol] = hosvd_lpv(data, dep, gridsize, svtol, keep)

    [S1 U1 sv tol] = hosvd_lpv(lpvdata1, dep1, gridsize1, svtol, keep);

    % generating tight polytopic representation
    U1 = genhull(U1, hull);
    S1 = coretensor(U1, lpvdata1, dep1);

    % check model approximation error
    [maxerr meanerr] = tperror(lpvux, S1, U1, domain1, 100);
    % disp('max and mean error:'); 
    % disp(maxerr);
    % disp(meanerr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K1,P1,Kim,diagnostic]=LMIsolve(S1,U1,n1,m,alpha,umax,phi)
    I = size(S1);
    sizes = I(1:end-2);
    R = prod(sizes);

    S1 = reshape(S1, [R I(end-1) I(end)]);

    m = I(end) - n1;
    p = I(end-1) - n1;
    A = S1(:, 1:n1, 1:n1);
    B = S1(:, 1:n1, n1+1:n1+m);
    C = S1(:, n1+1:n1+p, 1:n1);
    D = S1(:, n1+1:n1+p, n1+1:n1+m);

    X = sdpvar(n1+1, n1+1, 'symmetric');
    M = cell(1, R);
    for r = 1:R
        M{r} = sdpvar(m, n1+1, 'full');
    end

    % X > 0
    lmi.F = set(X > 0, 'poz def');

    %lmi = lmi_asym_decay(lmi1, 0.05);  % |x(t)| < c exp(-0.05 t)
    R = size(A, 1);
    % X*Ar' + Ar*X - Br*Mr - Mr'*Br' + 2*alpha*X < 0
    for r = 1:R
        Ar = [reshape(A(r,:,:), [n1 n1]),zeros(n1,1);-[1,0,0,0,0]];
        Br = [reshape(B(r,:,:), [n1 m]);0];
        lmi.F = lmi.F + set(X*Ar' + Ar*X - Br*M{r} - M{r}'*Br' + 2*alpha*X < 0, sprintf('type1 lmi %d', r));
    end

    % X*Ar' + Ar*X + X*As' + As*X - Br*Ms - Ms'*Br' - Bs*Mr - Mr'*Bs' + 4*alpha*X <= 0
    for r = 1:R
        for s = r+1:R
            Ar = [reshape(A(r,:,:), [n1 n1]),zeros(n1,1);-[1,0,0,0,0]];
            As = [reshape(A(s,:,:), [n1 n1]),zeros(n1,1);-[1,0,0,0,0]];
            Br = [reshape(B(r,:,:), [n1 m]);0];
            Bs = [reshape(B(s,:,:), [n1 m]);0];
            lmi.F = lmi.F + set(X*Ar' + Ar*X + X*As' + As*X - Br*M{s} - M{s}'*Br'...
                      - Bs*M{r} - M{r}'*Bs' + 4*alpha*X <= 0, sprintf('type2 lmi %d', r));
        end
    end

    %lmi1 = lmi_input(lmi, umax, phi);  % if |x| < phi then |u| < umax
    R = size(A, 1);
    % constraints on the control value

    % phi^2 I < X
    lmi.F = lmi.F + set(phi^2 * eye(n1+1) < X, 'phi^2 I < X');

    % [X, Mr'; Mr, mu^2 I] > 0
    for r = 1:R
        lmi.F = lmi.F + set([X M{r}'; M{r} umax^2*eye(m)] > 0, sprintf('type3 lmi %d', r));
    end
    lmi.sizes = sizes;
    lmi.n = n1+1;
    lmi.m = m;
    lmi.p = p;
    lmi.A = A;
    lmi.B = B;
    lmi.C = C;
    lmi.D = D;

    lmi.X = X;
    lmi.M = M;

    [K1,P1,Kim,diagnostic] = lmi_solve(lmi);
    % best value of t: 
    %diagnostic.solveroutput.copt
end
