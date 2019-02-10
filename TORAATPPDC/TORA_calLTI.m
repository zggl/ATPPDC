%  print LTI model
clc,clear
% Model specification: 
% Model constants
% constant: coupling coefficient between translational and rotational motion
epsilon = 0.2;
% functions in the system matrix
F = @(p) epsilon * p(2) * sin(p(1));
G = @(p) epsilon * cos(p(1));
H = @(p) 1 - epsilon^2 * cos(p(1))^2;
% LPV model
% system matrix: lpv = [A(p) B(p)]
lpv = {...
    @(p)0         @(p)1    @(p)0    @(p)0               @(p)0; 
    @(p)-1/H(p)   @(p)0    @(p)0    @(p)F(p)/H(p)       @(p)-G(p)/H(p);
    @(p)0         @(p)0    @(p)0    @(p)1               @(p)0;
    @(p)G(p)/H(p) @(p)0    @(p)0    @(p)-F(p)*G(p)/H(p) @(p)1/H(p);
};
% parameter dependencies
dep = zeros([size(lpv) 2]);
dep(2,1,:) = [1 0];
dep(2,4,:) = [1 1];
dep(2,5,:) = [1 0];
dep(4,1,:) = [1 0];
dep(4,4,:) = [1 1];
dep(4,5,:) = [1 0];
%%%%TP model transformation: 
% TP model transformation:

filename='tora_data';
A_size = 4;

% intervals
domain = [-45/180*pi 45/180*pi; -45/180*pi 45/180*pi];
% grid size
gridsize = [101, 101];
% weight function type: 'cno', 'snnn', 'canonic'
wtype = {'cno' 'cno'};

% TP model transformation
reply = input('TP model transformation? Y/N [Y]: ', 's');
if isempty(reply) || lower(reply)=='y'
	disp('Step 1: Sampling the LPV model');
	tic
    % sampling
    lpvdata = sampling_lpv(lpv, dep, domain, gridsize);
    % hosvd
    [S U sv tol] = hosvd_lpv(lpvdata, dep, gridsize, 0.001);
    % generating tight polytopic representation
    hull = 'cno';
    U = genhull(U, hull);
    S = coretensor(U, lpvdata, dep);
	toc
	disp('L_2 norm and mean square error');
    % check model approximation error
    [maxerr meanerr] = tperror(lpv, S, U, domain, 100);
    disp('max and mean error:'); disp(maxerr); disp(meanerr);
	save(filename, 'S', 'U', 'domain', 'gridsize');
else
	load(filename, 'S', 'U');
end

% plot the results
% plothull(U, domain);

% Plotting
reply = input('Draw weight functions? Y/N [Y]: ', 's');
if isempty(reply) || lower(reply)=='y'
    labels{1,1}='Angular position: x_3 (rad)';
    labels{1,2}='Weighting functions';
    labels{2,1}='Angular speed: x_4 (rad/sec)';
    labels{2,2}='Weighting functions';
	plotw(U, gridsize, domain, labels);
end

reply = input('Print TP model? Y/N [Y]: ', 's');
if isempty(reply) || lower(reply)=='y'
    Lti=printtp(S);
end
%% plot  the two dimensional weighting funciton
[XI,YI] = meshgrid(linspace(domain(1,1),domain(1,2),gridsize(1)),linspace(domain(2,1),domain(2,2),gridsize(2)));
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(4)/3 50 scrsz(3)/2 scrsz(4)/1.2])
for i=1:5
    for j=1:2
        [XII,YII ]= meshgrid(U{1}(:,i),U{2}(:,j));
        subplot(5,2,2*(i-1)+j);
        mesh(XI,YI,abs(XII.*YII));xlabel('x_3'); ylabel('x_4'); 
    end
end