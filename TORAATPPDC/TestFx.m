% LPV model

lpv = {@(p)(1+p(1)^-2+p(2)^-1.5)^2};
% number of states (size of the A matrix)
n = 1;
% parameter dependencies:
% dep(i,j,k) is 1 if Sp{i,j} depends on p(k)
dep = zeros([size(lpv) 2]);
dep(1,1,:) = [1 1];

% sampling intervals for each parameter
% domain = [-45/180*pi 45/180*pi; -0.5 0.5];
domain = [1,5;1,5];
% grid size: number of grid points for each parameter
% gridsize = [23 23];
gridsize = [800,800];
%% TP transformation, same as:
%   [S U] = tptrans(lpv, dep, domain, gridsize, 'close');

% sampling
lpvdata = sampling_lpv(lpv, dep, domain, gridsize);

[U1,D,U2]=svd(lpvdata{1})
U1=U1(:,[1:3])
figure(1)
plot(U1)
U2=U2(:,[1:3])
figure(2)
plot(U2)

% x = rand(1,50);
% y = rand(1,50);
% z = peaks(6*x-3,6*x-3);
% tri = delaunay(x,y);
% trisurf(tri,x,y,z)


% ti = linspace(1,5,5);
ti = linspace(1,5,800);
[XI,YI] = meshgrid(ti,ti);
ZI=(1+XI.^-2+YI.^-1.5).^2;
C=zeros(800,800,3);
C(:,:,1)=0; C(:,:,2)=0;C(:,:,3)=200;
mesh(XI,YI,ZI,C), hold
plot3(XI,YI,ZI,'*'), hold off



% rand('seed',0)
% x = rand(100,1)*4-2;  y = rand(100,1)*4-2;
% z = x.*exp(-x.^2-y.^2);
% ti = -2:.25:2; 
% [XI,YI] = meshgrid(ti,ti);
% ZI = griddata(x,y,z,XI,YI);
% mesh(XI,YI,ZI), hold
% plot3(x,y,z,'o'), hold off




% hosvd
[S U sv tol] = hosvd_lpv(lpvdata, dep, gridsize, 0.05);

% generating tight polytopic representation
hull = 'cno';
U = genhull(U, hull);
S = coretensor(U, lpvdata, dep);

% plot the results
plothull(U, domain);

% check model approximation error
[maxerr meanerr] = tperror(lpv, S, U, domain, 100);
disp('max and mean error:'); disp(maxerr); disp(meanerr);

save('tora_data', 'S', 'U', 'n', 'domain', 'gridsize');