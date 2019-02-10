% LPV model
tora_lpv

% sampling intervals for each parameter
domain = [-45/180*pi 45/180*pi; -45/180*pi 45/180*pi];
% grid size: number of grid points for each parameter
gridsize = [101 101];

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
plothull(U, domain);

% check model approximation error
[maxerr meanerr] = tperror(lpv, S, U, domain, 100);
disp('max and mean error:'); disp(maxerr); disp(meanerr);

% [XI,YI] = meshgrid(linspace(domain(1,1),domain(1,2),gridsize(1)),linspace(domain(2,1),domain(2,2),gridsize(2)));
% ZI = griddata(XI,YI,XI.*YI,XI,YI);
% [XII,YII ]= meshgrid(U{1}(:,1),U{2}(:,1));
% mesh(XI,YI,XII.*YII)
% [XII,YII ]= meshgrid(U{1}(:,1),U{2}(:,2));
% mesh(XI,YI,XII.*YII)
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(4)/3 50 scrsz(3)/2 scrsz(4)/1.2])
for j=1:5
    for i=1:2
        [XII,YII]= meshgrid(U{1}(:,j),U{2}(:,i));
        subplot(5,2,2*(j-1)+i);
        mesh(XI,YI,XII.*YII)
    end
end
% [XI,YI] = meshgrid(linspace(domain(1,:),gridsize(1)),linspace(domain(2,:),gridsize(2)));
% [XI,YI] = meshgrid(linspace(domain(1,1),domain(1,2),gridsize(1)),linspace(domain(2,1),domain(2,2:),gridsize(2)));
% ZI = griddata(XI,YI,XI.*YI,XI,YI);
