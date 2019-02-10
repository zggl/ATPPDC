clc
clear

load ATPPDC_X.mat
load TPDC_X.mat
load UDTPC_X.mat
load Hammersley_X.mat
% figure 1
Tspan=0.01;
SimTime=50;
IterationTimes=floor(SimTime/Tspan);
for i=1:IterationTimes
    time(i)=i*Tspan;
end
pos_screen=get(0,'screensize');
fig_pos=linspace(20,pos_screen(3),4);
figure('position',[fig_pos(1),10,600,860])
linsty={'--',':','-.','-'};
markesty={'*','o','v','d'};
%% figures
subplot(2,1,1)
indxda=1;
indxstep=50;
data1=[ATPPDC_X(indxda,:);TPDC_X(indxda,:);UDTPC_X(indxda,:);Hammersley_X(indxda,:)];
for iii=1:4
    plot(time, data1(iii,:),'LineStyle',linsty{iii},'Marker',markesty{iii},'MarkerIndices',1:indxstep:size(data1,2),'LineWidth',1) 
    hold on
end
xlabel('time (sec)');ylabel('Position');
legend('ATPDC','TPDC','UDTPDC','HTPDC')
% legend({'$q_1$','$q_{d1}$','$\hat{q}_1$'},'Interpreter','latex');

subplot(2,1,2)
indxda=2;
data1=[ATPPDC_X(indxda,:);TPDC_X(indxda,:);UDTPC_X(indxda,:);Hammersley_X(indxda,:)];
for iii=1:4
    plot(time, data1(iii,:),'LineStyle',linsty{iii},'Marker',markesty{iii},'MarkerIndices',1:indxstep:size(data1,2),'LineWidth',1)  
    hold on
end
xlabel('time (sec)');ylabel('Velocity');
legend('ATPDC','TPDC','UDTPDC','HTPDC')
% legend({'$q_2$','$q_{d2}$','$\hat{q}_2$'},'Interpreter','latex');
epsname = strcat('TORACompare1', '.eps' );
saveas(gcf,epsname,'epsc2')
% figure 2 
figure('position',[fig_pos(2),10,600,860])

subplot(2,1,1), 
indxda=3;
data1=[ATPPDC_X(indxda,:);TPDC_X(indxda,:);UDTPC_X(indxda,:);Hammersley_X(indxda,:)];
for iii=1:4
    plot(time, data1(iii,:),'LineStyle',linsty{iii},'Marker',markesty{iii},'MarkerIndices',1:indxstep:size(data1,2),'LineWidth',1)  
    hold on
end
xlabel('time (sec)');ylabel('Angular');
legend('ATPDC','TPDC','UDTPDC','HTPDC')

subplot(2,1,2), 
indxda=4;
data1=[ATPPDC_X(indxda,:);TPDC_X(indxda,:);UDTPC_X(indxda,:);Hammersley_X(indxda,:)];
for iii=1:4
    plot(time, data1(iii,:),'LineStyle',linsty{iii},'Marker',markesty{iii},'MarkerIndices',1:indxstep:size(data1,2),'LineWidth',1)  
    hold on
end
xlabel('time (sec)');ylabel('Angular Velocity');
legend('ATPDC','TPDC','UDTPDC','HTPDC')
epsname = strcat('TORACompare2','.eps' );
saveas(gcf,epsname,'epsc2')