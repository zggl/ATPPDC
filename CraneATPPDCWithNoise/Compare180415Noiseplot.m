clc
clear
% TPPDC
% TPPDC_180415Noise
% 
% HTPPDC
% HTPPDC_180415Noise

% ATPPDC
% ATPPDC_180415Noise
clear
clc
% UDTPPDC
% UDTPPDC_180415Noise
dt=1.0;
x_0=[0.3,0,pi/4,0]';
y_r=0.4;

ParaStructU=UDTPPC180415Noise(dt,x_0,y_r);
ParaStructH=HTPPD180415Noise(dt,x_0,y_r);
ParaStructA=ATPPDC180415Noise(dt,x_0,y_r);
ParaStructT=TPPDC180415Noise(dt,x_0,y_r);

load('TPPDC180415Noise.mat')
load('HTPPDC180415Noise.mat')
load('ATPPDC180415Noise.mat')
load('UDTPPDC180415Noise.mat')
Tspan=6E-3;
SimTime=50;IterationTimes=floor(SimTime/Tspan);
for i=1:IterationTimes
    time(i)=i*Tspan;
end

linsty={'--',':','-.','-'};
markesty={'*','o','v','d'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1) 
plot(time,tpdc_u,':',time,htpdc_u,'-.',time,atpdc_u,'--',time, udtpdc_u,'LineWidth',1.2)
legend('TPPDC:u(t)','HTPPDC:u(t)','ATPPDC:u(t)','UDTPPDC:u(t)')
xlabel('Time (s)');
ylabel('u(t)')
name='SPGPDCNoiseCompareu';
epsname1=strcat(name,'.eps' );
saveas(gcf,epsname1,'epsc2')
%set(gca,'XTick',0:SimTime);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
indxstep=500;
indxda=1;
data1=[X_PDC(indxda,:);X_htppdc(indxda,:);X_udtppdc(indxda,:);X_atppdc(indxda,:)];
for iii=1:4
    plot(time, data1(iii,:),'LineStyle',linsty{iii},'Marker',markesty{iii},'MarkerIndices',1:indxstep:size(data1,2),'LineWidth',1)  
    hold on
end
legend('TPPDC','HTPPDC','UDTPPDC','ATPPDC','Location','northeast')
ylabel('x_1(t)')
xlabel('Time (s)');
name='SPGPDCNoiseComparex1';
epsname1=strcat(name,'.eps' );
saveas(gcf,epsname1,'epsc2')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
indxstep=500;
indxda=3;
data1=[X_PDC(indxda,:);X_htppdc(indxda,:);X_udtppdc(indxda,:);X_atppdc(indxda,:)];
for iii=1:4
    plot(time, data1(iii,:),'LineStyle',linsty{iii},'Marker',markesty{iii},'MarkerIndices',1:indxstep:size(data1,2),'LineWidth',1)  
    hold on
end
legend('TPPDC','HTPPDC','UDTPPDC','ATPPDC')
ylabel('x_3(t)')
xlabel('Time (s)');
name='SPGPDCNoiseComparex3';
epsname1=strcat(name,'.eps' );
saveas(gcf,epsname1,'epsc2')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
indxstep=500;
indxda=2;
data1=[X_PDC(indxda,:);X_htppdc(indxda,:);X_udtppdc(indxda,:);X_atppdc(indxda,:)];
for iii=1:4
    plot(time, data1(iii,:),'LineStyle',linsty{iii},'Marker',markesty{iii},'MarkerIndices',1:indxstep:size(data1,2),'LineWidth',1)  
    hold on
end
legend('TPPDC','HTPPDC','UDTPPDC','ATPPDC')
ylabel('x_2(t)')
xlabel('Time (s)');
name='SPGPDCNoiseComparex2';
epsname1=strcat(name,'.eps' );
saveas(gcf,epsname1,'epsc2')
figure(5)
indxstep=500;
indxda=4;
data1=[X_PDC(indxda,:);X_htppdc(indxda,:);X_udtppdc(indxda,:);X_atppdc(indxda,:)];
for iii=1:4
    plot(time, data1(iii,:),'LineStyle',linsty{iii},'Marker',markesty{iii},'MarkerIndices',1:indxstep:size(data1,2),'LineWidth',1)  
    hold on
end
legend('TPPDC','HTPPDC','UDTPPDC','ATPPDC')
ylabel('x_4(t)')
xlabel('Time (s)');
name='SPGPDCNoiseComparex4';
epsname1=strcat(name,'.eps' );
saveas(gcf,epsname1,'epsc2')
% save workspace
pause(35)
close all



