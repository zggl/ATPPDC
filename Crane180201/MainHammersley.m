Datanum=200;
% POINT_LIST = DIST_HMERSLEY([-1,1;-1,1], Datanum);
% plot(POINT_LIST(:,1),POINT_LIST(:,2),'r*')
% hold on
p = primes(410); %  prime base
pbaseindx=1;%  prime base vector index
p=p(pbaseindx:end); % p_1=5
POINT_LIST1= DIST_HAMMERSLEYDiy([-1,1;-1,1], Datanum,p);
plot(POINT_LIST1(:,1),POINT_LIST1(:,2),'bo')
pstr2=strcat('p_1= ',num2str(p(1)));
legend('p_1 = 2',pstr2)
dfnamestr=strcat('HammersleyData',num2str(p(1)),'.mat')
save(dfnamestr,'POINT_LIST1')
name='hammersley2d200Nham';
epsname1=strcat(name,'.eps' );
saveas(gcf,epsname1,'epsc2')
pause(2)
close all
