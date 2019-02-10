clc
clear all
coli=1;
n=2500;
s=2;
%domain1 = [-5/180*pi 5/180*pi; -0.45 0.45];
if s==2
    min_ranges_p=[-5/180*pi -0.45]  
    max_ranges_p=[5/180*pi 0.45]  
elseif s==3
    min_ranges_p=[-1,-1,-1]  
    max_ranges_p=[1, 1,1]  
end
udflag=0;
[X_scaled,Xij]=UniformDesignWithScale(n,s,coli,min_ranges_p,max_ranges_p,udflag);
if udflag==0  
    dfnamestr=strcat('UDData',num2str(n),'.mat')
else
    dfnamestr=strcat('UDDatan',num2str(n),'.mat')
end
save(dfnamestr,'Xij')
plot(X_scaled(:,1),X_scaled(:,2),'r*')