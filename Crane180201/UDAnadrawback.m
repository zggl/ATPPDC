clc
clear all
coli=1;
n1=100;
s=2;
lp1=-1;
rp1=1;
if s==2
    min_ranges_p=[lp1,lp1]  
    max_ranges_p=[rp1, rp1]  
elseif s==3
    min_ranges_p=[lp1,lp1,lp1]  
    max_ranges_p=[rp1, rp1,rp1]  
end
udflag=0;
[X_scaled1,Xij]=UniformDesignWithScale(n1,s,coli,min_ranges_p,max_ranges_p,udflag);
plot(X_scaled1(:,1),X_scaled1(:,2),'r*') % n1=50
hold on
n2=200;
s=2;
[X_scaled2]=UniformDesignWithScale(n2,s,coli,min_ranges_p,max_ranges_p,udflag);
plot(X_scaled2(:,1),X_scaled2(:,2),'bo') % n1=100
legend({'$n=50$','$n=100$'}, 'Interpreter','latex' )
xlabel('$X_1$','Interpreter', 'latex');
ylabel('$X_2$','Interpreter', 'latex');
epsname = strcat('UDdistributePoints', '.eps' );
saveas(gcf,epsname,'epsc2')
pause(3)
epsname = strcat('UDdistributePoints.png');
saveas(gcf,epsname)
close all
%% plot for comparison
epsilon=0.2;
chebfnew = chebfun2(@(x,y) -epsilon.*y*cos(x)./(1 - epsilon^2 .* cos(x)^2),[lp1,rp1,lp1,rp1]);
[ignored,chebfextremaxy] = minandmax2(chebfnew);
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
Fun1=@(x,y) -0.2.*y.*cos(x)./(1 - epsilon.^2 .* cos(x).^2);
fmesh(Fun1,[-1 1 -1 1],'Parent',axes1);
%view(axes1,[-48.3 26.8]);
view(axes1,[56.1 15.6]);
grid(axes1,'on');
axis(axes1,'tight');
hold on 
scatter3(chebfextremaxy(:,1),chebfextremaxy(:,2),[0;0],'ro','filled') % extrema 
hold on 
Funcheb1=Fun1(chebfextremaxy(:,1),chebfextremaxy(:,2));
scatter3(chebfextremaxy(:,1),chebfextremaxy(:,2),Funcheb1,'bo','filled') % extrema 
hold on 
line([chebfextremaxy([1,1],1)],[chebfextremaxy([1,1],2)],[0,Funcheb1(1)],'LineStyle','--')
hold on 
line([chebfextremaxy([2,2],1)],[chebfextremaxy([2,2],2)],[0,Funcheb1(2)],'LineStyle','--')
hold on 
scatter3(X_scaled1(:,1),X_scaled1(:,2),zeros(length(X_scaled1(:,1)),1),'r*') % extrema 
hold on 
scatter3(X_scaled2(:,1),X_scaled2(:,2),zeros(length(X_scaled2(:,1)),1),'bo') % extrema 

Datanum=200;
p = primes(410); %  prime base
pbaseindx=1;%  prime base vector index
p=p(pbaseindx:end); % p_1=5
POINT_LIST1= DIST_HAMMERSLEYDiy([-1,1;-1,1], Datanum,p);
hold on 
scatter3(POINT_LIST1(:,1),POINT_LIST1(:,2),zeros(length(X_scaled2(:,1)),1),'filled') % extrema 

epsname = strcat('UDAnadrawback', '.eps' );
saveas(gcf,epsname,'epsc2')
