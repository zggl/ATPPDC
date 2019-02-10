clear
datafile={'ParaStructDataUD2500.mat','ParaStructDataUD3600.mat','ParaStructDataUD4900.mat','ParaStructDataUD6400.mat','ParaStructDataUD8100.mat'}
linsty={'--',':','-.','-','--'};
markesty={'*','o','v','d','^'};
for jj=1:length(datafile)
    load(datafile{jj}, 'ParaStruct')
    DR=cell2mat(ParaStruct.feasible);
    [feaflag,feaindx]=find(DR(:,3)<0);
    DRfea=DR(feaflag,:);
    DRValue=unique(DR(feaflag,1));
    FeaReaMaxData=zeros(size(DRValue,1),3);
    for ii=1:length(DRValue)
        feaflagii=find(DRValue(ii)==DRfea(:,1));
        tmp=DRfea(feaflagii,2);
        tmp1=DRfea(feaflagii,3);
        [M,I] =max(tmp);
        FeaReaMaxData(ii,:)=[ii,M,tmp1(I,1)];
    end
    plot(FeaReaMaxData(:,1),FeaReaMaxData(:,2),'LineStyle',linsty{jj})
    hold on   
end
hold off
legend('n=2500','n=3600','n=4900','n=6400','n=8100','Location','southeast')
name='UDFeasibleRegionAlg';
epsname1=strcat(name,'.eps' );
saveas(gcf,epsname1,'epsc2')
