% Reg curves of partition functions
A=6.02*10^23;
kPR=[0;... %1
    0;0;0;0;0;0;... %2-7
    0;0.0140;... %8-9
    0;0;0;0;0;0;... %10-15
    0.0140;... %16
    0;0;0;0;0;0;... %17-22
    0.0140;0;0;... %23-25
    0.0140;0;0;... %26-28
    0;0;0;0;0;0;0;0;0;0;0;0]; %29-40
kPRM=[0;... %1
    0;0;0;0;0;0;... %2-7
    0.00100;0;... %8-9
    0;0;0;0;0;0;... %10-15
    0.00100;... %16
    0;0;0;0;0;0;... %17-22
    0;0.0110;0.00100;... %23-25
    0;0.00100;0.00100;... %26-28
    0;0;0;0;0;0;0;0;... %29-36
    0.0110;0.00100;0.00100;0.0110]; %37-40
kPRE=[0;0.00004;0;0.015];
kPL=[0;0;0;0;0;0.011;0;0;0;0];

x=logspace(-9,-5,100);
% x2=logspace(-10,-5,100);
x2=0;
% PRstate=zeros(numel(x),numel(x));
T=[0,25,30,37];
cRNAP=[0.1,1,10,30,100,200,300];
%%
% Parameter search
for iterT=1:numel(T);
    for iterRNAP=1:numel(cRNAP)

        PRactivity=zeros(numel(x),numel(x),numel(x));
        PRMactivity=zeros(numel(x),numel(x),numel(x));
        PREactivity=zeros(numel(x),numel(x),numel(x));
        PLactivity=zeros(numel(x),numel(x),numel(x));
        for i1=1:numel(x)
            for i2=1:numel(x)
                for i3=1:numel(x2)
                    [a,b,c]=partition_functions(x(i1),x(i2),cRNAP(iterRNAP)*10^-9,x2(i3),T(iterT)+273.15);
                    PRactivity(i1,i2,i3)=sum(a.*kPR);
                    PRMactivity(i1,i2,i3)=sum(a.*kPRM);
                    PREactivity(i1,i2,i3)=sum(b.*kPRE);
                    PLactivity(i1,i2,i3)=sum(c.*kPL);
                end
            end
        end

       
        %%
        maxPR(iterT,iterRNAP)=max(PRactivity(:));
        maxPRM(iterT,iterRNAP)=max(PRMactivity(:));

        PRZeroCIsactivity=PRactivity(1,:,1);
        PR_ZCI_HalfActivity_Cro(iterT,iterRNAP)=x(find(PRZeroCIsactivity < maxPR(iterT,iterRNAP)/2,1));
        PRZeroCrosactivity=PRactivity(:,1,1);
        PR_ZCro_HalfActivity_CI(iterT,iterRNAP)=x(find(PRZeroCrosactivity < maxPRM(iterT,iterRNAP)/2,1));
        PRMZeroCrosactivity=PRMactivity(:,1,1);
        PRM_ZCro_FullActivity_CI(iterT,iterRNAP)=x(find(PRMZeroCrosactivity == max(PRMZeroCrosactivity(:)),1));

        disp(['Max PR activity= ' num2str(maxPR(iterT,iterRNAP))])
        disp(['Max PRM activity= ' num2str(maxPRM(iterT,iterRNAP))])

        disp(['cI=0,50% PR at cro= ' ...
            num2str(PR_ZCI_HalfActivity_Cro(iterT,iterRNAP)*A*10^(-15),4) ' accurate, or ' ...
            num2str(PR_ZCI_HalfActivity_Cro(iterT,iterRNAP)*10^9,4) ' inaccurate'])

        disp(['cro=0,50% PR at cI= ' ...
        num2str(PR_ZCro_HalfActivity_CI(iterT,iterRNAP)*A*10^(-15),4) ' accurate, or ' ...
        num2str(PR_ZCro_HalfActivity_CI(iterT,iterRNAP)*10^9,4) ' inaccurate'])

        disp(['cro=0,100% PRM at cI= ' ...
        num2str(PRM_ZCro_FullActivity_CI(iterT,iterRNAP)*A*10^(-15),4) ' accurate, or ' ...
        num2str(log10(PRM_ZCro_FullActivity_CI(iterT,iterRNAP)),4) ' log10'])
    end
end





%%

% maxPR/0.013-1
% maxPRM/.0063-1
% 
% PR_ZCI_HalfActivity_Cro*10^9/170-1
% PR_ZCro_HalfActivity_CI*10^9/36-1
% log10(PRM_ZCro_FullActivity_CI)/-6.3-1

mse1=(maxPR/0.013-1).^2;
mse2=(maxPRM/.0063-1).^2;

mse3=(PR_ZCI_HalfActivity_Cro*10^9/170-1).^2;
mse4=(PR_ZCro_HalfActivity_CI*10^9/36-1).^2;
mse5=(PRM_ZCro_FullActivity_CI*10^9/500-1).^2;

rmse=sqrt(mse1+mse2+mse3+mse4+mse5);
bestRNAP=cRNAP(min(rmse)==min(rmse(:)))
bestT=T(min(rmse,[],2)==min(rmse(:)))
%%
bestT=37;
bestRNAP=30;
for i1=1:numel(x)
    for i2=1:numel(x)
        for i3=1:numel(x2)
            [a,b,c]=partition_functions(x(i1),x(i2),bestRNAP*10^-9,x2(i3),bestT+273.15);
            PRactivity(i1,i2,i3)=sum(a.*kPR);
            PRMactivity(i1,i2,i3)=sum(a.*kPRM);
            PREactivity(i1,i2,i3)=sum(b.*kPRE);
            PLactivity(i1,i2,i3)=sum(c.*kPL);
        end
    end
end
%%
maxPR2=max(PRactivity(:));
maxPRM2=max(PRMactivity(:));

PRZeroCIsactivity=PRactivity(1,:,1);
PR_ZCI_HalfActivity_Cro2=x(find(PRZeroCIsactivity < maxPR2/2,1));
PRZeroCrosactivity=PRactivity(:,1,1);
PR_ZCro_HalfActivity_CI2=x(find(PRZeroCrosactivity < maxPRM2/2,1));
PRMZeroCrosactivity=PRMactivity(:,1,1);
PRM_ZCro_FullActivity_CI2=x(find(PRMZeroCrosactivity == max(PRMZeroCrosactivity(:)),1));

disp(['Max PR activity= ' num2str(maxPR2)])
disp(['Max PRM activity= ' num2str(maxPRM2)])

disp(['cI=0,50% PR at cro= ' ...
    num2str(PR_ZCI_HalfActivity_Cro2*A*10^(-15),4) ' accurate, or ' ...
    num2str(PR_ZCI_HalfActivity_Cro2*10^9,4) ' inaccurate'])

disp(['cro=0,50% PR at cI= ' ...
num2str(PR_ZCro_HalfActivity_CI2*A*10^(-15),4) ' accurate, or ' ...
num2str(PR_ZCro_HalfActivity_CI2*10^9,4) ' inaccurate'])

disp(['cro=0,100% PRM at cI= ' ...
num2str(PRM_ZCro_FullActivity_CI2*A*10^(-15),4) ' accurate, or ' ...
num2str(log10(PRM_ZCro_FullActivity_CI2),4) ' log10'])

%% Plot Figs
% close all
figure('Name','PR activity')
% surf(log10(x),log10(x),permute(PRactivity(:,:,1),[2,1,3]))%,'EdgeColor','none')
[con1,hc1]=contourf(log10(x),log10(x),permute(PRactivity(:,:,1),[2,1,3]),19);
% alpha(hc1,0.1)
xlabel('log10 (cI_2 (M))')
ylabel('log10 (cro_2 (M))')
zlabel('PR Activity (rxns per second per MOI)')
% grid off
%
h=figure('Name','PRM activity');
% surf(log10(x),log10(x),permute(PRMactivity(:,:,1),[2,1,3]))%,'EdgeColor','none')
[con2,hc2]=contourf(log10(x),log10(x),permute(PRMactivity(:,:,1),[2,1,3]),19);
% alpha(hc2,0.1)
xlabel('log10 (cI_2 (M))')
ylabel('log10 (cro_2 (M))')
zlabel('PRM Activity (rxns per second per MOI)')
% grid off
%%
figure('Name','PRE activity')
plot(log10(x),permute(PREactivity(1,1,:),[1,3,2]))
xlabel('log10 (cII (M))')
ylabel('PRE Activity (rxns per second per MOI)')
% grid off

figure('Name','PRE activity')
% surf(log10(x2),log10(x),permute(PREactivity(1,:,:),[2,3,1]))%,'EdgeColor','none')
contourf(log10(x2),log10(x),permute(PREactivity(1,:,:),[2,3,1]),19)
xlabel('log10 (cII (M))')
ylabel('log10 (cI_2 (M))')
zlabel('PRE Activity (rxns per second per MOI)')
% grid off

figure('Name','PL activity')
% surf(log10(x),log10(x),permute(PLactivity(:,:,1),[2,1,3]))%,'EdgeColor','none')
contourf(log10(x),log10(x),permute(PLactivity(:,:,1),[2,1,3]),19)
xlabel('log10 (cI_2 (M))')
ylabel('log10 (cro_2 (M))')
zlabel('PL Activity (rxns per second per MOI)')
% grid off

%%
% best_score1=1000;
% best_score2=1000;
% best_score3=1000;
% 
% for cRNAP=1:300
%     for T=[25,30,35,37]
%         for i4=1:numel(x)
%             [a,b,c]=partition_functions(0,x(i4),cRNAP*10^-9,x2(i3),T+273.15);
%             PRZeroCIsactivity(i4)=sum(a.*kPR);
% 
% 
%             [a,b,c]=partition_functions(x(i4),0,cRNAP*10^-9,x2(i3),T+273.15);
%             PRZeroCrosactivity(i4)=sum(a.*kPR);        
%         end
%         score1=abs(x(find(PRZeroCIsactivity < maxPRZCI/2,1))*A*10^(-15)-170);
%         if score1 < best_score1
%             best_score1=score1;
%             bestRNAP1=cRNAP;
%         end
%         score2=abs(x(find(PRZeroCrosactivity < maxPRZCro/2,1))*A*10^(-15)-36);
%         if score2 < best_score2
%             best_score2=score2;
%             bestRNAP2=cRNAP;
%         end
%         score3=score1+score2;
%         if score3 < best_score3
%             best_score3=score3;
%             bestRNAP3=cRNAP;
%             bestTemp3=T;
%         end
%     end
% end
% bestRNAP1
% bestRNAP2
% bestRNAP3
% bestTemp3
% %%
% 
% for i4=1:numel(x)
%     [a,b,c]=partition_functions(0,x(i4),bestRNAP3*10^-9,x2(i3),bestTemp3+273.15);
%     PRZeroCIsactivity(i4)=sum(a.*kPR);    
% 
%     [a,b,c]=partition_functions(x(i4),0,bestRNAP3*10^-9,x2(i3),bestTemp3+273.15);
%     PRZeroCrosactivity(i4)=sum(a.*kPR);
% end
% 
% maxPRZCI=max(PRZeroCIsactivity(:));
% maxPRZCro=max(PRZeroCrosactivity(:));
% 
% disp(['cI=0,50% PR at cro=' num2str(...
% x(find(PRZeroCIsactivity < maxPRZCI/2,1))*A*10^(-15))])
% % disp(['cI=0,50% PR at cro=' num2str(...
% % x(find(PRZeroCIsactivity < maxPRZCI/2,1))*10^9)])
% 
% disp(['cro=0,50% PR at cI=' num2str(...
% x(find(PRZeroCrosactivity < maxPRZCro/2,1))*A*10^(-15))])
% % disp(['cro=0,50% PR at cI=' num2str(...
% % x(find(PRZeroCrosactivity < maxPRZCro/2,1))*10^9)])
% 


