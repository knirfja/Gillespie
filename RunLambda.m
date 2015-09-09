% Run LambdaModel
%%
clear all %#ok<*CLSCR>
close all
%%
defaultprefix = datestr(date,['yyyy','mm','dd']);
defaultDIR = [pwd '\'];
numruns=1;
rannum=floor(rand(1)*10000000);
for runnum=1:numruns
    tic
    clear run
    disp(['Run ' num2str(runnum)])
    run(runnum) = LambdaModel(rannum);        
    toc
end

%% Save
save([defaultDIR defaultprefix '_' num2str(i)], run);
%% Plot
figure('Name','PR')
hold on
subplot(4,1,1)
    plot(run.t/60,run.PRmRNA,'r')
    title('PR mRNA')
    % xlabel('min'); 
    ylabel('molecules')
subplot(4,1,2)
    plot(run.t/60,run.PREmRNA,'r')
    title('PRE mRNA')
%     xlabel('min'); 
    ylabel('molecules')
subplot(4,1,3)
    plot(run.t/60,run.cro_2,'k')
    title('cro_2')
%     xlabel('min'); 
    ylabel('molecules')
subplot(4,1,4)
    plot(run.t/60,run.cII,'k')
    title('cII')
    xlabel('min'); 
    ylabel('molecules')
%
figure('Name','PRM')
hold on
subplot(2,1,1)
    plot(run.t/60,run.PRMmRNA,'r')
    title('PRM mRNA')
    % xlabel('min'); 
    ylabel('molecules')
subplot(2,1,2)
    plot(run.t/60,run.cI_2,'k')
    title('cI_2')
    xlabel('min'); 
    ylabel('molecules')
%
figure('Name','PL')
hold on
subplot(3,1,1)
    plot(run.t/60,run.PLmRNA,'r')
    title('PL mRNA')
    % xlabel('min'); 
    ylabel('molecules')
subplot(3,1,2)
    plot(run.t/60,run.N,'k')
    title('N')
%     xlabel('min'); 
    ylabel('molecules')
subplot(3,1,3)
    plot(run.t/60,run.cIII,'k')
    title('cIII')
    xlabel('min'); 
    ylabel('molecules')





