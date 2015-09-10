% Run LambdaModel
%%
clear all %#ok<*CLSCR>
close all
rannum=floor(rand(1)*10000000);
%%
numruns=0;
runnum=1;
%
numruns=numruns+1;
defaultprefix = datestr(date,['yyyy','mm','dd']);
defaultDIR = [pwd '\'];
defaultDIR = 'C:\Users\frink\OneDrive\OneNoteMatlab\Golding\Gillespie\Outputs\';
while runnum<=numruns
    tic
    clear run
    disp(['Run ' num2str(runnum)])
    run(runnum) = LambdaModel(rannum);
    % Save
    save([defaultDIR defaultprefix '_' num2str(runnum,'%.3u')], 'run');
    runnum=runnum+1;
    toc
end


%% 
runnum=1;
load([defaultDIR defaultprefix '_' num2str(runnum,'%.3u')], 'run');
%% Plot
close all
figure('Name','PR')
hold on
subplot(3,1,1)
    plot(run.t/60,run.PRmRNA,'r')
    title('PR mRNA')
    % xlabel('min'); 
    ylabel('molecules')
subplot(3,1,2)
    plot(run.t/60,run.cro_2,'k')
    title('cro_2')
%     xlabel('min'); 
    ylabel('molecules')
subplot(3,1,3)
    plot(run.t/60,run.cII,'k')
    title('cII')
    xlabel('min'); 
    ylabel('molecules')
%
figure('Name','PRM')
hold on
subplot(3,1,1)
    plot(run.t/60,run.PRMmRNA,'r')
    title('PRM mRNA')
    % xlabel('min'); 
    ylabel('molecules')
subplot(3,1,2)
    plot(run.t/60,run.PREmRNA,'r')
    title('PRE mRNA')
%     xlabel('min'); 
    ylabel('molecules')
subplot(3,1,3)
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





