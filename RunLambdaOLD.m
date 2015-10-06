% Run LambdaModel OLD
%%
clear all %#ok<*CLSCR>
close all
rannum=floor(rand(1)*10000000);
A=6.02*10^23;
AV=6.02*10^23*1.5*10^-15;
%%
numruns=1;
plot_trace=1;
runnum=1;
%
% numruns=numruns+1;
defaultprefix = datestr(date,['yyyy','mm','dd']);
% defaultDIR = [pwd '\'];
defaultDIR = 'C:\Users\frink\OneDrive\OneNoteMatlab\Golding\Gillespie\Outputs\';
while runnum<=numruns
    tic
    clear run
    disp(['Run ' num2str(runnum)])
    rannum=floor(rand(1)*10000000);
    run = LambdaModel(rannum);
    % Save
    save([defaultDIR defaultprefix '_' num2str(runnum,'%.3u')], 'run');
    runnum=runnum+1;
    toc
end

%% Load all runs
% Be careful! each run takes ~100 Mb RAM. Don't exceed your capacity.
clear runs
numruns=1;
for runnum=39
    clear run
    disp(['Run ' num2str(runnum)])
    load([defaultDIR defaultprefix '_' num2str(runnum,'%.3u')], 'run');
    runs(runnum) = run;
end
clear run

%% Test for size of time step
figure; hist(log(diff(runs(runnum).t)),100)
%% Bin each variable
nbins=3000; % 1 sec bins for 50 min traces
t=runs(runnum).t./60;
binrange = linspace(min(t),max(t),nbins);
[n,bin] = histc(t,binrange);
% n = histc(t,binrange);
for iter1 = 1:nbins
    ind=bin==iter1;
    runsbin(runnum).t(iter1)=mean(runs.t(ind));
    runsbin(runnum).V(iter1)=mean(runs.V(ind));
    runsbin(runnum).cI(iter1)=mean(runs.cI(ind));
    runsbin(runnum).cro(iter1)=mean(runs.cro(ind));    
    runsbin(runnum).cI_2(iter1)=mean(runs.cI_2(ind));
    runsbin(runnum).cro_2(iter1)=mean(runs.cro_2(ind));
    runsbin(runnum).N(iter1)=mean(runs.N(ind));
    runsbin(runnum).cII(iter1)=mean(runs.cII(ind));
    runsbin(runnum).cIII(iter1)=mean(runs.cIII(ind));
    runsbin(runnum).PRmRNA(iter1)=mean(runs.PRmRNA(ind));
    runsbin(runnum).PREmRNA(iter1)=mean(runs.PREmRNA(ind));
    runsbin(runnum).PRMmRNA(iter1)=mean(runs.PRMmRNA(ind));
    runsbin(runnum).PLmRNA(iter1)=mean(runs.PLmRNA(ind));
    runsbin(runnum).randseed=runs.randseed;
end
%% Plot binned
figure('Name',['Run_' num2str(runnum) '_Binned'])
subplot(3,3,1)
    plot(runsbin(runnum).t/60,runsbin(runnum).PRmRNA,'r')
    title('PR mRNA')
    % xlabel('min'); 
    ylabel('molecules')
    ylim([0 10])
subplot(3,3,4)
    plot(runsbin(runnum).t/60,runsbin(runnum).cro_2,'k')
    title('cro_2')
%     xlabel('min'); 
    ylabel('molecules')
    ylim([0 200])
subplot(3,3,7)
    plot(runsbin(runnum).t/60,runsbin(runnum).cII,'k')
    title('cII')
    xlabel('min'); 
    ylabel('molecules')
    ylim([0 150])

subplot(3,3,2)
    plot(runsbin(runnum).t/60,runsbin(runnum).PREmRNA,'r')
    title('PRE mRNA')
    % xlabel('min'); 
    ylabel('molecules')
    ylim([0 10])
subplot(3,3,5)
    plot(runsbin(runnum).t/60,runsbin(runnum).PRMmRNA,'r')
    title('PRM mRNA')
    % xlabel('min'); 
    ylabel('molecules')
    ylim([0 10])
subplot(3,3,8)
    plot(runsbin(runnum).t/60,runsbin(runnum).cI_2,'k')
    title('cI_2')
    xlabel('min'); 
    ylabel('molecules')
    ylim([0 100])
    
subplot(3,3,3)
    plot(runsbin(runnum).t/60,runsbin(runnum).PLmRNA,'r')
    title('PL mRNA')
    % xlabel('min'); 
    ylabel('molecules')
    ylim([0 10])
subplot(3,3,6)
    plot(runsbin(runnum).t/60,runsbin(runnum).N,'k')
    title('N')
%     xlabel('min'); 
    ylabel('molecules')
    ylim([0 60])
subplot(3,3,9)
    plot(runsbin(runnum).t/60,runsbin(runnum).cIII,'k')
    title('cIII')
    xlabel('min'); 
    ylabel('molecules')
    ylim([0 150])
    
% plot(1:nbins,runsbin(runnum).t)



%% sample each variable
nsamp=3000; % 1 sec bins for 50 min traces
t=runs(runnum).t./60;
% binrange = linspace(min(t),max(t),nbins);
% [n,bin] = histc(t,binrange);
% n = histc(t,binrange);
for iter1 = 1:nsamp
    ind=find(t<iter1,1,'last');
    runssamp(runnum).t(iter1)=mean(runs.t(ind));
    runssamp(runnum).V(iter1)=mean(runs.V(ind));
    runssamp(runnum).cI(iter1)=mean(runs.cI(ind));
    runssamp(runnum).cro(iter1)=mean(runs.cro(ind));    
    runssamp(runnum).cI_2(iter1)=mean(runs.cI_2(ind));
    runssamp(runnum).cro_2(iter1)=mean(runs.cro_2(ind));
    runssamp(runnum).N(iter1)=mean(runs.N(ind));
    runssamp(runnum).cII(iter1)=mean(runs.cII(ind));
    runssamp(runnum).cIII(iter1)=mean(runs.cIII(ind));
    runssamp(runnum).PRmRNA(iter1)=mean(runs.PRmRNA(ind));
    runssamp(runnum).PREmRNA(iter1)=mean(runs.PREmRNA(ind));
    runssamp(runnum).PRMmRNA(iter1)=mean(runs.PRMmRNA(ind));
    runssamp(runnum).PLmRNA(iter1)=mean(runs.PLmRNA(ind));
    runssamp(runnum).randseed=runs.randseed;
end

%% Plot Sampled
figure('Name',['Run_' num2str(runnum) '_Sampled'])
subplot(3,3,1)
    plot(runssamp(runnum).t/60,runssamp(runnum).PRmRNA,'r')
    title('PR mRNA')
    % xlabel('min'); 
    ylabel('molecules')
    ylim([0 10])
subplot(3,3,4)
    plot(runssamp(runnum).t/60,runssamp(runnum).cro_2,'k')
    title('cro_2')
%     xlabel('min'); 
    ylabel('molecules')
    ylim([0 200])
subplot(3,3,7)
    plot(runssamp(runnum).t/60,runssamp(runnum).cII,'k')
    title('cII')
    xlabel('min'); 
    ylabel('molecules')
    ylim([0 150])

subplot(3,3,2)
    plot(runssamp(runnum).t/60,runssamp(runnum).PREmRNA,'r')
    title('PRE mRNA')
    % xlabel('min'); 
    ylabel('molecules')
    ylim([0 10])
subplot(3,3,5)
    plot(runssamp(runnum).t/60,runssamp(runnum).PRMmRNA,'r')
    title('PRM mRNA')
    % xlabel('min'); 
    ylabel('molecules')
    ylim([0 10])
subplot(3,3,8)
    plot(runssamp(runnum).t/60,runssamp(runnum).cI_2,'k')
    title('cI_2')
    xlabel('min'); 
    ylabel('molecules')
    ylim([0 100])
    
subplot(3,3,3)
    plot(runssamp(runnum).t/60,runssamp(runnum).PLmRNA,'r')
    title('PL mRNA')
    % xlabel('min'); 
    ylabel('molecules')
    ylim([0 10])
subplot(3,3,6)
    plot(runssamp(runnum).t/60,runssamp(runnum).N,'k')
    title('N')
%     xlabel('min'); 
    ylabel('molecules')
    ylim([0 60])
subplot(3,3,9)
    plot(runssamp(runnum).t/60,runssamp(runnum).cIII,'k')
    title('cIII')
    xlabel('min'); 
    ylabel('molecules')
    ylim([0 150])
    
% plot(1:nbins,runssamp(runnum).t)

%% Plot raw data
if plot_trace    
    for runnum=39
        close all
%         load([defaultDIR defaultprefix '_' num2str(runnum,'%.3u')], 'run');        
        figure('Name',['Run ' num2str(runnum)])
        hold on
        subplot(3,3,1)
            plot(runs(runnum).t/60,runs(runnum).PRmRNA,'r')
            title('PR mRNA')
            % xlabel('min'); 
            ylabel('molecules')
            ylim([0 10])
        subplot(3,3,4)
            plot(runs(runnum).t/60,runs(runnum).cro_2 ./ (A*runs(runnum).V) * 10^9,'k')
            title('cro_2')
        %     xlabel('min'); 
            ylabel('molecules')
            ylim([0 200])
        subplot(3,3,7)
            plot(runs(runnum).t/60,runs(runnum).cII ./ (A*runs(runnum).V) * 10^9,'k')
            title('cII')
            xlabel('min'); 
            ylabel('molecules')
            ylim([0 150])
        %
%         figure('Name',['PRM' num2str(runnum)])
%         hold on
        subplot(3,3,2)
            plot(runs(runnum).t/60,runs(runnum).PREmRNA,'r')
            title('PRE mRNA')
        %     xlabel('min'); 
            ylabel('molecules')
            ylim([0 10])                    
        subplot(3,3,5)
            plot(runs(runnum).t/60,runs(runnum).PRMmRNA,'r')
            title('PRM mRNA')
            % xlabel('min'); 
            ylabel('molecules')
            ylim([0 10])
        subplot(3,3,8)
            plot(runs(runnum).t/60,runs(runnum).cI_2 ./ (A*runs(runnum).V) * 10^9,'k')
            title('cI_2')
            xlabel('min'); 
            ylabel('molecules')
            ylim([0 100])
        %
%         figure('Name',['PL' num2str(runnum)])
%         hold on
        subplot(3,3,3)
            plot(runs(runnum).t/60,runs(runnum).PLmRNA,'r')
            title('PL mRNA')
            % xlabel('min'); 
            ylabel('molecules')
            ylim([0 10])
        subplot(3,3,6)
            plot(runs(runnum).t/60,runs(runnum).N ./ (A*runs(runnum).V) * 10^9,'k')
            title('N')
        %     xlabel('min'); 
            ylabel('molecules')
            ylim([0 60])
        subplot(3,3,9)
            plot(runs(runnum).t/60,runs(runnum).cIII ./ (A*runs(runnum).V) * 10^9,'k')
            title('cIII')
            xlabel('min'); 
            ylabel('molecules')
            ylim([0 150])
            
            % Save Figs
        FigDir='C:\Users\frink\OneDrive\OneNoteMatlab\Golding\Gillespie\Outputs\20150911_Figs\';
        af_saveallcurrentfig(FigDir)
        af_saveallcurrentfig(FigDir,'','png')        
    end
end



%% Plot by run
if plot_trace    
    for runnum=1
        close all
%         load([defaultDIR defaultprefix '_' num2str(runnum,'%.3u')], 'run');        
        figure('Name',['Run ' num2str(runnum)])
        hold on
        subplot(3,3,1)
            plot(runs(runnum).t/60,runs(runnum).PRmRNA,'r')
            title('PR mRNA')
            % xlabel('min'); 
            ylabel('molecules')
            ylim([0 10])
        subplot(3,3,4)
            plot(runs(runnum).t/60,runs(runnum).cro_2 ./ (A*runs(runnum).V) * 10^9,'k')
            title('cro_2')
        %     xlabel('min'); 
            ylabel('nM')
            ylim([0 200])
        subplot(3,3,7)
            plot(runs(runnum).t/60,runs(runnum).cII ./ (A*runs(runnum).V) * 10^9,'k')
            title('cII')
            xlabel('min'); 
            ylabel('nM')
            ylim([0 150])
        %
%         figure('Name',['PRM' num2str(runnum)])
%         hold on
        subplot(3,3,2)
            plot(runs(runnum).t/60,runs(runnum).PREmRNA,'r')
            title('PRE mRNA')
        %     xlabel('min'); 
            ylabel('molecules')
            ylim([0 10])                    
        subplot(3,3,5)
            plot(runs(runnum).t/60,runs(runnum).PRMmRNA,'r')
            title('PRM mRNA')
            % xlabel('min'); 
            ylabel('molecules')
            ylim([0 10])
        subplot(3,3,8)
            plot(runs(runnum).t/60,runs(runnum).cI_2 ./ (A*runs(runnum).V) * 10^9,'k')
            title('cI_2')
            xlabel('min'); 
            ylabel('nM')
            ylim([0 100])
        %
%         figure('Name',['PL' num2str(runnum)])
%         hold on
        subplot(3,3,3)
            plot(runs(runnum).t/60,runs(runnum).PLmRNA,'r')
            title('PL mRNA')
            % xlabel('min'); 
            ylabel('molecules')
            ylim([0 10])
        subplot(3,3,6)
            plot(runs(runnum).t/60,runs(runnum).N ./ (A*runs(runnum).V) * 10^9,'k')
            title('N')
        %     xlabel('min'); 
            ylabel('nM')
            ylim([0 60])
        subplot(3,3,9)
            plot(runs(runnum).t/60,runs(runnum).cIII ./ (A*runs(runnum).V) * 10^9,'k')
            title('cIII')
            xlabel('min'); 
            ylabel('nM')
            ylim([0 150])
            
            % Save Figs
%         FigDir='C:\Users\frink\OneDrive\OneNoteMatlab\Golding\Gillespie\Outputs\20150911_Figs\';
%         af_saveallcurrentfig(FigDir)
%         af_saveallcurrentfig(FigDir,'','png')        
    end
end


%% Plot by promoter
if plot_trace
    close all
    for runnum=1:10
%         load([defaultDIR defaultprefix '_' num2str(runnum,'%.3u')], 'run');        
        figure('Name',['PR' num2str(runnum)])
        hold on
        subplot(3,1,1)
            plot(runs(runnum).t/60,runs(runnum).PRmRNA,'r')
            title('PR mRNA')
            % xlabel('min'); 
            ylabel('molecules')
        subplot(3,1,2)
            plot(runs(runnum).t/60,runs(runnum).cro_2,'k')
            title('cro_2')
        %     xlabel('min'); 
            ylabel('molecules')
        subplot(3,1,3)
            plot(runs(runnum).t/60,runs(runnum).cII,'k')
            title('cII')
            xlabel('min'); 
            ylabel('molecules')
        %
        figure('Name',['PRM' num2str(runnum)])
        hold on
        subplot(3,1,1)
            plot(runs(runnum).t/60,runs(runnum).PRMmRNA,'r')
            title('PRM mRNA')
            % xlabel('min'); 
            ylabel('molecules')
        subplot(3,1,2)
            plot(runs(runnum).t/60,runs(runnum).PREmRNA,'r')
            title('PRE mRNA')
        %     xlabel('min'); 
            ylabel('molecules')
        subplot(3,1,3)
            plot(runs(runnum).t/60,runs(runnum).cI_2,'k')
            title('cI_2')
            xlabel('min'); 
            ylabel('molecules')
        %
        figure('Name',['PL' num2str(runnum)])
        hold on
        subplot(3,1,1)
            plot(runs(runnum).t/60,runs(runnum).PLmRNA,'r')
            title('PL mRNA')
            % xlabel('min'); 
            ylabel('molecules')
        subplot(3,1,2)
            plot(runs(runnum).t/60,runs(runnum).N,'k')
            title('N')
        %     xlabel('min'); 
            ylabel('molecules')
        subplot(3,1,3)
            plot(runs(runnum).t/60,runs(runnum).cIII,'k')
            title('cIII')
            xlabel('min'); 
            ylabel('molecules')
    end
end


%% Unused

% for runnum=1:100
%     clear run
%     disp(['Run ' num2str(runnum)])
%     load([defaultDIR defaultprefix '_' num2str(runnum,'%.3u')], 'run');
%     temp = run(runnum);
%     clear run
%     run = temp;
%     save([defaultDIR defaultprefix '_' num2str(runnum,'%.3u')], 'run');
% end




