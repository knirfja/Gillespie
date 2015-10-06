% Run LambdaModel OLD
%% Initialize
% clear all %#ok<*CLSCR>
close all
rannum=floor(rand(1)*10000000);
A=6.02*10^23;
% AV=6.02*10^23*1.5*10^-15;
%%
numruns=200;

plot_meanstd=1;
plot_rawtrace=0;
plot_cIvscro=1;
% runnum=1;
%
% numruns=numruns+1;
defaultprefix = datestr(date,['yyyy','mm','dd']);
% defaultDIR = [pwd '\'];
% defaultDIR = ['C:\Users\frink\OneDrive\OneNoteMatlab\Golding\Gillespie\Outputs\' defaultprefix '\'];
defaultDIR = ['C:\Users\frink\Documents\MATLAB\Lambda\LambdaModel\Outputs\' defaultprefix '\'];
%%
% for runnum=138:numruns
% while runnum<=numruns
    tic
    clear run
    disp(['Run ' num2str(runnum)])
    rannum=floor(rand(1)*10000000);
    run = LambdaModel(rannum);
    % Save
    if ~exist(defaultDIR,'dir')
        mkdir(defaultDIR)
    end
    save([defaultDIR defaultprefix '_' num2str(runnum,'%.3u')], 'run');
%     runnum=runnum+1;
    runtime(runnum)=toc/60

%% sample each variable
for runnum=1:numruns
    clear run
    load([defaultDIR defaultprefix '_' num2str(runnum,'%.3u') '.mat'])

    maxt=round(run.t(end));
    nsamp=maxt; % 1 sec bins
    % t=run.t./60;
    % binrange = linspace(min(t),max(t),nbins);
    % [n,bin] = histc(t,binrange);
    % n = histc(t,binrange);
    for iter1 = 1:nsamp
        ind=find(run.t<iter1,1,'last');
        runsamp(runnum).t(iter1)=mean(run.t(ind));
        runsamp(runnum).V(iter1)=mean(run.V(ind));
        runsamp(runnum).cI(iter1)=mean(run.cI(ind));
        runsamp(runnum).cro(iter1)=mean(run.cro(ind));    
        runsamp(runnum).cI_2(iter1)=mean(run.cI_2(ind));
        runsamp(runnum).cro_2(iter1)=mean(run.cro_2(ind));
        runsamp(runnum).N(iter1)=mean(run.N(ind));
        runsamp(runnum).cII(iter1)=mean(run.cII(ind));
        runsamp(runnum).cIII(iter1)=mean(run.cIII(ind));
        runsamp(runnum).PRmRNA(iter1)=mean(run.PRmRNA(ind));
        runsamp(runnum).PREmRNA(iter1)=mean(run.PREmRNA(ind));
        runsamp(runnum).PRMmRNA(iter1)=mean(run.PRMmRNA(ind));
        runsamp(runnum).PLmRNA(iter1)=mean(run.PLmRNA(ind));   
        
    end
    runsamp(runnum).CroPerPRmRNA=run.CroPerPRmRNA;
    runsamp(runnum).NPerPLmRNA=run.NPerPLmRNA;
    runsamp(runnum).PRmRNALife=run.PRmRNALife;
    runsamp(runnum).PREmRNALife=run.PREmRNALife;
    runsamp(runnum).PRMmRNALife=run.PRMmRNALife;
    runsamp(runnum).PLmRNALife=run.PLmRNALife;
    
%     runsamp(runnum).PRlifetime=run.PRmRNALife(:,2)-run.PRmRNALife(:,1);
%     runsamp(runnum).PRElifetime=run.PREmRNALife(:,2)-run.PREmRNALife(:,1);
%     runsamp(runnum).PRMlifetime=run.PRMmRNALife(:,2)-run.PRMmRNALife(:,1);
%     runsamp(runnum).PLlifetime=run.PLmRNALife(:,2)-run.PLmRNALife(:,1);
    
    runsamp(runnum).randseed=run.randseed;
    
    
    save([defaultDIR defaultprefix '_' num2str(runnum,'%.3u')], 'runsamp','-append');
%     %%
%     figure; hist(runsamp(runnum).PRlifetime)
%     figure; hist(runsamp(runnum).PRElifetime)
%     figure; hist(runsamp(runnum).PRMlifetime)
%     figure; hist(runsamp(runnum).PLlifetime)

    %% Plot and save figs
    close all
    maxt=round(runsamp(runnum).t(end));
    figure('Name',['Run_' num2str(runnum) '_Sampled'])
    ssp(1)=subplot(3,3,1);
        plot(runsamp(runnum).t/60,runsamp(runnum).PRmRNA,'r')
        title('PR mRNA')
        % xlabel('min'); 
        ylabel('molecules')
        ylim([0 10])
        xlim([0 maxt/60])
    ssp(4)=subplot(3,3,4);
        plot(runsamp(runnum).t/60,runsamp(runnum).cro_2 ./ (A*runsamp(runnum).V) * 10^9,'k')
        title('cro_2')
    %     xlabel('min'); 
        ylabel('nM')
        ylim([0 200])
        xlim([0 maxt/60])
    ssp(7)=subplot(3,3,7);
        plot(runsamp(runnum).t/60,runsamp(runnum).cII ./ (A*runsamp(runnum).V) * 10^9,'k')
        title('cII')
        xlabel('min'); 
        ylabel('nM')
        ylim([0 150])
        xlim([0 maxt/60])

    ssp(2)=subplot(3,3,2);
        plot(runsamp(runnum).t/60,runsamp(runnum).PREmRNA,'r')
        title('PRE mRNA')
        % xlabel('min'); 
        ylabel('molecules')
        ylim([0 10])
        xlim([0 maxt/60])
    ssp(5)=subplot(3,3,5);
        plot(runsamp(runnum).t/60,runsamp(runnum).PRMmRNA,'r')
        title('PRM mRNA')
        % xlabel('min'); 
        ylabel('molecules')
        ylim([0 10])
        xlim([0 maxt/60])
    ssp(8)=subplot(3,3,8);
        plot(runsamp(runnum).t/60,runsamp(runnum).cI_2 ./ (A*runsamp(runnum).V) * 10^9,'k')
        title('cI_2')
        xlabel('min'); 
        ylabel('nM')
        ylim([0 100])
        xlim([0 maxt/60])

    ssp(3)=subplot(3,3,3);
        plot(runsamp(runnum).t/60,runsamp(runnum).PLmRNA,'r')
        title('PL mRNA')
        % xlabel('min'); 
        ylabel('molecules')
        ylim([0 10])
        xlim([0 maxt/60])
    ssp(6)=subplot(3,3,6);
        plot(runsamp(runnum).t/60,runsamp(runnum).N ./ (A*runsamp(runnum).V) * 10^9,'k')
        title('N')
    %     xlabel('min'); 
        ylabel('nM')
        ylim([0 80])
        xlim([0 maxt/60])
    ssp(1)=subplot(3,3,9);
        plot(runsamp(runnum).t/60,runsamp(runnum).cIII ./ (A*runsamp(runnum).V) * 10^9,'k')
        title('cIII')
        xlabel('min'); 
        ylabel('nM')
        ylim([0 150])
        xlim([0 maxt/60])
                
        FigDir=[defaultDIR 'Figs\'];
        if ~exist(FigDir,'dir')
            mkdir(FigDir)
        end
        af_saveallcurrentfig(FigDir,defaultprefix)
        af_saveallcurrentfig(FigDir,defaultprefix,'png')    
        
end
close all

% load([defaultDIR defaultprefix '_' num2str(runnum,'%.3u') '.mat'])

%% Analyze mean / std of each
t=mean(vertcat(runsamp(:).t))/60;
V=mean(vertcat(runsamp(:).V));
all.cI=vertcat(runsamp(:).cI);
all.cro=vertcat(runsamp(:).cro);
all.cI_2=vertcat(runsamp(:).cI_2);
all.cro_2=vertcat(runsamp(:).cro_2);
all.N=vertcat(runsamp(:).N);
all.cII=vertcat(runsamp(:).cII);
all.cIII=vertcat(runsamp(:).cIII);

if plot_meanstd
%     close all
    figure('Name','meancro2cI2'); 
    hold on
    h_cro=plot(t,(10^9*mean(all.cro_2)./(A*V)),'k');
    plot(t,(10^9*(mean(all.cro_2)+std(all.cro_2))./(A*V)),'Color',[.9 .9 .9],'MarkerSize',0.1);
    plot(t,(10^9*(mean(all.cro_2)-std(all.cro_2))./(A*V)),'Color',[.9 .9 .9],'MarkerSize',0.1);
    
    h_cI=plot(t,(10^9*mean(all.cI_2)./(A*V)),'r');
    plot(t,(10^9*(mean(all.cI_2)+std(all.cI_2))./(A*V)),'Color',[.9 .8 .8],'MarkerSize',0.1);
    plot(t,(10^9*(mean(all.cI_2)-std(all.cI_2))./(A*V)),'Color',[.9 .8 .8],'MarkerSize',0.1);
    
    ylim([0 100])
    
    xlabel('Time (min)'); ylabel('nanomolar')    
    legend([h_cro,h_cI],{'cro2','cI2'},'Location','best')
   
    figure('Name','meancIIcIII'); 
    hold on
    h_cII=plot(t,(10^9*mean(all.cII)./(A*V)),'Color',[0 0 1]);
    plot(t,(10^9*(mean(all.cII)+std(all.cII))./(A*V)),'Color',[.8 .8 .9],'MarkerSize',0.1);
    plot(t,(10^9*(mean(all.cII)-std(all.cII))./(A*V)),'Color',[.8 .8 .9],'MarkerSize',0.1);

    h_cIII=plot(t,(10^9*mean(all.cIII)./(A*V)),'Color',[1 0 1]);
    plot(t,(10^9*(mean(all.cIII)+std(all.cIII))./(A*V)),'Color',[.9 .8 .9],'MarkerSize',0.1);
    plot(t,(10^9*(mean(all.cIII)-std(all.cIII))./(A*V)),'Color',[.9 .8 .9],'MarkerSize',0.1);
    
    h_N=plot(t,(10^9*mean(all.N)./(A*V)),'Color',[1 1 0]);
    plot(t,(10^9*(mean(all.N)+std(all.N))./(A*V)),'Color',[.9 .9 .8],'MarkerSize',0.1);
    plot(t,(10^9*(mean(all.N)-std(all.N))./(A*V)),'Color',[.9 .9 .8],'MarkerSize',0.1);
    
    ylim([0 70])
    
    xlabel('Time (min)'); ylabel('nanomolar')
    legend([h_cII,h_cIII,h_N],{'free cII','free cIII','N'},'Location','best')
    
%     FigDir=[defaultDIR 'Figs\'];
%     if ~exist(FigDir,'dir')
%         mkdir(FigDir)
%     end
%     af_saveallcurrentfig(FigDir,defaultprefix)
%     af_saveallcurrentfig(FigDir,defaultprefix,'png') 
end

%% Percent Lysogen
if plot_cIvscro
    figure('Name','cI_2_vs_cro_2'); 
    plot(all.cI_2(:,end),all.cro_2(:,end),'.')
    xlim([0 110])
    ylim([0 110])
    xlabel('cI_2 (nM)')
    ylabel('cro_2 (nM)')

    % cI2 vs cro2 Histograms
    cro2cI2_binctrs={[(0:10:105)+5]',[(0:10:105)+5]'};
    figure('Name','heatmap_cI_2_vs_cro_2');
    hist3([all.cI_2(:,end) all.cro_2(:,end)],'Ctrs',cro2cI2_binctrs,'FaceAlpha',.65);
    xlabel('cI_2 (nM)'), ylabel('cro_2 (nM)')
    xlim([0 110])
    ylim([0 110])
    set(gcf,'renderer','opengl');
    set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
end



%% Plot raw curves  
if plot_rawtrace    
    for runnum=1
%         close all
%         load([defaultDIR defaultprefix '_' num2str(runnum,'%.3u')], 'run');        
        figure('Name',['Run ' num2str(runnum)])
        hold on
        subplot(3,3,1)
            plot(run.t/60,run.PRmRNA,'r')
            title('PR mRNA')
            % xlabel('min'); 
            ylabel('molecules')
            ylim([0 10])
        subplot(3,3,4)
            plot(run.t/60,run.cro_2 ./ (A*run.V) * 10^9,'k')
            title('cro_2')
        %     xlabel('min'); 
            ylabel('molecules')
            ylim([0 200])
        subplot(3,3,7)
            plot(run.t/60,run.cII ./ (A*run.V) * 10^9,'k')
            title('cII')
            xlabel('min'); 
            ylabel('molecules')
            ylim([0 150])
        %
%         figure('Name',['PRM' num2str(runnum)])
%         hold on
        subplot(3,3,2)
            plot(run.t/60,run.PREmRNA,'r')
            title('PRE mRNA')
        %     xlabel('min'); 
            ylabel('molecules')
            ylim([0 10])                    
        subplot(3,3,5)
            plot(run.t/60,run.PRMmRNA,'r')
            title('PRM mRNA')
            % xlabel('min'); 
            ylabel('molecules')
            ylim([0 10])
        subplot(3,3,8)
            plot(run.t/60,run.cI_2 ./ (A*run.V) * 10^9,'k')
            title('cI_2')
            xlabel('min'); 
            ylabel('molecules')
            ylim([0 100])
        %
%         figure('Name',['PL' num2str(runnum)])
%         hold on
        subplot(3,3,3)
            plot(run.t/60,run.PLmRNA,'r')
            title('PL mRNA')
            % xlabel('min'); 
            ylabel('molecules')
            ylim([0 10])
        subplot(3,3,6)
            plot(run.t/60,run.N ./ (A*run.V) * 10^9,'k')
            title('N')
        %     xlabel('min'); 
            ylabel('molecules')
            ylim([0 60])
        subplot(3,3,9)
            plot(run.t/60,run.cIII ./ (A*run.V) * 10^9,'k')
            title('cIII')
            xlabel('min'); 
            ylabel('molecules')
            ylim([0 150])
            
            % Save Figs
        FigDir='C:\Users\frink\OneDrive\OneNoteMatlab\Golding\Gillespie\Outputs\20150911_Figs\';
%         af_saveallcurrentfig(FigDir)
%         af_saveallcurrentfig(FigDir,'','png')        
    end
end
    

%% Sample from saved data
% for runnum=1:numruns
%     clear run
%     load([defaultDIR defaultprefix '_' num2str(runnum,'%.3u') '.mat'])
% 
%     maxt=round(run.t(end));
%     nsamp=maxt; % 1 sec bins
%     % t=run.t./60;
%     % binrange = linspace(min(t),max(t),nbins);
%     % [n,bin] = histc(t,binrange);
%     % n = histc(t,binrange);
%     for iter1 = 1:nsamp
%         ind=find(run.t<iter1,1,'last');
%         runsamp(runnum).t(iter1)=mean(run.t(ind));
%         runsamp(runnum).V(iter1)=mean(run.V(ind));
%         runsamp(runnum).cI(iter1)=mean(run.cI(ind));
%         runsamp(runnum).cro(iter1)=mean(run.cro(ind));    
%         runsamp(runnum).cI_2(iter1)=mean(run.cI_2(ind));
%         runsamp(runnum).cro_2(iter1)=mean(run.cro_2(ind));
%         runsamp(runnum).N(iter1)=mean(run.N(ind));
%         runsamp(runnum).cII(iter1)=mean(run.cII(ind));
%         runsamp(runnum).cIII(iter1)=mean(run.cIII(ind));
%         runsamp(runnum).PRmRNA(iter1)=mean(run.PRmRNA(ind));
%         runsamp(runnum).PREmRNA(iter1)=mean(run.PREmRNA(ind));
%         runsamp(runnum).PRMmRNA(iter1)=mean(run.PRMmRNA(ind));
%         runsamp(runnum).PLmRNA(iter1)=mean(run.PLmRNA(ind));
%         runsamp(runnum).randseed=run.randseed;
%     end
% end
%     
% 
%% Time trace cI2 vs cro2
%% single traces
for runnum=1:10;
    figure('Name',['CIcroTrace_Run_' num2str(runnum)])
    x=runsamp(runnum).cI_2;
    y=runsamp(runnum).cro_2;
    z=zeros(size(x));
    col=runsamp(runnum).t;        

    % x = 0:.05:2*pi;
    % y = sin(x);
    % z = zeros(size(x));
    % col = x;  % This is the color, vary with x in this case.
    surface([x;x],[y;y],[z;z],[col;col],...
            'facecol','no',...
            'edgecol','interp',...
            'Marker','.','MarkerSize',10,...
            'LineStyle',':','LineWidth',0.05);
    xlim([0 120])
    ylim([0 120])
    xlabel('cI_2 (nM)')
    ylabel('cro_2 (nM)')
     
end

%% Pop mean trace
t=mean(vertcat(runsamp(:).t))/60;
V=mean(vertcat(runsamp(:).V));
all.cI=vertcat(runsamp(:).cI);
all.cro=vertcat(runsamp(:).cro);
all.cI_2=vertcat(runsamp(:).cI_2);
all.cro_2=vertcat(runsamp(:).cro_2);
all.N=vertcat(runsamp(:).N);
all.cII=vertcat(runsamp(:).cII);
all.cIII=vertcat(runsamp(:).cIII);

figure('Name','MeanCIcroTrace')
x=(10^9*mean(all.cI_2)./(A*V));
y=(10^9*mean(all.cro_2)./(A*V));
z=zeros(size(x));
col=t;  
    
%      h_cro=plot(t,(10^9*mean(all.cro_2)./(A*V)),'k');
%     plot(t,(10^9*(mean(all.cro_2)+std(all.cro_2))./(A*V)),'Color',[.9 .9 .9],'MarkerSize',0.1);
%     plot(t,(10^9*(mean(all.cro_2)-std(all.cro_2))./(A*V)),'Color',[.9 .9 .9],'MarkerSize',0.1);
    
%     h_cI=plot(t,(10^9*mean(all.cI_2)./(A*V)),'r');

surface([x;x],[y;y],[z;z],[col;col],...
            'facecol','no',...
            'edgecol','interp',...
            'Marker','.','MarkerSize',10,...
            'LineStyle',':','LineWidth',0.05);
xlim([0 120])
ylim([0 120])
xlabel('cI_2 (nM)')
ylabel('cro_2 (nM)')

%% Protein per mRNA

CroPerPRmRNA=[runsamp(:).CroPerPRmRNA];
NPerPLmRNA=[runsamp(:).NPerPLmRNA];

figure('Name','CroPerM'); hold on
%     hist(CroPerPRmRNA,20)
h= histogram(CroPerPRmRNA,'Normalization','probability');
xlabel('Cro produced per mRNA'); ylabel('fraction of population')
str1=['Mean= ' num2str(mean(CroPerPRmRNA)) ', Std= ' num2str(std(CroPerPRmRNA))];    
% annotation('textbox',[0.5 0.5 0.2 0.2],'String',str1,'FitBoxToText','on');


%     num=0:35;
%     f=@(x) (0.2*exp(-0.2)*3^x )/factorial(x);
%     plot(num,arrayfun(f,num),'r')
%     fNt=@(t,k) 0.2*0.6^k/factorial(k)*t^k*exp(-(0.2+0.6)*t);
fN=@(k) 0.2*0.9^k/((0.2+0.9)^(k+1));
fplot(fN,[0 35],'r')        
str2=['Analytical Mean= ' num2str(0.9/0.2) ' ± 1.5'];
% annotation('textbox',[0.5 0.4 0.2 0.2],'String',str2,'FitBoxToText','on');
legend(str1,str2);

figure('Name','NPerM'); hold on
%     hist(CroPerPRmRNA,20)
h= histogram(NPerPLmRNA,'Normalization','probability');
xlabel('N produced per mRNA'); ylabel('fraction of population')    
str1=['Mean= ' num2str(mean(NPerPLmRNA)) ', Std= ' num2str(std(NPerPLmRNA))];    
% annotation('textbox',[0.5 0.5 0.2 0.2],'String',str1,'FitBoxToText','on');


 fN=@(k) 0.2*0.9^k/((0.2+0.9)^(k+1));
fplot(fN,[0 35],'r')        
str2=['Analytical Mean= ' num2str(0.9/0.2) ' ± 1.5'];
% annotation('textbox',[0.5 0.4 0.2 0.2],'String',str2,'FitBoxToText','on');
legend(str1,str2);

%% mRNA lifetime
%%% Genetic locations 
% Rightward genes
PR_trans=[38023 38653+20];
% PR_trans=[38023 40624]; % -> Actual transcript length
% tr0 terminator 38135
cro_span=[38041,38241]; % ->
NutR=38265;
% NutR 38265-38281
TR1=38337;
% tr1a terminator 38315
% tr1b terminator 38337
% tr1c terminator 38370
cII_span=[38360,38653]; % ->

% Leftward genes
PRE_trans=[38343 37227-20]; % <-
% PRE_trans=[38343 35798]; % <- Actual transcript length
% anticro_span=[38241,38041]; % <-
PRM_trans=[37958 37227-20];
% PRM_trans=[37958 35798]; % <- Actual transcript length
cI_span=[37940,37227]; % <-

PL_trans=[35600 33299-20];
% PL_trans=[35600 33100]; % <- Actual transcript length
NutL=35518;
% NutL 35534-35518 % <-
N_span=[35582,34560]; % <-
TL1=34560;
% tl1 terminator 34560
% tl2a terminator 33930
% tl2b terminator 33494
cIII_span=[33463,33299]; % <-
% tl2c terminator 33141
% tl2d terminator 33100
% tl3 terminator 31262
% ti terminator 27538
% j1 terminator 18671
RNAPrate=30; %nt/s
% meanshift=
38041-PR_trans(1);

for iter=1:numel(runsamp)
    tempPR{iter} = runsamp(iter).PRmRNALife(:,2)-runsamp(iter).PRmRNALife(:,1);
    tempPR{iter}=tempPR{iter}';       
    tempPRE{iter} = runsamp(iter).PREmRNALife(:,2)-runsamp(iter).PREmRNALife(:,1);
    tempPRE{iter}=tempPRE{iter}';     
    tempPRM{iter} = runsamp(iter).PRMmRNALife(:,2)-runsamp(iter).PRMmRNALife(:,1);
    tempPRM{iter}=tempPRM{iter}';     
    tempPL{iter} = runsamp(iter).PLmRNALife(:,2)-runsamp(iter).PLmRNALife(:,1);
    tempPL{iter}=tempPL{iter}';     
end

PRLifetime=[tempPR{:}];
PRELifetime=[tempPRE{:}];
PRMLifetime=[tempPRM{:}];
PLLifetime=[tempPL{:}];
clear tempPR tempPRE tempPRM tempPL

fp=@(lam,c,x) lam*exp(-lam*(x-c));

%
g1=fittype(@(lam,x) fp(lam,18/30,x)); %
% figure;
sax(1)=subplot(2,2,1);
h(1)= histogram(PRLifetime,'Normalization','probability','BinMethod','integers');
% h= histogram(PRLifetime,'BinMethod','integers');
hold on
hx=(1:h(1).NumBins)';
% f1=fit(hx(2:end),h.Values(2:end)',g1,'StartPoint',0.3);
f1=fit(hx(2:end),h(1).Values(2:end)',g1,'StartPoint',[0.3]);
err=(f1(0)+f1(1)-h(1).Values(1))/(numel(hx)-1);
% plot(2:hx(end),f1(2:hx(end))+err,'r')
plot(hx(2:end),f1(hx(2:end)),'k')
halflife= log(2)/coeffvalues(f1)+18/30;
halflifeCI= log(2)./confint(f1)+18/30;
plot([halflife halflife],[0 max(h(1).Values)],'r--')
legend('Sim',...
    ['Exp fit,\lambda = ' num2str(coeffvalues(f1),3) ' ± ' num2str(sum(abs(confint(f1)-coeffvalues(f1)))/2,2)],...
       ['Fit Half-life = ' num2str(halflife,3) ' ± ' num2str(sum(abs(halflifeCI-halflife))/2,2) ' s'])
%     fit([cumsum(diff(h.BinEdges))./2]',h.Values','exp1')
title('PR mRNA')    
%%
sax(2)=subplot(2,2,2);
g1=fittype(@(lam,x) fp(lam,403/30,x)); %
h(2)= histogram(PRELifetime,30,'Normalization','probability','BinMethod','integers');
hold on
hx=(1:h(2).NumBins)';
% f1=fit(hx(2:end),h.Values(2:end)',g1,'StartPoint',0.3);
f1=fit(hx(2:end),h(2).Values(2:end)',g1,'StartPoint',[0.3]);
err=(f1(0)+f1(1)-h(2).Values(1))/(numel(hx)-1);
% plot(2:hx(end),f1(2:hx(end))+err,'r')
% plot(hx(2:end)+403/30,f1(hx(2:end)),'k')
halflife= log(2)/coeffvalues(f1)+403/30;
halflifeCI= log(2)./confint(f1)+403/30;
% plot([halflife halflife],[0 max(h(2).Values)],'r--')
legend('Sim',...
    ['Exp fit,\lambda = ' num2str(coeffvalues(f1),3) ' ± ' num2str(sum(abs(confint(f1)-coeffvalues(f1)))/2,2)],...
       ['Fit Half-life = ' num2str(halflife,3) ' ± ' num2str(sum(abs(halflifeCI-halflife))/2,2) ' s'])
%     fit([cumsum(diff(h.BinEdges))./2]',h.Values','exp1')
title('PRE mRNA')
%%

sax(3)=subplot(2,2,3);
h= histogram(PRMLifetime,30,'Normalization','probability');
title('PRM mRNA')

sax(4)=subplot(2,2,4);
h= histogram(PLLifetime,30,'Normalization','probability');
title('PL mRNA')
linkaxes(sax)
    
%% Save Figs
FigDir=[defaultDIR 'Figs\'];
if ~exist(FigDir,'dir')
    mkdir(FigDir)
end
af_saveallcurrentfig(FigDir,defaultprefix)
af_saveallcurrentfig(FigDir,defaultprefix,'png') 
   
%% Unused
% %% Rename files
% for i=138:200
%     movefile([defaultprefix '_Run_' num2str(i,3) '_Sampled_1' '.png'],[OLDdefaultprefix '_Run_' num2str(i,3) '_Sampled_1' '.png'])
% end
    
    
    
    