% Gillespie 1 Gene No Feedback Txn/Tln
close all
clear all

%% Setup Parameters
%   k_r k_p g_r g_p
% k = [.5 .6 .05 .1];
k(3) = 0.1824/60; % Nucleic Acids Res. 2013 Jul; 41(13): 6381–6390.
% k(3) = 0.0009/60;
k(4) = 0.0087/60; % Nucleic Acids Res. 2013 Jul; 41(13): 6381–6390.
k(1) = 10 * k(3); % Arbitrary
k(2) = 5 * k(4); % Arbitrary
rxnm = [1 0 0; 0 1 0; 0 1 0; 0 0 1];
scm = [0 1 0; 0 0 1; 0 -1 0; 0 0 -1];
inspc = [1 0 0];
rpts = 10;
maxstep = 10000;

%% Analytical solutions

calc_mean(1) = inspc(1);
calc_var(1) = 0;

calc_mean(2) = k(1)/k(3);
calc_var(2) = k(1)/k(3);

calc_mean(3) = k(1)*k(2)/(k(3)*k(4));
calc_var(3) = k(1)*k(2)/(k(3)*k(4)) + k(1)*k(2)^2/(k(3)*k(4)^2+k(3)^2*k(4));

%% Run SSA

[tout, spcout ] = Gillespie(k,rxnm, scm, inspc, rpts, maxstep);


%% Display traces
figure('Name','1 Gene No Feedback')
hold on
title('1 Gene No Feedback')
for iter1 = 1:rpts    
    h_sp1(iter1) = subplot(rpts/5,rpts/2,iter1);
    hold on
    h_trace1(iter1) = plot(tout(1,:,1),spcout(2,:,iter1),'r');
    h_mean1(iter1) = plot(tout(1,:,1),calc_mean(2).*ones(1,maxstep),'r--');
    h_posstd1(iter1) = plot(tout(1,:,1),(calc_mean(2)+2*sqrt(calc_var(2))) .*ones(1,maxstep),'r:');
    h_negstd1(iter1) = plot(tout(1,:,1),(calc_mean(2)-2*sqrt(calc_var(2))) .*ones(1,maxstep),'r:');
    
    h_trace2(iter1) = plot(tout(1,:,1),spcout(3,:,iter1),'k');
    h_mean2(iter1) = plot(tout(1,:,1),calc_mean(3).*ones(1,maxstep),'k--');
    h_posstd2(iter1) = plot(tout(1,:,1),(calc_mean(3)+2*sqrt(calc_var(3))) .*ones(1,maxstep),'k:');
    h_negstd2(iter1) = plot(tout(1,:,1),(calc_mean(3)-2*sqrt(calc_var(3))) .*ones(1,maxstep),'k:')   ; 
    
    xlabel('time (A.U.)')
    ylabel('# molecules')
end
linkaxes(h_sp1)
ylim([-5 max(spcout(:))*1.01])
xlim([0 max(tout(:))])
legend([h_trace1(end),h_trace2(end)],'mRNA','Protein')
annotation('textbox',[.01 0 .1 .1],'String',{'k = ', num2str(k(1:2)), num2str(k(3:4))})



%% Calculations

timedur=[diff(tout(1,:,:),1) zeros(1,1,rpts)];
timedur=repmat(timedur,[numel(inspc) 1]);
for iter1 = 1:rpts
    timenorm(:,:,iter1) = timedur(:,:,iter1) ./ tout(1,end,iter1);
end

% Calc mean of entire trace
spc_means = sum(spcout.*timenorm,2);
m_spc_mean = mean(spc_means,3)';
std_spc_mean = std(spc_means,0,3)';
ext_spc_means = repmat(spc_means,[1 maxstep 1]);
spc_var = sum(timenorm.*(spcout-ext_spc_means).^2,2);
m_spc_var = mean(spc_var,3)';
std_spc_var = std(spc_var,0,3)';

% Calc mean of trace ignoring first % of it
% ignore_per = 0.3;
% spc_means = sum(spcout(:,maxstep*ignore_per+1:end,:).*timenorm(:,maxstep*ignore_per+1:end,:),2);
% m_spc_mean = mean(spc_means,3)';
% std_spc_mean = std(spc_means,0,3)';
% ext_spc_means = repmat(spc_means,[1 maxstep*(1-ignore_per) 1]);
% spc_var = sum(timenorm(:,maxstep*ignore_per+1:end,:).*(spcout(:,maxstep*ignore_per+1:end,:)-ext_spc_means).^2,2);
% m_spc_var = mean(spc_var,3)';
% std_spc_var = std(spc_var,0,3)';


m_spc_mean
calc_mean

m_spc_var
calc_var

output(1,:) = numel(inspc);
% output(2,:) = 
figure('Name','Steady State Means')
hold on
plot([1:2],calc_mean(2:3),'r.','markersize',15)
errorbar(m_spc_mean(2:3),std_spc_mean(2:3),'k+')
% errorbar([1:3]+.1,m_spc_var,std_spc_var,'v')
ylabel('mean # molecules')
spc_names = {'DNA', 'mRNA','Protein'};
set(gca,'xtick',[1 2])
set(gca,'xticklabel',spc_names(2:3))
legend('Predicted Mean','Measured Mean')
%%
figure('Name','Steady State Variance')
hold on
plot([1:2],calc_var(2:3),'r.','markersize',15)
errorbar(m_spc_var(2:3),std_spc_var(2:3),'k+')
% errorbar([1:3]+.1,m_spc_var,std_spc_var,'v')
ylabel('variance in # molecules')
spc_names = {'DNA', 'mRNA','Protein'};
set(gca,'xtick',[1 2])
set(gca,'xticklabel',spc_names(2:3))
legend('Predicted \sigma ^2','Measured \sigma ^2')


%% Dynamic
A=-k(1)/k(3);
B=-k(1)*k(2)/(k(3)*k(4)-k(4)^2);
m_fun_t = @(t) k(1)/k(3)+A*exp(-k(3)*t);
p_fun_t = @(t) k(1)*k(2)/(k(3)*k(4))+k(2)*A*exp(-k(3)*t)/(k(4)-k(3))+B*exp(-k(4)*t);

% Display species
figure('Name','1 Gene No Feedback Dynamic')
hold on
title('1 Gene No Feedback')
for iter1 = 1:rpts    
    h_sp2(iter1) = subplot(rpts/5,rpts/2,iter1);
    hold on
    plot(tout(1,:,1),spcout(2,:,iter1),'r')
    plot(tout(1,:,1),spcout(3,:,iter1),'k')
    fplot(m_fun_t,[0 10000],'r:'); 
    fplot(p_fun_t,[0 10000],'k:');
    xlabel('time (A.U.)')
    ylabel('# molecules')
end

linkaxes(h_sp2)
ylim([-5 max(spcout(:))*1.01])
legend('mRNA','Protein')
annotation('textbox',[.01 0 .1 .1],'String',{'k = ', num2str(k(1:2)), num2str(k(3:4))})

%%
% Display species
figure('Name','1 Gene No Feedback Dynamic Linear Fit')
hold on
title('1 Gene No Feedback')

% m_fun_t_ln = @(t) log( k(1)/k(3)+A*exp(-k(3)*t));
% p_fun_t_ln = @(t) log( k(1)*k(2)/(k(3)*k(4))+k(2)*A*exp(-k(3)*t)/(k(4)-k(3))+B*exp(-k(4)*t));

eq_ind(1)=1;
for iter1 = 1:rpts    
%     h_sp2(iter1) = subplot(2,3,iter1);    
    eq_ind(iter1+1) = find(spcout(2,:,iter1)>=10,1)-1;    
    m_dyn_lin(1,eq_ind(iter1)+1:eq_ind(iter1+1)+eq_ind(iter1)) = tout(1,1:eq_ind(iter1+1),iter1);
    m_dyn_lin(2,eq_ind(iter1)+1:eq_ind(iter1+1)+eq_ind(iter1)) = -log( abs((calc_mean(2)-spcout(2,1:eq_ind(iter1+1),iter1))/calc_mean(2)) );    
%     ylabel('# molecules')
end
% m_dyn_lin_pooled = 
fit_op = fitoptions('Method','NonlinearLeastSquares',...
           'Lower',[],...
           'Upper',[],...
           'Startpoint',[0]);

%     f = fittype('a*2^(b*x)+c','options',s); 
f = fittype('a*x','options',fit_op);    

[c,gof] = fit(m_dyn_lin(1,:)',m_dyn_lin(2,:)',f);


p = [coeffvalues(c) 0];
    ci=confint(c);
p_low = [ci(1) 0];
p_high = [ci(2) 0];


%     [p,S] = polyfit(m_dyn_lin(1,:),m_dyn_lin(2,:),1);
yfit = polyval(p,m_dyn_lin(1,:));
yfit_low = polyval(p_low,m_dyn_lin(1,:));
yfit_high = polyval(p_high,m_dyn_lin(1,:));
yresid = m_dyn_lin(2,:) - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(m_dyn_lin(2,:))-1) * var(m_dyn_lin(2,:));
rsq = 1 - SSresid/SStotal; 
%     ste = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df);
%     p_h = p; p_h(1) = p(1)+ste(1);

hold on
plot(m_dyn_lin(1,:),m_dyn_lin(1,:).*k(3),'k-')
plot(m_dyn_lin(1,:),m_dyn_lin(2,:),'r.')
plot(m_dyn_lin(1,:),yfit,'r.-')
plot(m_dyn_lin(1,:),yfit_low,'r:')
plot(m_dyn_lin(1,:),yfit_high,'r:')


lgd_str = {['Predicted k_3 = ' num2str(k(3))] ,'Simulated',['Best Fit k_3 = ' num2str(p(1),3) ' ± ' num2str(ci(2)-p(1),3) ]};

legend(lgd_str,'location','NW')

xlabel('time (A.U.)')
% linkaxes(h_sp2)
% ylim([-5 max(spcout(:))*1.01])
% ylim([-7 1])
% xlim([0 800])

annotation('textbox',[.01 0 .1 .1],'String',{'k = ', num2str(k(1:2)), num2str(k(3:4))})


