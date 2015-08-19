% Gillespie 1 Gene No Feedback Txn/Tln
close all
clear all

%% Setup Parameters
%   k_r k_p g_r g_p
k = [.5 .6 .05 .1];
rxnm = [1 0 0; 0 1 0; 0 1 0; 0 0 1];
scm = [0 1 0; 0 0 1; 0 -1 0; 0 0 -1];
inspc = [1 0 0];
rpts = 6;
maxstep = 50000;

%% Analytical solutions

calc_mean(1) = inspc(1);
calc_var(1) = 0;

calc_mean(2) = k(1)/k(3);
calc_var(2) = k(1)/k(3);

calc_mean(3) = k(1)*k(2)/(k(3)*k(4));
calc_var(3) = k(1)*k(2)/(k(3)*k(4)) + k(1)*k(2)^2/(k(3)*k(4)^2+k(3)^2*k(4));

%% Run SSA

[tout, spcout ] = Gillespie3(k,rxnm, scm, inspc, rpts, maxstep);


%% Display traces
figure('Name','1 Gene No Feedback')
hold on
title('1 Gene No Feedback')
for iter1 = 1:rpts    
    h_sp1(iter1) = subplot(2,3,iter1);
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
ext_spc_means = repmat(spc_means,[1 maxstep 1]);
spc_var = sum(timenorm.*(spcout-ext_spc_means).^2,2);
m_spc_var = mean(spc_var,3)';

% Calc mean of trace ignoring first % of it
ignore_per = 0.3;
spc_means = sum(spcout(:,maxstep*ignore_per+1:end,:).*timenorm(:,maxstep*ignore_per+1:end,:),2);
m_spc_mean = mean(spc_means,3)';
ext_spc_means = repmat(spc_means,[1 maxstep*(1-ignore_per) 1]);
spc_var = sum(timenorm(:,maxstep*ignore_per+1:end,:).*(spcout(:,maxstep*ignore_per+1:end,:)-ext_spc_means).^2,2);
m_spc_var = mean(spc_var,3)';


m_spc_mean
calc_mean

m_spc_var
calc_var

output(1,:) = numel(inspc);
% output(2,:) = 

%% Dynamic
A=-k(1)/k(3);
B=-k(1)*k(2)/(k(3)*k(4)-k(4)^2);
m_fun_t = @(t) k(1)/k(3)+A*exp(-k(3)*t);
p_fun_t = @(t) k(1)*k(2)/(k(3)*k(4))+k(2)*A*exp(-k(3)*t)/(k(4)-k(3))+B*exp(-k(4)*t);

% Display species
figure('Name','1 Gene No Feedback')
hold on
title('1 Gene No Feedback')
for iter1 = 1:rpts    
    h_sp2(iter1) = subplot(2,3,iter1);
    hold on
    plot(tout(1,:,1),spcout(2,:,iter1),'r')
    plot(tout(1,:,1),spcout(3,:,iter1),'k')
    fplot(m_fun_t,[0 400],'r:'); 
    fplot(p_fun_t,[0 400],'k:');
    xlabel('time (A.U.)')
    ylabel('# molecules')
end

linkaxes(h_sp2)
ylim([-5 max(spcout(:))*1.01])
legend('mRNA','Protein')
annotation('textbox',[.01 0 .1 .1],'String',{'k = ', num2str(k(1:2)), num2str(k(3:4))})


