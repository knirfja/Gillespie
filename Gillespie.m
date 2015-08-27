function [tout, spcout, varargout ] = Gillespie( rates, rxn_matrix, state_change_matrix, varargin )
%Gillespie1 Summary of this function goes here
%   Detailed explanation goes here
% AFrink 20150809

% [t n] = Gillespie2([1 0.1],[0 ; 1],[1 ; -1]);

if nargin < 4
    in_spc = zeros(1,size(rxn_matrix,2));
else
    in_spc = varargin{1};
end

if nargin < 5
    numrpts = 1;
else
    numrpts = varargin{2};
end

if nargin < 6.
    maxstep = 1000;
else
    maxstep = varargin{3};
end

if nargin < 7.
    stepsize = 1;
else
    stepsize = varargin{4};
end

if nargin < 8.
    rseed = 159265357;
else
    rseed = varargin{5};
end
rng(rseed)

numrxns = numel(rates);
numspc = numel(in_spc);
spcout = zeros(numspc, maxstep, numrpts);
tout = zeros(1,maxstep, numrpts);
rnums = rand(2,maxstep,numrpts);

% figure
% hold on

tic
for runnum = 1:numrpts          
    t = 0;
    curstep = 0;    
    spc = in_spc;
    curprob = zeros(1,numrxns);
    while curstep < maxstep
        curstep = curstep+1;
        
        for iter1 = 1:numrxns        
            curprob(iter1) = rates(iter1)*prod(spc .^ rxn_matrix(iter1,:));
        end

%       curprob = zeros(1,numrxns);
%         if ~exist('curprob','var')
%             curprob = zeros(1,numrxns);
%             for iter1 = 1:numrxns        
%                 curprob(iter1) = rates(iter1)*prod(spc .^ rxn_matrix(iter1,:));
%             end
%         else
%             % Determine which rates need updating
%             updated_prob = logical(rxn_matrix(:,logical(selected_rxn)));
%             % update probabilites of each reaction occuring ( curprob )
%             curprob(updated_prob) = rates(updated_prob)*prod(spc .^ rxn_matrix(updated_prob,:));
%         end
                
%         try
%             % Determine which rates need updating
%             updated_prob = logical(rxn_matrix(:,logical(selected_rxn)));
%             % update probabilites of each reaction occuring ( curprob )
%             curprob(updated_prob) = rates(updated_prob)*prod(spc .^ rxn_matrix(updated_prob,:));
%         catch err1
%             % Catch error when run the first time
%             if ~exist('curprob','var')
%                 curprob = zeros(1,numrxns);
%                 for iter1 = 1:numrxns        
%                     curprob(iter1) = rates(iter1)*prod(spc .^ rxn_matrix(iter1,:));
%                 end
%             else
%                 throw(err1)
%             end
%         end

        % waiting time randomly selected from exponetial distribution with
        % lambda = sum of all reaction probabilites
        delta_t = -log(rnums(1,curstep,runnum)) ./ sum(curprob);

        % Filter reactions with 0 probability to avoid errors
        valid_curprobs = curprob(curprob>0);
        valid_state_changes = state_change_matrix(curprob>0,:);
        % Determine intervals with size proportional to the
        % probability of that reaction occuring
        intvls = cumsum(valid_curprobs) ./ sum(valid_curprobs);
        % use a random number to select the interval
        selected_ind = find(intvls > rnums(2,curstep,runnum),1);
        selected_rxn = valid_state_changes(selected_ind,:);

        % update species
        spcout(:,curstep,runnum)=spc;
        spc = spc + selected_rxn;

        % update time
        tout(1,curstep,runnum)=t;
        t = t + delta_t;        

    end    
    toc
    
%     plot(tout(1,:,runnum),spcout(1,:,runnum),'k')
end

% figure
% plot(tout(1,:,1),spcout(1,:,1))

% inaccurate measure of mean concentration because of the different 
% delta t values. Higher values update more often and so they
% are sampled more often. What is the equation for this biased
% distribution?
% figure; hist(mean(spcout,3),40)

timedur=[diff(tout(1,:,:),1) zeros(1,1,numrpts)];
% figure; hist(permute(timedur,[2,3,1]),40)

% permute(sum(timedur,2),[3,1,1])
% figure; hist(permute(sum(timedur,2),[2,3,1]),20)
% n.*timedur.*repmat(sum(timedur,2),1,maxstep)


end
