function [ varargout ] = partition_functions( cI_2,cro_2,RNAP,cII, varargin )
%Calculates the propensity factors for transcript initiation at the 4 promoters
%   AF20150827
%   Based on 
%   Shea, M. A. & Ackers, G. K. The OR control system of bacteriophage lambda. A physical-chemical model for gene regulation.�J. Mol. Biol.�181,�211�30 (1985).

% Gas constant
R = 1.987*10^-3; %kcal / mol K

if nargin>4
    T = varargin{1};
else
%     T = 30+273.15;
    T = 37+273.15;
end

if nargin>5   
    new_delG=varargin{2};
%     new_delG_ind=varargin{3};
else
    new_delG=[];
%     new_delG_ind=0;
end

% if nargin>6
%     rand_num=varargin{2};
% else
%     rand_num= rand(3,1);
% end

%% Determine transcription initiation rate of PR / PRM promoters
% delta G of each of the possible 40 states
delG_PR=[0;-11.7;-10.1;-10.1;-10.8;-10.8;-12.1;-11.5;-12.5;...
        -23.7;-21.8;-22.2;-21.6;-22.9;-22.9;-24;-22.5;-20.9;...
        -20.9;-23.8;-20.9;-22.2;-22.6;-21.6;-23.2;-24.6;-22.3;...
        -22.3;-33.8;-33.7;-35.8;-32.6;-33;-31.7;-33;-34.6;-35.2;...
        -33.1;-34;-32.4];
if ~isempty(new_delG)
    delG_PR=new_delG;
end
% Number of each element bound in each state
spc_bound_PR=[0,0,0;...
    1,0,0;1,0,0;1,0,0;...
    0,1,0;0,1,0;0,1,0;...
    0,0,1;0,0,1;...
    2,0,0;2,0,0;2,0,0;...
    0,2,0;0,2,0;0,2,0;...
    0,0,2;...
    1,1,0;1,1,0;1,1,0;1,1,0;1,1,0;1,1,0;...
    1,0,1;1,0,1;1,0,1;...
    0,1,1;0,1,1;0,1,1;...
    3,0,0;0,3,0;...
    2,1,0;2,1,0;2,1,0;...
    1,2,0;1,2,0;1,2,0;...
    2,0,1;0,2,1;...
    1,1,1;1,1,1];
    
% partial partition function
Z_PR= exp(-delG_PR./(R*T)).*...
     ((cI_2) .^ spc_bound_PR(:,1)).*...
     ((cro_2) .^ spc_bound_PR(:,2)).*...
     ((RNAP) .^ spc_bound_PR(:,3));
% Z / Ztotal for each state
Prob_PR= Z_PR/sum(Z_PR);  

% DEBUG
% ind=3;
% exp(-delG_PR(ind)./(R*T)).*...
%      ((cI_2) .^ spc_bound_PR(ind,1)).*...
%      ((cro_2) .^ spc_bound_PR(ind,2)).*...
%      ((RNAP) .^ spc_bound_PR(ind,3))


% find index of the state biased randomly chosen
% state_PR=find(cumsum(Prob_PR) > rand_num(1),1);

% rate of transcrpion activation of each state
% kPRM=[0;0;0;0;0;0;0;0.00100;0;0;0;0;0;0;0;0.00100;0;0;0;0;0;0;0;...
%     0.0110;0.00100;0;0.00100;0.00100;0;0;0;0;0;0;0;0;0.0110;...
%     0.00100;0.00100;0.0110];
% kPR=[0;0;0;0;0;0;0;0;0.0140;0;0;0;0;0;0;0.0140;0;0;0;0;0;0;0.0140;...
%     0;0;0.0140;0;0;0;0;0;0;0;0;0;0;0;0;0;0];

% output the rates of the chosen state
% kPRMout = kPRM(state_PR);
% kPRout = kPR(state_PR);




%% Determine transcription initiation rate of PL promoter

delG_PL=[0;-10.9;-12.1;-11.7;...
        -10.1;-12.5;-22.9;-20.9;...
        -22.8;-23.7];

spc_bound_PL=[0,0,0;...
              0,1,0;0,1,0;...
              1,0,0;1,0,0;...
              0,0,1;...
              0,2,0;...
              1,1,0;1,1,0;...
              2,0,0];

Z_PL= exp(-delG_PL./(R*T)).*...
     ((cI_2) .^ spc_bound_PL(:,1)).*...
     ((cro_2) .^ spc_bound_PL(:,2)).*...
     ((RNAP) .^ spc_bound_PL(:,3));
Prob_PL= Z_PL/sum(Z_PL); 

% state_PL=find(cumsum(Prob_PL) > rand_num(2),1);
% 
% kPL=[0;0;0;0;0;0.011;0;0;0;0];

% kPLout=kPL(state_PL);


%% Determine transcription initiation rate of PRE promoter

delG_PRE=[0;-9.9;-9.7;-21.5];

spc_bound_PRE= [0,0;0,1;...
                1,0;1,1];

Z_PRE= exp(-delG_PRE./(R*T)).*...
     ((cII) .^ spc_bound_PRE(:,1)).*...
     ((RNAP) .^ spc_bound_PRE(:,2));
Prob_PRE= Z_PRE/sum(Z_PRE); 

% state_PRE=find(cumsum(Prob_PRE) > rand_num(3),1);
% 
% kPRE=[0;0.00004;0;0.015];

% kPREout=kPRE(state_PRE);

%% Outputs

% kPout(1) = kPR(state_PR);
% kPout(2) = kPRE(state_PRE);
% kPout(3) = kPRM(state_PR);
% kPout(4) = kPL(state_PL);

varargout{1} = Prob_PR;
varargout{2} = Prob_PRE;
varargout{3} = Prob_PL;


end

