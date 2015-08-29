function [ varargout ] = LambdaModel( varargin )
%Gillespie model of pahge lambda infection of e. coli
% 20150828
% Andrew Frink
% Golding Lab
% based on
% Arkin, A., Ross, J. & McAdams, H. H. 
% Stochastic kinetic analysis of developmental pathway bifurcation in phage
% lambda-infected Escherichia coli cells. Genetics 149, 1633–48 (1998)
close all
%% Initial conditions and parameters
% Number of cells to simulate
numruns = 1;

% Set random seed for reproducibility
rng(159265357)

% temp
T = 30+273.15; 

% Initial Volume (L)
Vin= 10^-15;

% size (nt) of the region blocked by RNAP
nt_RNAPblocked=15; 
% Komissarova, N., and M. Kashlev. 1998. Functional topography of 
% nascent RNA in elongation intermediates of RNA polymerase. 
% Proc Natl Acad Sci USA 95:14699-704.

% size (nt) of the region blocked by Ribosome
nt_riboblocked=15;
% For now just using RNAP occulsion.
% Need to find a reference for this.

% Ending condition of the simulation
    % Final Volume (L)
    Vmax= 2*10^-15;
    % Final Time (s)
%     tmax=35*60;
    % Maximum number of steps
%     stepmax=10000;

% initial MOI of cells
n.phage=1;

% No phage proteins before infection
n.cro=0;
n.cro_2=0;
n.cI=0;
n.cI_2=0;
n.cII=0;
n.cIII=0;
n.N=0;
n.P1cII=0;
n.P1cIII=0;
n.P2cII=0;
n.P2cIII=0;

% No phage mRNA or bound ribosomes
PRmRNA.pos=[];
PRmRNA.ribopos=[];
PRmRNA.N=0; % is N bound?
PRmRNA.T=0; % is the mRNA terminated?
PRmRNA.km=zeros(1,5);

PREmRNA.pos=[];
PREmRNA.ribopos=[];
PREmRNA.T=0;
PREmRNA.km=zeros(1,5);

PRMmRNA.pos=[];
PRMmRNA.ribopos=[];
PRMmRNA.T=0;
PRMmRNA.km=zeros(1,5);

PLmRNA.pos=[];
PLmRNA.ribopos=[];
PLmRNA.N=0; % is N bound?
PLmRNA.T=0; % is the mRNA terminated?
PLmRNA.km=zeros(1,5);

% Reaction rate constants (see table)
k= [0.0007;0.050;0.50;0.0025;0.050;0.50;0.00231;0.010;0.010;0.0020;0.010;0.0010;0.00010;0.00025;0.065;0.60;0.010;0.010;0.0010;NaN;NaN;30;5;0.145;0.10;30;15;15;30;NaN;5;25;30;0.0020;100;0.20];

% Growth rate constant from Arkin paper
% They claim it is to achieve 10^-15 L linear growth in 35 min
% BUT 10^-15 L / (35*60) s = 4.76*10^-19
k0=4.76*10^-18; 

% Fixed concentrations (M) of host proteins
c.RNAP=30*10^-9;
c.ribo=500*10^-9;
c.P1=35*10^-9;
c.P2=140*10^-9;

% Constants
A = 6.0221413*10^23; %avagodro's
R = 1.987*10^-3; % gas constant (kcal / mol)

%% Create some molecules for DeBugging

n.cro=10;
n.cro_2=20;
n.cI=10;
n.cI_2=20;
n.cII=10;
n.cIII=10;
n.N=10;
n.P1cII=10;
n.P1cIII=10;
n.P2cII=10;
n.P2cIII=10;

PRmRNA.pos=[38023+1000];
PRmRNA.ribopos=[38023+100 38023+300];
PRmRNA.N=0; % is N bound?
PRmRNA.T=0; % is the mRNA terminated?
PRmRNA.km=zeros(1,5);

PRmRNA(2).pos=[38265];
PRmRNA(2).ribopos=[38023+100];
PRmRNA(2).N=0; % is N bound?
PRmRNA(2).T=0; % is the mRNA terminated?
PRmRNA(2).km=zeros(1,5);

PREmRNA.pos=[38343-1000];
PREmRNA.ribopos=[38343-300 38343-400 38343-700];
PREmRNA.km=zeros(1,5);

PRMmRNA.pos=[];
PRMmRNA.ribopos=[];
PRMmRNA.km=zeros(1,5);

PLmRNA.pos=[];
PLmRNA.ribopos=[];
PLmRNA.N=0; % is N bound?
PLmRNA.T=0; % is the mRNA terminated?
PLmRNA.km=zeros(1,5);

%% Genetic locations 
% Rightward genes
PR_trans = [38023 40624]; % ->
% tr0 terminator 38135
cro_span = [38041,38241]; % ->
NutR= 38265;
% NutR 38265-38281
TR1 = 38337;
% tr1a terminator 38315
% tr1b terminator 38337
% tr1c terminator 38370
cII_span = [38360,38653]; % ->

% Leftward genes
PRE_trans = [38343 35798]; % <-
% anticro_span = [38241,38041]; % <-
PRM_trans = [37958 35798]; % <-
cI_span = [37940,37227]; % <-

PL_trans = [35600 33100]; % <-
NutL = 35518;
% NutL 35534-35518 % <-
N_span = [35582,34560]; % <-
TL1 = 33930; % technically TL2 but I needed TX to stop more than 15 bp from N
% tl1 terminator 34560
% tl2a terminator 33930
% tl2b terminator 33494
cIII_span = [33463,33299]; % <-
% tl2c terminator 33141
% tl2d terminator 33100
% tl3 terminator 31262
% ti terminator 27538
% j1 terminator 18671

%% Simulation Start
tic
for runnum = 1:numruns
    t=0;
    curstep=0;
    propen=zeros(1,40);
    f_continue=1;
    %%
    while f_continue
    % update volume based on growth rate
        V=Vin+k0*t;
       
%     % update housekeeping molecules
%         % Number of RNAP required for 30 nM concentration
%         n.RNAP= c.RNAP*A*V;
%         % Number of ribo required for 500 nM concentration
%         n.ribo= c.ribo*A*V;
%         % Number of P1 required for 35 nM concentration
%         n.P1= c.P1*A*V;
%         % Number of RNAP required for 30 nM concentration
%         n.P2= c.P2*A*V;

    % update concentrations for rate calculations
        c.cro=	n.cro/A/V;
        c.cro_2= n.cro_2/A/V;
        c.cI=   n.cI/A/V;
        c.cI_2=  n.cI_2/A/V;
        c.cII=  n.cII/A/V;
        c.cIII= n.cIII/A/V;
        c.N=    n.N/A/V;
        c.P1cII= n.P1cII/A/V;
        c.P2cII= n.P2cII/A/V;
        c.P1cIII= n.P1cIII/A/V;
        c.P2cIII= n.P2cIII/A/V;

    % Calc propensity of TX initiation at each promoter
    propen(20:23)=partition_functions(c.cI_2,c.cro_2,c.RNAP,c.cII,T,rand(3,1));
    % Check for RNAP occlusion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

    
    
    
    %% Calc propensity of TX reactions
    % Calc Propensity vector (km) for each mRNA
    % km(1)= DNA+1, km(2)= ribo binds, km(3)= RNAse binds,...
    % km(4) = terminate TX, km(5)= toggle N state    
    
    % Clear vector so that we can use it's size later
    tx_propen=[];    
    % PRmRNA
    if ~isempty(PRmRNA(1).pos)
        for iter1=1:numel(PRmRNA)
            if PRmRNA(iter1).T % Is the mRNA terminated
                PRmRNA(iter1).km([1 4 5])=0;                                                               %#ok<*AGROW>
            elseif PRmRNA(iter1).pos == NutR % at N site
                if PRmRNA(iter1).N          % N bound to RNAP
                    PRmRNA(iter1).km(1) = k(26);
                    PRmRNA(iter1).km(4) = 0;
                    PRmRNA(iter1).km(5) = k(25);                    
                else
                    PRmRNA(iter1).km(1) = k(23);
                    PRmRNA(iter1).km(4) = 0;
                    PRmRNA(iter1).km(5) = k(24)*c.N;                    
                end
            elseif PRmRNA(iter1).pos == TR1; % at TR1
                PRmRNA(iter1).km(1) = k(27);
                PRmRNA(iter1).km(4) = k(28);  
                PRmRNA(iter1).km(5) = 0;                          
            else                             % anywhere else
                PRmRNA(iter1).km(1) = k(22);
                PRmRNA(iter1).km(4) = 0;
                PRmRNA(iter1).km(5) = 0;               
            end
            % Ribo and RNAse RBS binding
            if min(PRmRNA(iter1).ribopos) < PR_trans(1)+nt_riboblocked ... % RBS is blocked by Ribosome
                    ||  PRmRNA(iter1).pos < PR_trans(1)+nt_RNAPblocked % RBS is blocked by RNAP
                PRmRNA(iter1).km(2)=0;
                PRmRNA(iter1).km(3)=0;
            else
                PRmRNA(iter1).km(2)=k(34)*c.ribo;
                PRmRNA(iter1).km(3)=k(36);
            end
        end        
    end
    % PREmRNA
    if ~isempty(PREmRNA(1).pos)
        for iter2=1:numel(PREmRNA)
            if PREmRNA(iter2).T
                PREmRNA(iter2).km([1 4 5])=0;
            else
                PREmRNA(iter2).km(1) = k(22);
                PREmRNA(iter2).km(4) = 0;
                PREmRNA(iter2).km(5) = 0; 
            end
            % Ribo and RNAse RBS binding
            if max(PREmRNA(iter2).ribopos) > PRE_trans(1)-nt_riboblocked ...
                    || PREmRNA(iter2).pos > PRE_trans(1)-nt_RNAPblocked
                PREmRNA(iter2).km(2)=0;
                PREmRNA(iter2).km(3)=0;
            else
                PREmRNA(iter2).km(2)=k(34)*c.ribo;
                PREmRNA(iter2).km(3)=k(36);                
            end
        end
    end
    % PRMmRNA
    if ~isempty(PRMmRNA(1).pos)
        for iter3=1:numel(PRMmRNA)
            if PRMmRNA(iter3).T
                PRMmRNA(iter3).km([1 4 5])=0;
            else
                PRMmRNA(iter3).km(1) = k(22);
                PRMmRNA(iter3).km(4) = 0;
                PRMmRNA(iter3).km(5) = 0; 
            end
            % Ribo and RNAse RBS binding
            if max(PRMmRNA(iter3).ribopos) > PRM_trans(1)-nt_riboblocked ...
                    || PRMmRNA(iter3).pos > PRM_trans(1)-nt_RNAPblocked
                PRMmRNA(iter3).km(2)=0;
                PRMmRNA(iter3).km(3)=0;
            else
                PRMmRNA(iter3).km(2)=k(34)*c.ribo;
                PRMmRNA(iter3).km(3)=k(36);                
            end
        end
    end
    %PLmRNA
    if ~isempty(PLmRNA(1).pos)
        for iter4=1:numel(PLmRNA)
            if PLmRNA(iter4).T % Is the mRNA terminated
                PLmRNA(iter4).km([1 4 5])=0;                                                               %#ok<*AGROW>
            elseif PLmRNA(iter4).pos == NutL % at N site
                if PLmRNA(iter4).N          % N bound to RNAP
                    PLmRNA(iter4).km(1) = k(26);
                    PLmRNA(iter4).km(4) = 0;
                    PLmRNA(iter4).km(5) = k(25);                    
                else
                    PLmRNA(iter4).km(1) = k(23);
                    PLmRNA(iter4).km(4) = 0;
                    PLmRNA(iter4).km(5) = k(24)*c.N;                    
                end
            elseif PLmRNA(iter4).pos == TL1; % at TL1
                PLmRNA(iter4).km(1) = k(27);
                PLmRNA(iter4).km(4) = k(28);  
                PLmRNA(iter4).km(5) = 0;                          
            else                             % anywhere else
                PLmRNA(iter4).km(1) = k(22);
                PLmRNA(iter4).km(4) = 0;
                PLmRNA(iter4).km(5) = 0;               
            end
            % Ribo and RNAse RBS binding
            if min(PLmRNA(iter4).ribopos) > PL_trans(1)-nt_riboblocked ... % RBS is blocked by Ribosome
                    ||  PLmRNA(iter4).pos > PL_trans(1)-nt_RNAPblocked % RBS is blocked by RNAP
                PLmRNA(iter4).km(2)=0;
                PLmRNA(iter4).km(3)=0;
            else
                PLmRNA(iter4).km(2)=k(34)*c.ribo;
                PLmRNA(iter4).km(3)=k(36);
            end
        end        
    end
        
    tx_propen=[PRmRNA(:).km PREmRNA(:).km PRMmRNA(:).km PLmRNA(:).km];
    propen(24:23+numel(tx_propen)) = tx_propen;
    
    
    
    %% Calc propensity of TL reactions
    
    tl_propen=[];
    
            
    %% Calc propensity of Protein reactions
    % See table for details on these reactions
    % cI
    propen(1)=k(1)*c.cI;
    propen(2)=k(2)*c.cI^2;
    propen(3)=k(3)*c.cI_2;
    % cro
    propen(4)=k(4)*c.cro;
    propen(5)=k(5)*c.cro^2;    
    propen(6)=k(6)*c.cro_2;
    % N deg
    propen(7)=k(7)*c.N;
    % P1 deg pathway
    propen(8)=k(8)*c.cII*c.P1;
    propen(9)=k(9)*c.P1cII;
    propen(10)=k(10)*c.P1cII;
    propen(11)=k(11)*c.cIII*c.P1;
    propen(12)=k(12)*c.P1cIII;
    propen(13)=k(13)*c.P1cIII;
    % P2 deg pathway
    propen(14)=k(8)*c.cII*c.P2;
    propen(15)=k(9)*c.P2cII;
    propen(16)=k(10)*c.P2cII;
    propen(17)=k(11)*c.cIII*c.P2;
    propen(18)=k(12)*c.P2cIII;
    propen(19)=k(13)*c.P2cIII;  
    
    
    
    %% Choose reaction 
    rannum=rand(1);
    rxn_ind=[1: 23+numel(tx_propen)+numel(tl_propen)];
    v_rxn_ind=rxn_ind(propen>0);
    v_propen= propen(propen>0);
    v_ind=find(cumsum(v_propen)./ sum(v_propen) > rannum,1);
    selected_rxn = v_rxn_ind(v_ind)
    %% Execute reaction
    
%     tx_propen(selected_rxn-24)
    
    
    
    %% Choose timestep
    
    delta_t = -log(rand(1)) /sum(propen);
    
    f_continue=0;
        
    end
    toc
end





