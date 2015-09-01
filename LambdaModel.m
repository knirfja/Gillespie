function [ varargout ] = LambdaModel( varargin )
%Gillespie model of pahge lambda infection of e. coli V1.0
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
Tmp = 30+273.15; 

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
    tmax=10*60;
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
PRmRNA.D=0; % is the mRNA being degraded by RNAse?
PRmRNA.km=zeros(1,5);
PRmRNA.kr=[];

PREmRNA.pos=[];
PREmRNA.ribopos=[];
PREmRNA.T=0;
PREmRNA.D=0; % is the mRNA being degraded by RNAse?
PREmRNA.km=zeros(1,5);
PREmRNA.kr=[];

PRMmRNA.pos=[];
PRMmRNA.ribopos=[];
PRMmRNA.T=0;
PRMmRNA.D=0; % is the mRNA being degraded by RNAse?
PRMmRNA.km=zeros(1,5);
PRMmRNA.kr=[];

PLmRNA.pos=[];
PLmRNA.ribopos=[];
PLmRNA.N=0; % is N bound?
PLmRNA.T=0; % is the mRNA terminated?.
PLmRNA.D=0; % is the mRNA being degraded by RNAse?
PLmRNA.km=zeros(1,5);
PLmRNA.kr=[];

% Reaction rate constants (see table)
k= [0.0007;0.050;0.50;0.0025;0.050;0.50;0.00231;0.010;0.010;0.0020;...
    0.010;0.0010;0.00010;0.00025;0.065;0.60;0.010;0.010;0.0010;NaN;...
    NaN;30;5;0.145;0.10;30;15;15;30;NaN;5;25;30;0.0020;100;0.20];

% reaction rates of each promoter state
kPR=[0;0;0;0;0;0;0;0;0.0140;0;0;0;0;0;0;0.0140;0;0;0;0;0;0;0.0140;...
    0;0;0.0140;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
kPRM=[0;0;0;0;0;0;0;0.00100;0;0;0;0;0;0;0;0.00100;0;0;0;0;0;0;0;...
    0.0110;0.00100;0;0.00100;0.00100;0;0;0;0;0;0;0;0;0.0110;...
    0.00100;0.00100;0.0110];
kPRE=[0;0.00004;0;0.015];
kPL=[0;0;0;0;0;0.011;0;0;0;0];

% Growth rate constant from Arkin paper
% They claim it is to achieve 10^-15 L linear growth in 35 min
% BUT 10^-15 L / (35*60) s = 4.76*10^-19
% k0=4.76*10^-18; 
k0=4.76*10^-19; 

% Fixed concentrations (M) of host proteins
c.RNAP=30*10^-9;
c.ribo=500*10^-9;
c.P1=35*10^-9;
c.P2=140*10^-9;

% Constants
A = 6.0221413*10^23; %avagodro's
R = 1.987*10^-3; % gas constant (kcal / mol)


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
PRmRNA.D=0; % is the mRNA being degraded by RNAse?
PRmRNA.km=zeros(1,5);

PRmRNA(2).pos=[38265];
PRmRNA(2).ribopos=[38023+100];
PRmRNA(2).N=0; % is N bound?
PRmRNA(2).T=0; % is the mRNA terminated?
PRmRNA(2).D=0; % is the mRNA being degraded by RNAse?
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

%% Simulation Start
timeout=zeros(numruns,100000);
maxstep=1;
tic
for runnum = 1:numruns
    disp(['Run ' num2str(runnum)])
    t=0;
    curstep=0;
    V=Vin;
    propen=zeros(1,40);
    f_continue=1;    
    while f_continue
        curstep=curstep+1;
        if curstep == 3223
            disp('')
        end

       
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

        %% Propensity of TX initiation at each promoter
    %     propen(20:23)=partition_functions(c.cI_2,c.cro_2,c.RNAP,c.cII,T,rand(3,1));
        [Prob_PR,Prob_PRE,Prob_PL]=partition_functions(c.cI_2,c.cro_2,c.RNAP,c.cII,Tmp,rand(3,1));

        state_PR=find(cumsum(Prob_PR) > rand(1),1);
        state_PRE=find(cumsum(Prob_PRE) > rand(1),1);
        state_PL=find(cumsum(Prob_PL) > rand(1),1);

        if min([PRmRNA.pos]) < PR_trans(1)+nt_RNAPblocked % Check for RNAP occlusion of the start site.
            propen(20)=0;
        else
            propen(20)=kPR(state_PR);

        end

        if max([PREmRNA.pos]) > PRE_trans(1)-nt_RNAPblocked
            propen(21)=0;
        else
            propen(21)=kPRE(state_PRE);
        end

        if max([PRMmRNA.pos]) > PRM_trans(1)-nt_RNAPblocked
            propen(22)=0;
        else
            propen(22)= kPRM(state_PR);
        end

        if max([PLmRNA.pos]) > PL_trans(1)-nt_RNAPblocked
            propen(23)=0;
        else
        propen(23)= kPL(state_PL);  
        end

        %% Propensity of TX reactions
        % Propensity vector (km) for each mRNA
        % km(1)= DNA+1
        % km(2)= ribo binds
        % km(3)= RNAse binds,...
        % km(4) = terminate TX
        % km(5)= toggle N state    

        % Clear vector so that we can use it's size later
        tx_propen=[];     %#ok<NASGU>
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
                if PRmRNA(iter1).pos < PR_trans(1)+nt_RNAPblocked % RBS is blocked by RNAP
                    PRmRNA(iter1).km(2)=0;
                    PRmRNA(iter1).km(3)=0;
                elseif min(PRmRNA(iter1).ribopos) < PR_trans(1)+nt_riboblocked % RBS is blocked by a Ribosome
                    PRmRNA(iter1).km(2)=0;
                    PRmRNA(iter1).km(3)=0;
                elseif PRmRNA(iter1).D                            % RNAse has already degraded the RBS 
                    PRmRNA(iter1).km(2)=0;
                    PRmRNA(iter1).km(3)=0;
                else                                              % RBS is free to bind
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
                if PREmRNA(iter2).pos < PRE_trans(1)-nt_RNAPblocked % RBS is blocked by RNAP
                    PREmRNA(iter2).km(2)=0;
                    PREmRNA(iter2).km(3)=0;
                elseif min(PREmRNA(iter2).ribopos) < PRE_trans(1)-nt_riboblocked % RBS is blocked by a Ribosome
                    PREmRNA(iter2).km(2)=0;
                    PREmRNA(iter2).km(3)=0;
                elseif PREmRNA(iter2).D                            % RNAse has already degraded the RBS 
                    PREmRNA(iter2).km(2)=0;
                    PREmRNA(iter2).km(3)=0;
                else                                              % RBS is free to bind
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
                if PRMmRNA(iter3).pos > PRM_trans(1)-nt_RNAPblocked
                    PRMmRNA(iter3).km(2)=0;
                    PRMmRNA(iter3).km(3)=0;
                elseif PRMmRNA(iter3).D 
                    PRMmRNA(iter3).km(2)=0;
                    PRMmRNA(iter3).km(3)=0;
                elseif max(PRMmRNA(iter3).ribopos) > PRM_trans(1)-nt_riboblocked
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
                if PLmRNA(iter4).pos > PL_trans(1)-nt_RNAPblocked % RBS is blocked by RNAP
                    PLmRNA(iter4).km(2)=0;
                    PLmRNA(iter4).km(3)=0;
                elseif PLmRNA(iter4).D % RBS is degraded
                    PLmRNA(iter4).km(2)=0;
                    PLmRNA(iter4).km(3)=0;
                elseif min(PLmRNA(iter4).ribopos) > PL_trans(1)-nt_riboblocked
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

        %% Propensity of TL reactions   
        if ~isempty(PRmRNA(1).pos) % If there are mRNA
            for iter5=1:numel(PRmRNA) % For each mRNA
                if ~isempty(PRmRNA(iter5).ribopos) % If there are bound ribo
                    n_ribo=numel([PRmRNA(iter5).ribopos]);
                    PRmRNA(iter5).kr=ones(1,n_ribo).*k(35);
                    % Can't get too close to each other
                    if n_ribo>1
                        PRmRNA(iter5).kr(find(-diff(PRmRNA(iter5).ribopos) < nt_riboblocked)+1)=0;
                    end
                    % Can't get too close to the RNAP
                    if PRmRNA(iter5).ribopos(1) > PRmRNA(iter5).pos - nt_RNAPblocked 
                        PRmRNA(iter5).kr(1)=0;
                    end           
                end
            end
        end
        
        if ~isempty(PREmRNA(1).pos)
            for iter6=1:numel(PREmRNA)
                if ~isempty(PREmRNA(iter6).ribopos)
                    n_ribo=numel([PREmRNA(iter6).ribopos]);
                    PREmRNA(iter6).kr=ones(1,n_ribo).*k(35);
                    % Can't get too close to each other
                    if n_ribo>1
                        PREmRNA(iter6).kr(diff(PREmRNA(iter6).ribopos) < nt_riboblocked)=0;
                    end
                    % Can't get too close to the RNAP
                    if PREmRNA(iter6).ribopos(1) < PREmRNA(iter6).pos + nt_RNAPblocked 
                        PREmRNA(iter6).kr(1)=0;
                    end           
                end
            end
        end
        
        if ~isempty(PRMmRNA(1).pos)
            for iter7=1:numel(PRMmRNA)
                if ~isempty(PRMmRNA(iter7).ribopos)
                    n_ribo=numel([PRMmRNA(iter7).ribopos]);
                    PRMmRNA(iter7).kr=ones(1,n_ribo).*k(35);
                    % Can't get too close to each other
                    if n_ribo>1
                        PRMmRNA(iter7).kr(diff(PRMmRNA(iter7).ribopos) < nt_riboblocked)=0;
                    end
                    % Can't get too close to the RNAP
                    if PRMmRNA(iter7).ribopos(1) < PRMmRNA(iter7).pos + nt_RNAPblocked 
                        PRMmRNA(iter7).kr(1)=0;
                    end           
                end
            end
        end
        
        if ~isempty(PLmRNA(1).pos)
            for iter8=1:numel(PLmRNA)
                if ~isempty(PLmRNA(iter8).ribopos)
                    n_ribo=numel([PLmRNA(iter8).ribopos]);
                    PLmRNA(iter8).kr=ones(1,n_ribo).*k(35);
                    % Can't get too close to each other
                    if n_ribo>1
                        PLmRNA(iter8).kr(diff(PLmRNA(iter8).ribopos) < nt_riboblocked)=0;
                    end
                    % Can't get too close to the RNAP
                    if PLmRNA(iter8).ribopos(1) < PLmRNA(iter8).pos + nt_RNAPblocked 
                        PLmRNA(iter8).kr(1)=0;
                    end           
                end
            end
        end
                        
        tl_propen=[PRmRNA(:).kr PREmRNA(:).kr PRMmRNA(:).kr PLmRNA(:).kr];
        propen(24+numel(tx_propen):23+numel(tx_propen)+numel(tl_propen)) = tl_propen;
        %% Propensity of Protein reactions
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
        rxn_ind=1: 23+numel(tx_propen)+numel(tl_propen);
        v_rxn_ind=rxn_ind(propen>0);
        v_propen= propen(propen>0);
        v_ind=find(cumsum(v_propen)./ sum(v_propen) > rannum,1);
        selected_rxn = v_rxn_ind(v_ind);
        %% Execute reaction    
        switch selected_rxn
            %% TX reactions
            case num2cell(24:23+numel(tx_propen))
                % A mRNA based reaction

                tx_ind(1) = numel([PRmRNA(:).km]);
                tx_ind(2) = numel([PREmRNA(:).km]);
                tx_ind(3) = numel([PRMmRNA(:).km]);
                tx_ind(4) = numel([PLmRNA(:).km]);
                tx_ind = cumsum(tx_ind);
                
                switch find(tx_ind>=selected_rxn-23,1)
                    case 1 %PR
                        m_ind=selected_rxn-23;
                        switch mod(selected_rxn-23,5)    
                            case 1
                                % km(1)= DNA+1                        
                                %                         disp('TX elongation')                            
                                PRmRNA(m_ind).pos = PRmRNA(m_ind).pos+1;
                                if PRmRNA(m_ind).pos >= PR_trans(2) % Check for end of transcript
                                    PRmRNA(m_ind).T = 1; 
                                end
                            case 2
                                % km(2)= ribo binds
                            %                         disp('ribosome binding')
                                PRmRNA(m_ind).ribopos = [PRmRNA(m_ind).ribopos PR_trans(1)]; % approximating that RBS is located at transcriptional start
                            case 3
                                % km(3)= RNAse binds
                            %                         disp('RNAse binding')
                                PRmRNA(m_ind).D= 1;                      
                            case 4
                                % km(4) = terminate TX
                            %                         disp('TX termination')
                                PRmRNA(m_ind).T= 1;
                            case 0
                                % km(5)= toggle N state 
                            %                         disp('toggle N state')
                                PRmRNA(m_ind).N= ~PRmRNA(m_ind).N;
                            otherwise
                                disp('Error, reaction not found')
                        end                        
                    case 2 %PRE
                        m_ind=selected_rxn-23-tx_ind(1);
                        switch mod(selected_rxn-23,5)    
                            case 1
                                % km(1)= DNA+1                        
                                %                         disp('TX elongation')                            
                                PREmRNA(m_ind).pos = PREmRNA(m_ind).pos-1; % reverse
                                if PREmRNA(m_ind).pos <= PRE_trans(2) % Check for end of transcript
                                    PREmRNA(m_ind).T = 1; 
                                end
                            case 2
                                % km(2)= ribo binds
                            %                         disp('ribosome binding')
                                PREmRNA(m_ind).ribopos = [PREmRNA(m_ind).ribopos PRE_trans(1)]; % approximating that RBS is located at transcriptional start
                            case 3
                                % km(3)= RNAse binds
                            %                         disp('RNAse binding')
                                PREmRNA(m_ind).D= 1;                      
                            case 4
                                % km(4) = terminate TX
                            %                         disp('TX termination')
                                PREmRNA(m_ind).T= 1;
                            case 0
                                % km(5)= toggle N state 
                            %                         disp('toggle N state')
                                PREmRNA(m_ind).N= ~PREmRNA(m_ind).N;
                            otherwise
                                disp('Error, reaction not found')
                        end
                    case 3 %PRM
                        m_ind=selected_rxn-23-tx_ind(2);
                        switch mod(selected_rxn-23,5)    
                            case 1
                                % km(1)= DNA+1                        
                                %                         disp('TX elongation')                            
                                PRMmRNA(m_ind).pos = PRMmRNA(m_ind).pos-1; % reverse
                                if PRMmRNA(m_ind).pos <= PRM_trans(2) % Check for end of transcript
                                    PRMmRNA(m_ind).T = 1; 
                                end
                            case 2
                                % km(2)= ribo binds
                            %                         disp('ribosome binding')
                                PRMmRNA(m_ind).ribopos = [PRMmRNA(m_ind).ribopos PRM_trans(1)]; % approximating that RBS is located at transcriptional start
                            case 3
                                % km(3)= RNAse binds
                            %                         disp('RNAse binding')
                                PRMmRNA(m_ind).D= 1;                      
                            case 4
                                % km(4) = terminate TX
                            %                         disp('TX termination')
                                PRMmRNA(m_ind).T= 1;
                            case 0
                                % km(5)= toggle N state 
                            %                         disp('toggle N state')
                                PRMmRNA(m_ind).N= ~PRMmRNA(m_ind).N;
                            otherwise
                                disp('Error, reaction not found')
                        end
                    case 4 %PL   
                        m_ind=selected_rxn-23-tx_ind(3);
                        switch mod(selected_rxn-23,5)    
                            case 1
                                % km(1)= DNA+1                        
                                %                         disp('TX elongation')                            
                                PLmRNA(m_ind).pos = PLmRNA(m_ind).pos-1; % reverse
                                if PLmRNA(m_ind).pos <= PL_trans(2) % Check for end of transcript
                                    PLmRNA(m_ind).T = 1; 
                                end
                            case 2
                                % km(2)= ribo binds
                            %                         disp('ribosome binding')
                                PLmRNA(m_ind).ribopos = [PLmRNA(m_ind).ribopos PL_trans(1)]; % approximating that RBS is located at transcriptional start
                            case 3
                                % km(3)= RNAse binds
                            %                         disp('RNAse binding')
                                PLmRNA(m_ind).D= 1;                      
                            case 4
                                % km(4) = terminate TX
                            %                         disp('TX termination')
                                PLmRNA(m_ind).T= 1;
                            case 0
                                % km(5)= toggle N state 
                            %                         disp('toggle N state')
                                PLmRNA(m_ind).N= ~PLmRNA(m_ind).N;
                            otherwise
                                disp('Error, reaction not found')
                        end
                    otherwise
                end
                
                
                    

            %% TL reactions
            case num2cell((24+numel(tx_propen)):(23+numel(tx_propen)+numel(tl_propen))



            %% Protein reactions
            case 1

                n.cI = n.cI - 1;
            case 2
                n.cI = n.cI - 2;
                n.cI_2 = n.cI_2 + 1;
            case 3
                n.cI = n.cI + 2;
                n.cI_2 = n.cI_2 - 1;
            case 4
                n.cro = n.cro - 1;           
            case 5
                n.cro = n.cro - 2;
                n.cro_2 = n.cro_2 + 1;
            case 6
                n.cro = n.cro + 2;
                n.cro_2 = n.cro_2 - 1;
            case 7
                n.N = n.N - 1;            
            case 8
                n.cII = n.cII - 1;
    %             n.P1 = n.P1 - 1;
                n.P1cII = n.P1cII + 1;
            case 9
                n.cII = n.cII + 1;
    %             n.P1 = n.P1 + 1;
                n.P1cII = n.P1cII - 1;
            case 10
    %             n.P1 = n.P1 + 1;  
                n.P1cII = n.P1cII - 1;          
            case 11
                n.cIII = n.cIII - 1;
    %             n.P1 = n.P1 - 1;
                n.P1cIII = n.P1cIII + 1;
            case 12
                n.cIII = n.cIII + 1;
    %             n.P1 = n.P1 + 1;
                n.P1cIII = n.P1cIII - 1;
            case 13
    %             n.P1 = n.P1 + 1;
                n.P1cIII = n.P1cIII - 1;
            case 14
                n.cII = n.cII - 1;
    %             n.P2 = n.P2 - 1;
                n.P2cII = n.P2cII + 1;            
            case 15
                n.cII = n.cII + 1;
    %             n.P2 = n.P2 + 1;
                n.P2cII = n.P2cII - 1;    
            case 16
    %             n.P2 = n.P2 + 1;
                n.P2cII = n.P2cII - 1;             
            case 17
                n.cIII = n.cIII - 1;
    %             n.P2 = n.P2 - 1;
                n.P2cIII = n.P2cIII + 1;  
            case 18
                n.cIII = n.cIII + 1;
    %             n.P2 = n.P2 + 1;
                n.P2cIII = n.P2cIII - 1;  
            case 19
    %             n.P2 = n.P2 + 1;
                n.P2cIII = n.P2cIII - 1;             
            case 20
                %PR
                if isempty(PRmRNA(1).pos)
                    ind1=1;
                else
                    ind1 = numel(PRmRNA)+1;
                end
                PRmRNA(ind1).pos = PR_trans(1);
                PRmRNA(ind1).ribopos=[];
                PRmRNA(ind1).N=0; % is N bound?
                PRmRNA(ind1).T=0; % is the mRNA terminated?
                PRmRNA(ind1).km=zeros(1,5);
            case 21
                %PRE
                if isempty(PREmRNA(1).pos)
                    ind2=1;
                else
                    ind2 = numel(PREmRNA)+1;
                end
                PREmRNA(ind2).pos = PRE_trans(1);
                PREmRNA(ind2).ribopos=[];
                PREmRNA(ind2).N=0; % is N bound?
                PREmRNA(ind2).T=0; % is the mRNA terminated?
                PREmRNA(ind2).km=zeros(1,5);            
            case 22
                %PRM
                if isempty(PRMmRNA(1).pos)
                    ind3=1;
                else
                    ind3 = numel(PRMmRNA)+1;
                end
                PRMmRNA(ind3).pos = PRM_trans(1);
                PRMmRNA(ind3).ribopos=[];
                PRMmRNA(ind3).N=0; % is N bound?
                PRMmRNA(ind3).T=0; % is the mRNA terminated?
                PRMmRNA(ind3).km=zeros(1,5);   
            case 23
                %PL
                if isempty(PLmRNA(1).pos)
                    ind4=1;
                else
                    ind4 = numel(PLmRNA)+1;
                end
                PLmRNA(ind4).pos = PL_trans(1);
                PLmRNA(ind4).ribopos=[];
                PLmRNA(ind4).N=0; % is N bound?
                PLmRNA(ind4).T=0; % is the mRNA terminated?
                PLmRNA(ind4).km=zeros(1,5);   
            otherwise
                disp('Error, reaction not found')
        end

    % kPout(1) = kPR(state_PR);
    % kPout(2) = kPRE(state_PRE);
    % kPout(3) = kPRM(state_PR);
    % kPout(4) = kPL(state_PL);           


    %     tx_propen(selected_rxn-24)


    %% Outputs
        timeout(runnum,curstep) = t;

    %% Choose next timestep

        delta_t = -log(rand(1)) /sum(propen);

        t=t+delta_t;
        % update volume based on growth rate
        V=Vin+k0*t;

        if V > Vmax;
            f_continue=0;
        end

        if t > tmax
            f_continue=0;
        end
    end
    toc
    timeplot=timeout(runnum,1:curstep);
    maxstep=max(maxstep,curstep);
    figure; plot(1:numel(timeplot),timeplot,'.')
    run(runnum).t=timeout(runnum,1:curstep);
end

varargout{1}=run;

end



%% Unused
% %% Aborted execute tx function code
%                 if selected_rxn-23 <= tx_ind(1)
%                     m_name='PR';
%                     m_ind = num2str(selected_rxn-23);
%                 elseif selected_rxn-23 <= tx_ind(2)
%                     m_name='PRE';
%                     m_ind = num2str(selected_rxn-23-tx_ind(1));
%                 elseif selected_rxn-23 <= tx_ind(3)
%                     m_name='PRM';
%                     m_ind = num2str(selected_rxn-23-tx_ind(2));
%                 elseif selected_rxn-23 <= tx_ind(4)
%                     m_name='PL';
%                     m_ind = num2str(selected_rxn-23-tx_ind(3));
%                 else
%                     disp('Error, reaction not found')
%                 end
%                 p_ind=find(tx_ind>=selected_rxn-23,1);
%                 m_ind=floor((selected_rxn-23)/5)+1;
%                 r_ind=mod(selected_rxn-23,5);
                
%                 eval(['tx_rxn(' m_name 'mRNA(' m_ind '),mod(selected_rxn-23,5))'])





