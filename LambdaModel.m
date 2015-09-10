function [ varargout ] = LambdaModel( varargin )
%Gillespie model of pahge lambda infection of e. coli V1.0
% 20150828
% Andrew Frink
% Golding Lab
% based on
% Arkin, A., Ross, J. & McAdams, H. H. 
% Stochastic kinetic analysis of developmental pathway bifurcation in phage
% lambda-infected Escherichia coli cells. Genetics 149, 1633–48 (1998)
% close all

%% Initial conditions and parameters



% Number of cells to simulate
numruns = 1;

%%% Ending condition of the simulation
    % Final Volume (L)
    Vmax= 3*10^-15;
    % Final Time (s)
    tmax=60*50;
    % Maximum number of steps
    stepmax=8*10^5;
    
%%% Set random seed for reproducibility
if nargin > 0
    rng(varargin{1})
else
    rng(159265357)
end

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
PRmRNA=create_mRNA(1);
PREmRNA=create_mRNA(2);
PRMmRNA=create_mRNA(3);
PLmRNA=create_mRNA(4);

% PRmRNA.pos=[];
% PRmRNA.ribopos=[];
% PRmRNA.N=0; % is N bound?
% PRmRNA.T=0; % is the mRNA terminated?
% PRmRNA.Dcro=0; % is the mRNA being degraded by RNAse at cro?
% PRmRNA.DcII=0;  % is the mRNA being degraded by RNAse at cII?
% PRmRNA.km=zeros(1,7);
% PRmRNA.kr=[];
% 
% PREmRNA.pos=[];
% PREmRNA.ribopos=[];
% PREmRNA.T=0;
% PREmRNA.D=0; % is the mRNA being degraded by RNAse?
% PREmRNA.km=zeros(1,7);
% PREmRNA.kr=[];
% 
% PRMmRNA.pos=[];
% PRMmRNA.ribopos=[];
% PRMmRNA.T=0;
% PRMmRNA.D=0; % is the mRNA being degraded by RNAse?
% PRMmRNA.km=zeros(1,7);
% PRMmRNA.kr=[];
% 
% PLmRNA.pos=[];
% PLmRNA.ribopos=[];
% PLmRNA.N=0; % is N bound?
% PLmRNA.T=0; % is the mRNA terminated?.
% PLmRNA.DN=0; % is the mRNA being degraded by RNAse at N?
% PLmRNA.DcIII=0; % is the mRNA being degraded by RNAse at cIII?
% PLmRNA.km=zeros(1,7);
% PLmRNA.kr=[];

% Reaction rate constants (see table)
k= [0.0007;0.050;0.50;0.0025;0.050;0.50;0.00231;0.010;0.010;0.0020;...
    0.010;0.0010;0.00010;0.00025;0.065;0.60;0.010;0.010;0.0010;NaN;...
    NaN;30;5;0.145;0.10;30;15;15;30;NaN;5;25;30;0.002;100;0.20];

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
TL1 = 34560;
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

% n.cro=10;
% n.cro_2=20;
% n.cI=10;
% n.cI_2=20;
% n.cII=10;
% n.cIII=10;
% n.N=10;
% n.P1cII=10;
% n.P1cIII=10;
% n.P2cII=10;
% n.P2cIII=10;
% 
% PRmRNA.pos=[38023+1000];
% PRmRNA.ribopos=[38023+100 38023+300];
% PRmRNA.N=0; % is N bound?
% PRmRNA.T=0; % is the mRNA terminated?
% PRmRNA.Dcro=0; % is the mRNA being degraded by RNAse?
% PRmRNA.DcII=0; % is the mRNA being degraded by RNAse?
% PRmRNA.km=zeros(1,5);
% 
% PRmRNA(2).pos=[38265];
% PRmRNA(2).ribopos=[38023+100];
% PRmRNA(2).N=0; % is N bound?
% PRmRNA(2).T=0; % is the mRNA terminated?
% PRmRNA(2).Dcro=0; % is the mRNA being degraded by RNAse?
% PRmRNA(2).DcII=0; % is the mRNA being degraded by RNAse?
% PRmRNA(2).km=zeros(1,5);
% 
% PREmRNA.pos=[38343-1000];
% PREmRNA.ribopos=[38343-300 38343-400 38343-700];
% PREmRNA.km=zeros(1,5);
% 
% PRMmRNA.pos=[];
% PRMmRNA.ribopos=[];
% PRMmRNA.km=zeros(1,5);
% 
% PLmRNA.pos=[];
% PLmRNA.ribopos=[];
% PLmRNA.N=0; % is N bound?
% PLmRNA.T=0; % is the mRNA terminated?
% PLmRNA.km=zeros(1,5);

%% Simulation Start
timeout=zeros(1,stepmax);
Vout=zeros(1,stepmax);
cI_2out = zeros(1,stepmax);
cro_2out = zeros(1,stepmax);
PRmRNAout = zeros(1,stepmax);
PREmRNAout = zeros(1,stepmax);
PRMmRNAout = zeros(1,stepmax);
PLmRNAout = zeros(1,stepmax);
Nout= zeros(1,stepmax);
cIIout= zeros(1,stepmax);
cIIIout= zeros(1,stepmax);


for runnum = 1:numruns
%     disp(['Run ' num2str(runnum)])
    t=0;
    curstep=0;
    V=Vin;
    propen=zeros(1,40);
    f_continue=1;    
    while f_continue
        curstep=curstep+1;
%         % Debug
%         if curstep == 44551
%             disp('')
%         end

        AV=A*V;
%         % DEBUG
%         if t>60*2
%             disp('')
%         end
    % update housekeeping molecules
        % Number of RNAP required for 30 nM concentration
        n.RNAP= c.RNAP*AV;
        % Number of ribo required for 500 nM concentration
        n.ribo= c.ribo*AV;
        % Number of P1 required for 35 nM concentration
        n.P1= c.P1*AV;
        % Number of RNAP required for 30 nM concentration
        n.P2= c.P2*AV;
        

    % update concentrations for partition function calculations
%         c.cro=	n.cro/AV;
        c.cro_2= n.cro_2/AV;
%         c.cI=   n.cI/AV;
        c.cI_2=  n.cI_2/AV;
        c.cII=  n.cII/AV;
%         c.cIII= n.cIII/AV;
%         c.N=    n.N/AV;
%         c.P1cII= n.P1cII/AV;
%         c.P2cII= n.P2cII/AV;
%         c.P1cIII= n.P1cIII/AV;
%         c.P2cIII= n.P2cIII/AV;
    %%% Store Outputs
        timeout(1,curstep) = t;
        Vout(1,curstep)= V;
        
        cIout(1,curstep) = n.cI;
        croout(1,curstep) = n.cro;
        cI_2out(1,curstep) = n.cI_2;
        cro_2out(1,curstep) = n.cro_2;
        Nout(1,curstep) = n.N;
        cIIout(1,curstep) = n.cII;
        cIIIout(1,curstep) = n.cIII;
        % DEBUG#
        if isempty(PRmRNA)
            disp('')
        end
        if isempty(PRmRNA(1).pos)
            PRmRNAout(1,curstep) = 0;
        else
            PRmRNAout(1,curstep) = numel(PRmRNA);
        end
        if isempty(PREmRNA(1).pos)
            PREmRNAout(1,curstep) = 0;
        else
            PREmRNAout(1,curstep) = numel(PREmRNA);
        end
        if isempty(PRMmRNA(1).pos)
            PRMmRNAout(1,curstep) = 0;
        else
            PRMmRNAout(1,curstep) = numel(PRMmRNA);
        end
        if isempty(PLmRNA(1).pos)
            PLmRNAout(1,curstep) = 0;
        else
            PLmRNAout(1,curstep) = numel(PLmRNA);
        end

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
        % km(3)= RNAse binds
        % km(4) = terminate TX
        % km(5)= toggle N state 
        % km(6)= ribo binds 2nd gene
        % km(7)= RNAse binds 2nd gene
        % 

        % Clear vector so that we can use it's size later
        tx_propen=[]; 
        % PRmRNA
        if ~isempty(PRmRNA(1).pos)
            for iter1=1:numel(PRmRNA)
                if PRmRNA(iter1).T % Is the mRNA terminated
                    PRmRNA(iter1).km([1 4 5])=0;                                                             
                elseif PRmRNA(iter1).pos == NutR % at N site
                    if PRmRNA(iter1).N          % N bound to RNAP
                        PRmRNA(iter1).km(1) = k(26);
                        PRmRNA(iter1).km(4) = 0;
                        PRmRNA(iter1).km(5) = k(25);                    
                    else
                        PRmRNA(iter1).km(1) = k(23);
                        PRmRNA(iter1).km(4) = 0;
                        PRmRNA(iter1).km(5) = k(24)*n.N;                    
                    end
                elseif PRmRNA(iter1).pos == TR1 && ~PRmRNA(iter1).N; % at TR1 and no N
                    PRmRNA(iter1).km(1) = k(27);
                    PRmRNA(iter1).km(4) = k(28);  
                    PRmRNA(iter1).km(5) = 0;                          
                else                             % anywhere else
                    PRmRNA(iter1).km(1) = k(22);
                    PRmRNA(iter1).km(4) = 0;
                    PRmRNA(iter1).km(5) = 0;               
                end                
                % Ribo and RNAse RBS binding
                % cro
                if PRmRNA(iter1).Dcro                            % RNAse has already degraded the cro RBS 
                    PRmRNA(iter1).km(2)=0;
                    PRmRNA(iter1).km(3)=0;
                elseif PRmRNA(iter1).pos < cro_span(1)+nt_RNAPblocked % RBS is not transcribed or is blocked by RNAP
                    PRmRNA(iter1).km(2)=0;
                    PRmRNA(iter1).km(3)=0;
                elseif min(PRmRNA(iter1).ribopos) < cro_span(1)+nt_riboblocked
                    PRmRNA(iter1).km(2)=0;
                    PRmRNA(iter1).km(3)=0;
                else                                              % RBS is free to bind
                    PRmRNA(iter1).km(2)=k(34)*n.ribo;
                    PRmRNA(iter1).km(3)=k(36);
                end                
                % cII
                if PRmRNA(iter1).DcII                            % RNAse has already degraded the cII RBS 
                    PRmRNA(iter1).km(6)=0;
                    PRmRNA(iter1).km(7)=0;
                elseif PRmRNA(iter1).pos < cII_span(1)+nt_RNAPblocked % RBS is not transcribed or is blocked by RNAP
                    PRmRNA(iter1).km(6)=0;
                    PRmRNA(iter1).km(7)=0;
                elseif PRmRNA(iter1).ribopos(PRmRNA(iter1).ribopos > cro_span(2)) < cII_span(1)+nt_riboblocked
                    % RBS is blocked by a Ribosome right of cro and left of
                    % the cII RBS
                    PRmRNA(iter1).km(6)=0;
                    PRmRNA(iter1).km(7)=0;                
                else                                              % RBS is free to bind
                    PRmRNA(iter1).km(6)=k(34)*n.ribo;
                    PRmRNA(iter1).km(7)=k(36);
                end
%                 PRmRNA(iter1).km= PRmRNA(iter1).km./AV; % adjust to cell concentration
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
                % cI
                if PREmRNA(iter2).pos > cI_span(1)-nt_RNAPblocked % RBS is not transcribed or is blocked by RNAP
                    PREmRNA(iter2).km(2)=0;
                    PREmRNA(iter2).km(3)=0;
                elseif max(PREmRNA(iter2).ribopos) > cI_span(1)-nt_riboblocked % RBS is blocked by a Ribosome
                    PREmRNA(iter2).km(2)=0;
                    PREmRNA(iter2).km(3)=0;
                elseif PREmRNA(iter2).D                            % RNAse has already degraded the RBS 
                    PREmRNA(iter2).km(2)=0;
                    PREmRNA(iter2).km(3)=0;
                else                                              % RBS is free to bind
                    PREmRNA(iter2).km(2)=k(34)*n.ribo;
                    PREmRNA(iter2).km(3)=k(36);
                end
%                 PREmRNA(iter2).km= PREmRNA(iter2).km./AV; % adjust to cell concentration
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
                % cI
                if PRMmRNA(iter3).pos > cI_span(1)-nt_RNAPblocked
                    PRMmRNA(iter3).km(2)=0;
                    PRMmRNA(iter3).km(3)=0;
                elseif PRMmRNA(iter3).D 
                    PRMmRNA(iter3).km(2)=0;
                    PRMmRNA(iter3).km(3)=0;
                elseif max(PRMmRNA(iter3).ribopos) > cI_span(1)-nt_riboblocked
                    PRMmRNA(iter3).km(2)=0;
                    PRMmRNA(iter3).km(3)=0;
                else
                    PRMmRNA(iter3).km(2)=k(34)*n.ribo;
                    PRMmRNA(iter3).km(3)=k(36);                
                end
%                 PRMmRNA(iter3).km= PRMmRNA(iter3).km./AV; % adjust to cell concentration
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
                        PLmRNA(iter4).km(5) = k(24)*n.N;                    
                    end
                elseif PLmRNA(iter4).pos == TL1 && ~PLmRNA(iter4).N; % at TL1 and no N
                    PLmRNA(iter4).km(1) = k(31);
                    PLmRNA(iter4).km(4) = k(32);  
                    PLmRNA(iter4).km(5) = 0;                          
                else                             % anywhere else
                    PLmRNA(iter4).km(1) = k(22);
                    PLmRNA(iter4).km(4) = 0;
                    PLmRNA(iter4).km(5) = 0;               
                end
                % Ribo and RNAse RBS binding
                % N
                if PLmRNA(iter4).pos > N_span(1)-nt_RNAPblocked % RBS is blocked by RNAP
                    PLmRNA(iter4).km(2)=0;
                    PLmRNA(iter4).km(3)=0;
                elseif PLmRNA(iter4).DN % RBS is degraded at N
                    PLmRNA(iter4).km(2)=0;
                    PLmRNA(iter4).km(3)=0;

                elseif max(PLmRNA(iter4).ribopos) > N_span(1)-nt_riboblocked
                    PLmRNA(iter4).km(2)=0;
                    PLmRNA(iter4).km(3)=0;
                else
                    PLmRNA(iter4).km(2)=k(34)*n.ribo;
                    PLmRNA(iter4).km(3)=k(36);
                end
                % cIII
                if PLmRNA(iter4).pos > cIII_span(1)-nt_RNAPblocked % RBS is blocked by RNAP
                    PLmRNA(iter4).km(6)=0;
                    PLmRNA(iter4).km(7)=0;
                elseif PLmRNA(iter4).DN % RBS is degraded at N
                    PLmRNA(iter4).km(6)=0;
                    PLmRNA(iter4).km(7)=0;
                elseif PLmRNA(iter4).ribopos(PLmRNA(iter4).ribopos < N_span(2)) > cIII_span(1)-nt_riboblocked
                    % If any ribo are left of N and right of cIII RBS
                    PLmRNA(iter4).km(6)=0;
                    PLmRNA(iter4).km(7)=0;
                else
                    PLmRNA(iter4).km(6)=k(34)*n.ribo;
                    PLmRNA(iter4).km(7)=k(36);
                end  
%                 PLmRNA(iter4).km= PLmRNA(iter4).km./AV; % adjust to cell concentration
            end        
        end
        tx_propen=[PRmRNA(:).km PREmRNA(:).km PRMmRNA(:).km PLmRNA(:).km];
%         tx_propen(tx_propen>0)
        propen(24:23+numel(tx_propen)) = tx_propen;

        %% Propensity of TL reactions   
        % PR transcript ribos        
        if ~isempty(PRmRNA(1).pos) % If there are mRNA
            for iter5=1:numel(PRmRNA) % For each mRNA
                if ~isempty(PRmRNA(iter5).ribopos) % If there are bound ribo
                    n_ribo=numel([PRmRNA(iter5).ribopos]);
                    PRmRNA(iter5).kr=ones(1,n_ribo).*k(35);
                    % Can't get too close to each other
                    if n_ribo>1
                        PRmRNA(iter5).kr(find(-diff(PRmRNA(iter5).ribopos) < nt_riboblocked)+1)=0; % Have to use find(...)+1 to return the correct index for the ribosome that is too close
                    end
                    % Can't get too close to the RNAP
                    if PRmRNA(iter5).ribopos(1) > PRmRNA(iter5).pos - nt_RNAPblocked 
                        PRmRNA(iter5).kr(1)=0;
                    end
                else
                    PRmRNA(iter5).kr=[];           
                end
%                 PRmRNA(iter5).kr= PRmRNA(iter5).kr./AV;
            end
        else
            PRmRNA(1).kr=[];
        end
        
        % PRE transcript ribos
        if ~isempty(PREmRNA(1).pos)
            for iter6=1:numel(PREmRNA)
                if ~isempty(PREmRNA(iter6).ribopos)
                    n_ribo=numel([PREmRNA(iter6).ribopos]);
                    PREmRNA(iter6).kr=ones(1,n_ribo).*k(35);
                    % Can't get too close to each other
                    if n_ribo>1
                        PREmRNA(iter6).kr(diff(PREmRNA(iter6).ribopos) < nt_riboblocked)=0; % reverse gene so diff(...) returns the correct index
                    end
                    % Can't get too close to the RNAP
                    if PREmRNA(iter6).ribopos(1) < PREmRNA(iter6).pos + nt_RNAPblocked 
                        PREmRNA(iter6).kr(1)=0;
                    end
                else
                    PREmRNA(iter6).kr=[];           
                end
%                 PREmRNA(iter6).kr= PREmRNA(iter6).kr./AV;
            end
        else
            PREmRNA(1).kr=[];
        end
        
        % PRM transcript ribos
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
                else
                    PRMmRNA(iter7).kr=[];           
                end
%                 PRMmRNA(iter7).kr= PRMmRNA(iter7).kr./AV;
            end
        else
            PRMmRNA(1).kr=[];
        end
        
        % PL transcript ribos
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
                else
                    PLmRNA(iter8).kr=[];
                end 
%             PLmRNA(iter8).kr= PLmRNA(iter8).kr./AV;
            end
        else
            PLmRNA(1).kr=[];
        end
                        
        tl_propen=[PRmRNA(:).kr PREmRNA(:).kr PRMmRNA(:).kr PLmRNA(:).kr];
%         tl_propen(tl_propen>0)
        propen(24+numel(tx_propen):23+numel(tx_propen)+numel(tl_propen)) = tl_propen;
        %% Propensity of Protein reactions
        % See table for details on these reactions
        % cI
        propen(1)=k(1)*n.cI;
        propen(2)=k(2)*(n.cI^2-n.cI)/2;
        propen(3)=k(3)*n.cI_2;
        % cro
        propen(4)=k(4)*n.cro;
        propen(5)=k(5)*(n.cro^2-n.cro)/2;    
        propen(6)=k(6)*n.cro_2;
        % N deg
        propen(7)=k(7)*n.N;
        % P1 deg pathway
        propen(8)=k(8)*n.cII*n.P1;
        propen(9)=k(9)*n.P1cII;
        propen(10)=k(10)*n.P1cII;
        propen(11)=k(11)*n.cIII*n.P1;
        propen(12)=k(12)*n.P1cIII;
        propen(13)=k(13)*n.P1cIII;
        % P2 deg pathway
        propen(14)=k(8)*n.cII*n.P2;
        propen(15)=k(9)*n.P2cII;
        propen(16)=k(10)*n.P2cII;
        propen(17)=k(11)*n.cIII*n.P2;
        propen(18)=k(12)*n.P2cIII;
        propen(19)=k(13)*n.P2cIII;  

        %% Choose reaction 
        rannum=rand(1);        
        rxn_ind=1: 23+numel(tx_propen)+numel(tl_propen);
        % Debug#
        if numel(propen) > numel(rxn_ind)
            disp('')
        end
        v_rxn_ind=rxn_ind(propen>0);
        v_propen= propen(propen>0);
        v_ind=find(cumsum(v_propen)./ sum(v_propen) > rannum,1);
        selected_rxn = v_rxn_ind(v_ind);
        
        % choose timestep
        delta_t = -log(rand(1)) /sum(propen);
        
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
                m_ind=0;
                switch find(tx_ind>=selected_rxn-23,1)
                    case 1 %PR
                        m_ind=selected_rxn-23;
                        m_ind=ceil(m_ind/7);
                        switch mod(selected_rxn-23,7)    
                            case 1
                                % km(1)= DNA+1                        
                                %                         disp('TX elongation')                                                                
                                PRmRNA(m_ind).pos = PRmRNA(m_ind).pos+1;
                                if PRmRNA(m_ind).pos >= PR_trans(2) % Check for end of transcript
                                    PRmRNA(m_ind).T = 1; 
                                end
%                                 % DEBUG
%                                 if m_ind == 1
%                                     disp(['PRmRNA(', num2str(m_ind), ').pos=', num2str(PRmRNA(m_ind).pos)])
%                                 end
                            case 2
                                % km(2)= ribo binds cro
                            %                         disp('ribosome binding')
                                PRmRNA(m_ind).ribopos = [PRmRNA(m_ind).ribopos cro_span(1)]; % approximating that RBS is located at start codon
%                                 % Debug
%                                 disp('Ribo bound cro!')
                            case 3
                                % km(3)= RNAse binds cro
                                % Debug
                            %                         disp('RNAse binding')
                                PRmRNA(m_ind).Dcro= 1; 
                                % if there are no ribosomes delete mRNA
                                if isempty(PRmRNA(m_ind).ribopos) && PRmRNA(m_ind).DcII
%                                     % DEBUG
%                                     PRmRNA
%                                     if m_ind ==1
%                                         disp(['deleting PRmRNA(1) at time=', num2str(t)])
%                                     else
%                                         disp(['deleting PRmRNA(', num2str(m_ind), ') at time=', num2str(t)])
%                                     end
                                    PRmRNA(m_ind)=[];
                                    if isempty(PRmRNA) % if the last mRNA is destroyed, create a blank mRNA with no pos
                                        PRmRNA(1)=create_mRNA(1);
                                    end                                    
                                    propen((0:6)+(7*(m_ind-1)+24))=[];
                                end
%                                 % Debug
%                                 if m_ind == 1
%                                     disp('RNAse bound cro!')
%                                 end
                            case 4
                                % km(4) = terminate TX
                                % Debug
                            %                         disp('TX termination')
                                PRmRNA(m_ind).T= 1;
                                % move the RNAP out of the way. 
                                % In reality it would unbind but this is a
                                % hack to keep the mRNA position value
                                PRmRNA(m_ind).pos=48000;                                
                                % DEG#
                                % signal cII degradation so terminated
                                % transcripts can still be degraded
                                % also prevents cII production despite
                                % above hack
                                PRmRNA(m_ind).DcII=1; 
                                if isempty(PRmRNA(m_ind).ribopos) && PRmRNA(m_ind).Dcro
%                                     % DEBUG
%                                     PRmRNA;
%                                     disp(['deleting PRmRNA(', num2str(m_ind), ') at time=', num2str(t)])
                                    PRmRNA(m_ind)=[];                                                                        
                                    if isempty(PRmRNA) % if the last mRNA is destroyed, create a blank mRNA with no pos
                                        PRmRNA(1)=create_mRNA(1);
                                    end
                                    propen((0:6)+(7*(m_ind-1)+24))=[];
                                end
                            case 5
                                % km(5)= toggle N state 
                                % Debug
                            %                         disp('toggle N state')
                                PRmRNA(m_ind).N= ~PRmRNA(m_ind).N;
                            case 6 
                                % km(6)= ribo binds cII
                                PRmRNA(m_ind).ribopos = [PRmRNA(m_ind).ribopos cII_span(1)]; % approximating that RBS is located at start codon                                
                            case 0
                                % km(7)= RNAse binds cII
                                PRmRNA(m_ind).DcII= 1;
                                % if there are no ribosomes delete mRNA
                                if isempty(PRmRNA(m_ind).ribopos) && PRmRNA(m_ind).Dcro
%                                     % DEBUG
%                                     PRmRNA
%                                     disp(['deleting PRmRNA(', num2str(m_ind), ') at time=', num2str(t)])                                    
                                    PRmRNA(m_ind)=[];
                                    if isempty(PRmRNA) % if the last mRNA is destroyed, create a blank mRNA with no pos
                                        PRmRNA(1)=create_mRNA(1);
                                    end
                                    propen((0:6)+(7*(m_ind-1)+24))=[];
                                end                                
                            otherwise
                                disp('Error, reaction not found')
                        end                        
                    case 2 %PRE
                        m_ind=selected_rxn-23-tx_ind(1);
                        m_ind=ceil(m_ind/7);
                        switch mod(selected_rxn-23,7)    
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
                                PREmRNA(m_ind).ribopos = [PREmRNA(m_ind).ribopos cI_span(1)]; % approximating that RBS is located at start codon
                            case 3
                                % km(3)= RNAse binds
                            %                         disp('RNAse binding')
                                PREmRNA(m_ind).D= 1;     
                                % if there are no ribosomes delete mRNA
                                if isempty(PREmRNA(m_ind).ribopos)
                                    PREmRNA(m_ind)=[];
                                    if isempty(PREmRNA)
                                        PREmRNA(1)=create_mRNA(2);
                                    end
                                propen((0:6)+(tx_ind(1)+7*(m_ind-1)+24))=[];
                                end
                            case 4
                                % km(4) = terminate TX
                                disp(['Error in step ' num2str(curstep), ', PRE cannot be terminated'])
                                PREmRNA(m_ind).T= 1;
                            case 5
                                % km(5)= toggle N state 
                            %                         disp('toggle N state')
                                disp(['Error in step ' num2str(curstep), ', N cannot bind PRE'])
                                PREmRNA(m_ind).N= ~PREmRNA(m_ind).N;
                            otherwise
                                disp(['Error in step ', num2str(curstep), ', reaction not found'])
                        end
                    case 3 %PRM
                        m_ind=selected_rxn-23-tx_ind(2);
                        m_ind=ceil(m_ind/7);
                        switch mod(selected_rxn-23,7)    
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
                                PRMmRNA(m_ind).ribopos = [PRMmRNA(m_ind).ribopos cI_span(1)]; % approximating that RBS is located at start codon
                            case 3
                                % km(3)= RNAse binds
                            %                         disp('RNAse binding')
                                PRMmRNA(m_ind).D= 1;
                                % if there are no ribosomes delete mRNA
                                if isempty(PRMmRNA(m_ind).ribopos)
                                    PRMmRNA(m_ind)=[];
                                    if isempty(PRMmRNA)
                                        PRMmRNA(1)=create_mRNA(3);
                                    end
                                propen((0:6)+(tx_ind(2)+7*(m_ind-1)+24))=[];
                                end
                            case 4
                                % km(4) = terminate TX
                            %                         disp('TX termination')
                                PRMmRNA(m_ind).T= 1;
                                disp(['Error in step ' num2str(curstep), ', PRM cannot be terminated'])
                            case 5
                                % km(5)= toggle N state 
                            %                         disp('toggle N state')
                                disp(['Error in step ' num2str(curstep), ', N cannot bind PRE'])                               
                                PRMmRNA(m_ind).N= ~PRMmRNA(m_ind).N;
                            otherwise
                                disp('Error, reaction not found')
                        end
                    case 4 %PL   
                        m_ind=selected_rxn-23-tx_ind(3);
                        m_ind=ceil(m_ind/7);
                        switch mod(selected_rxn-23,7)    
                            case 1
                                % km(1)= DNA+1                        
                                %                         disp('TX elongation')                            
                                PLmRNA(m_ind).pos = PLmRNA(m_ind).pos-1; % reverse
                                if PLmRNA(m_ind).pos <= PL_trans(2) % Check for end of transcript
                                    PLmRNA(m_ind).T = 1; 
                                    % move the RNAP out of the way. 
                                    % In reality it would unbind but this is a
                                    % hack to keep the mRNA position value
                                    PLmRNA(m_ind).pos=23000;
                                end
                            case 2
                                % km(2)= ribo binds N
                            %                         disp('ribosome binding')
                                PLmRNA(m_ind).ribopos = [PLmRNA(m_ind).ribopos N_span(1)]; % approximating that RBS is located at start codon
                            case 3
                                % km(3)= RNAse binds N RBS
                            %                         disp('RNAse binding')
                                PLmRNA(m_ind).DN= 1;    
                                % if there are no ribosomes delete mRNA
                                if isempty(PLmRNA(m_ind).ribopos) && PLmRNA(m_ind).DcIII; 
                                    PLmRNA(m_ind)=[];
                                    if isempty(PLmRNA)
                                        PLmRNA(1)=create_mRNA(4);
                                    end
                                    propen((0:6)+(tx_ind(3)+7*(m_ind-1)+24))=[];
                                end
                            case 4
                                % km(4) = terminate TX
                            %                         disp('TX termination')                                
                                PLmRNA(m_ind).T= 1;
                                % move the RNAP out of the way. 
                                % In reality it would unbind but this is a
                                % hack to keep the mRNA position value
                                PLmRNA(m_ind).pos=23000;
                                % DEG#
                                % Signal cIII degradation
                                PLmRNA(m_ind).DcIII=1;
                                % if there are no ribosomes delete mRNA
                                if isempty(PLmRNA(m_ind).ribopos) && PLmRNA(m_ind).DN
                                    PLmRNA(m_ind)=[];
                                    if isempty(PLmRNA)
                                        PLmRNA(1)=create_mRNA(4);
                                    end
                                    propen((0:6)+(tx_ind(3)+7*(m_ind-1)+24))=[];
                                end
                            case 5
                                % km(5)= toggle N state 
                            %                         disp('toggle N state')
                                PLmRNA(m_ind).N= ~PLmRNA(m_ind).N;
                            case 6
                                % km(6) = ribo binds cIII RBS
                                PLmRNA(m_ind).ribopos = [PLmRNA(m_ind).ribopos cIII_span(1)];
                            case 0
                                % km(7) = RNAse binds cIII RBS
                                PLmRNA(m_ind).DcIII= 1;
                                % if there are no ribosomes delete mRNA
                                if isempty(PLmRNA(m_ind).ribopos) && PLmRNA(m_ind).DN
                                    PLmRNA(m_ind)=[];
                                    if isempty(PLmRNA)
                                        PLmRNA(1)=create_mRNA(4);
                                    end
                                    propen((0:6)+(tx_ind(3)+7*(m_ind-1)+24))=[];
                                end
                            otherwise
                                disp('Error, reaction not found')
                        end
                    otherwise
                        disp('Error, reaction not found')
                end
                if isempty(PRmRNA)
                    disp(num2str(curstep))
                end
                                                    

            %% TL reactions
            case num2cell((24+numel(tx_propen)):(23+numel(tx_propen)+numel(tl_propen)))   
                tx_ind(1) = numel([PRmRNA(:).km]);
                tx_ind(2) = numel([PREmRNA(:).km]);
                tx_ind(3) = numel([PRMmRNA(:).km]);
                tx_ind(4) = numel([PLmRNA(:).km]);
                tx_ind = cumsum(tx_ind);  
                tl_ind(1) = numel([PRmRNA(:).kr]);
                tl_ind(2) = numel([PREmRNA(:).kr]);
                tl_ind(3) = numel([PRMmRNA(:).kr]);
                tl_ind(4) = numel([PLmRNA(:).kr]);
                tl_ind = cumsum(tl_ind);                
                switch find(tl_ind>=selected_rxn-23-numel(tx_propen),1)
                    case 1
                        % PR
                        m_ind=1;
                        while numel([PRmRNA(1:m_ind).kr]) < selected_rxn-23-numel(tx_propen)
                            m_ind=m_ind+1;
                        end
                        r_ind = selected_rxn-23-numel(tx_propen)-numel([PRmRNA(1:m_ind-1).kr]);
%                         % Debug
%                         if r_ind > numel(PRmRNA(m_ind).ribopos) 
%                             disp('') 
%                         end
                        % Advance the ribo one nt                        
                        PRmRNA(m_ind).ribopos(r_ind)=PRmRNA(m_ind).ribopos(r_ind)+1;
%                         disp([num2str(r_ind) 'at ' num2str(PRmRNA(m_ind).ribopos(r_ind))]) % Debug
                        if (PRmRNA(m_ind).ribopos(r_ind) > cro_span(2)) && (PRmRNA(m_ind).ribopos(r_ind) < cII_span(1))
                            % if the ribo is in between cro and cII then 
                            % the cro protein is finished
                            PRmRNA(m_ind).ribopos(r_ind) = []; % ribo unbinds
                            PRmRNA(m_ind).kr(r_ind)=[]; % clear propensity
                            propen(selected_rxn)=[]; % remove reaction
                            n.cro = n.cro + 1; % cro is produced
                            if PRmRNA(m_ind).Dcro && PRmRNA(m_ind).DcII && isempty(PRmRNA(m_ind).ribopos)
                                % if the mRNA is being degraded at both locations and there
                                % are no more active ribo
                                PRmRNA(m_ind)=[]; % mRNA is degraded
                                if isempty(PRmRNA) % if the last mRNA is destroyed, create a blank mRNA with no pos
                                    PRmRNA(1)=create_mRNA(1);
                                end
                                propen((0:6)+(7*(m_ind-1)+24))=[];
%                                 propen() = []; % delete propensities
                            end                        
                        elseif PRmRNA(m_ind).ribopos(r_ind) > cII_span(2)
                            % if the ribo is past cII, then the cII protein is
                            % finished
                            PRmRNA(m_ind).ribopos(r_ind) = []; % ribo unbinds
                            PRmRNA(m_ind).kr(r_ind)=[]; % clear propensity
                            propen(selected_rxn)=[];
                            n.cII = n.cII + 1; % cII is produced
                            if PRmRNA(m_ind).DcII && PRmRNA(m_ind).Dcro && isempty(PRmRNA(m_ind).ribopos)
                                % if the mRNA is being degraded at both locations and there
                                % are no more active ribo
                                PRmRNA(m_ind)=[]; % mRNA is degraded
                                if isempty(PRmRNA) % if the last mRNA is destroyed, create a blank mRNA with no pos
                                    PRmRNA(1)=create_mRNA(1);
                                end
                                propen((0:6)+(7*(m_ind-1)+24))=[];
                            end
                        end                            
                    case 2
                        % PRE
                        m_ind=1;                        
                        while numel([PREmRNA(1:m_ind).kr]) < selected_rxn-23-numel(tx_propen)-tl_ind(1)
                            m_ind=m_ind+1;
                        end
                        r_ind = selected_rxn-23-numel(tx_propen)-tl_ind(1)-numel([PREmRNA(1:m_ind-1).kr]);
                        % Advance the ribo one nt
                        PREmRNA(m_ind).ribopos(r_ind) = PREmRNA(m_ind).ribopos(r_ind)-1; % reverse
                        if PREmRNA(m_ind).ribopos(r_ind) < cI_span(2)
                            % if the ribo is "past" cI then
                            PREmRNA(m_ind).ribopos(r_ind)=[]; % ribo unbinds
                            PREmRNA(m_ind).kr(r_ind)=[]; % clear propensity
                            propen(selected_rxn)=[];
                            n.cI= n.cI + 1; % cI is produced
                            if PREmRNA(m_ind).D && isempty(PREmRNA(m_ind).ribopos)
                                % if the mRNA is being degraded and there 
                                % are no more active ribo
                                PREmRNA(m_ind)=[];
                                if isempty(PREmRNA)
                                    PREmRNA(1)=create_mRNA(2);
                                end
                                propen((0:6)+(tx_ind(1)+7*(m_ind-1)+24))=[];
                            end
                        end
                    case 3
                        % PRM
                        m_ind=1;                        
                        while numel([PRMmRNA(1:m_ind).kr]) < selected_rxn-23-numel(tx_propen)-tl_ind(2)
                            m_ind=m_ind+1;
                        end
                        r_ind = selected_rxn-23-numel(tx_propen)-tl_ind(2)-numel([PRMmRNA(1:m_ind-1).kr]);
                        % Advance the ribo one nt
                        PRMmRNA(m_ind).ribopos(r_ind) = PRMmRNA(m_ind).ribopos(r_ind)-1; % reverse
                        if PRMmRNA(m_ind).ribopos(r_ind) < cI_span(2)
                            % if the ribo is "past" cI then
                            PRMmRNA(m_ind).ribopos(r_ind)=[]; % ribo unbinds
                            PRMmRNA(m_ind).kr(r_ind)=[]; % clear propensity
                            propen(selected_rxn)=[];
                            n.cI= n.cI + 1; % cI is produced
                            if PRMmRNA(m_ind).D && isempty(PRMmRNA(m_ind).ribopos)
                                % if the mRNA is being degraded and there 
                                % are no more active ribo
                                PRMmRNA(m_ind)=[];
                                if isempty(PRMmRNA)
                                    PRMmRNA(1)=create_mRNA(3);
                                end
                                propen((0:6)+(tx_ind(2)+7*(m_ind-1)+24))=[];
                            end
                        end
                    case 4
                        % PL
                        m_ind=1;                        
                        while numel([PLmRNA(1:m_ind).kr]) < selected_rxn-23-numel(tx_propen)-tl_ind(3)
                            m_ind=m_ind+1;
                        end                        
                        r_ind = selected_rxn-23-numel(tx_propen)-tl_ind(3)-numel([PLmRNA(1:m_ind-1).kr]);
                        % FIXBUG run=LambdaModel(314159265)
                        % Debug
                        if r_ind > numel(PLmRNA(m_ind).ribopos)
                            disp('')
                        end
                        % Advance the ribo one nt
                        PLmRNA(m_ind).ribopos(r_ind) = PLmRNA(m_ind).ribopos(r_ind)-1;
                        if PLmRNA(m_ind).ribopos(r_ind) < N_span(2) && PLmRNA(m_ind).ribopos(r_ind) > cIII_span(1)
                            % if the ribo is left of N but right of cIII
                            PLmRNA(m_ind).ribopos(r_ind)=[]; % ribo unbinds
                            PLmRNA(m_ind).kr(r_ind)=[]; % clear propensity
                            propen(selected_rxn)=[];
                            n.N= n.N + 1; % N is produced
%                             % Debug
%                             try
%                                 PLmRNA(m_ind).DN && PLmRNA(m_ind).DcIII && isempty(PLmRNA(m_ind).ribopos)
%                             catch
%                                 disp('')
%                             end
                            if PLmRNA(m_ind).DN && PLmRNA(m_ind).DcIII && isempty(PLmRNA(m_ind).ribopos)
                                % if the mRNA is being degraded at both RBS and there 
                                % are no more active ribo
                                PLmRNA(m_ind)=[];
                                if isempty(PLmRNA) % if the last mRNA is destroyed, create a blank mRNA with no pos
                                    PLmRNA(1)=create_mRNA(4);
                                end
                                propen((0:6)+(tx_ind(3)+7*(m_ind-1)+24))=[];
                            end
                        elseif PLmRNA(m_ind).ribopos(r_ind) < cIII_span(2)
                            % if ribo is left of cIII
                            PLmRNA(m_ind).ribopos(r_ind)=[]; % ribo unbinds
                            PLmRNA(m_ind).kr(r_ind)=[]; % clear propensity
                            propen(selected_rxn)=[];
                            n.cIII= n.cIII + 1; % cIII is produced
                            if PLmRNA(m_ind).DN && PLmRNA(m_ind).DcIII && isempty(PLmRNA(m_ind).ribopos)
                                % if the mRNA is being degraded at both RBS and there 
                                % are no more active ribo
                                PLmRNA(m_ind)=[];
                                if isempty(PLmRNA) % if the last mRNA is destroyed, create a blank mRNA with no pos
                                    PLmRNA(1)=create_mRNA(4);
                                end
                                propen((0:6)+(tx_ind(3)+7*(m_ind-1)+24))=[];
                            end
                        end
                    otherwise
                        disp('Error, reaction not found')
                end
%                 tl_propen=[PRmRNA(:).kr PREmRNA(:).kr PRMmRNA(:).kr PLmRNA(:).kr];

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
                PRmRNA(ind1)=create_mRNA(1); % Create empty mRNA
                PRmRNA(ind1).pos = PR_trans(1); % start it at tss               
            case 21
                %PRE
                if isempty(PREmRNA(1).pos)
                    ind2=1;
                else
                    ind2 = numel(PREmRNA)+1;
                end
                PREmRNA(ind2)=create_mRNA(2);
                PREmRNA(ind2).pos = PRE_trans(1);
            case 22
                %PRM
                if isempty(PRMmRNA(1).pos)
                    ind3=1;
                else
                    ind3 = numel(PRMmRNA)+1;
                end
                PRMmRNA(ind3)=create_mRNA(3);
                PRMmRNA(ind3).pos = PRM_trans(1); 
            case 23
                %PL
                if isempty(PLmRNA(1).pos)
                    ind4=1;
                else
                    ind4 = numel(PLmRNA)+1;
                end
                PLmRNA(ind4)=create_mRNA(4);
                PLmRNA(ind4).pos = PL_trans(1);
            otherwise
                disp('Error, reaction not found')
        end

    % kPout(1) = kPR(state_PR);
    % kPout(2) = kPRE(state_PRE);
    % kPout(3) = kPRM(state_PR);
    % kPout(4) = kPL(state_PL);           


    %     tx_propen(selected_rxn-24)




    %% Advance next timestep

        

        t=t+delta_t;
        % update volume based on growth rate
        V=Vin+k0*t;
        
        %% Check end conditions
        if V > Vmax;
            f_continue=0;
        end

        if t > tmax
            f_continue=0;
        end
    end    
    
    %% Outputs
    % DEBUG
%     timeplot=timeout(1,1:curstep);
%     figure; plot(1:numel(timeplot),timeplot,'.')    
    
    curstep=max(curstep,curstep);

    run(runnum).t=timeout(runnum,1:curstep);
    run(runnum).V=Vout(runnum,1:curstep);
    
    run(runnum).cI=cIout(runnum,1:curstep);
    run(runnum).cro=croout(runnum,1:curstep);    
    run(runnum).cI_2=cI_2out(runnum,1:curstep);
    run(runnum).cro_2=cro_2out(runnum,1:curstep);
    run(runnum).N=Nout(runnum,1:curstep);
    run(runnum).cII=cIIout(runnum,1:curstep);
    run(runnum).cIII=cIIIout(runnum,1:curstep);
    
    run(runnum).PRmRNA=PRmRNAout(runnum,1:curstep);
    run(runnum).PREmRNA=PREmRNAout(runnum,1:curstep);
    run(runnum).PRMmRNA=PRMmRNAout(runnum,1:curstep);
    run(runnum).PLmRNA=PLmRNAout(runnum,1:curstep);
    run(runnum).randseed=varargin{1};
    
end

varargout{1}=run;

end

function mRNA = create_mRNA(promoter)
switch promoter
    case 1
        mRNA.pos=[];
        mRNA.ribopos=[];
        mRNA.N=0; % is N bound?
        mRNA.T=0; % is the mRNA terminated?
        mRNA.Dcro=0; % is the mRNA being degraded by RNAse at cro?
        mRNA.DcII=0;  % is the mRNA being degraded by RNAse at cII?
        mRNA.km=zeros(1,7);
        mRNA.kr=[];
    case 2
        mRNA.pos=[];
        mRNA.ribopos=[];
        mRNA.T=0;
        mRNA.D=0; % is the mRNA being degraded by RNAse?
        mRNA.km=zeros(1,7);
        mRNA.kr=[];
    case 3
        mRNA.pos=[];
        mRNA.ribopos=[];
        mRNA.T=0;
        mRNA.D=0; % is the mRNA being degraded by RNAse?
        mRNA.km=zeros(1,7);
        mRNA.kr=[];
    case 4
        mRNA.pos=[];
        mRNA.ribopos=[];
        mRNA.N=0; % is N bound?
        mRNA.T=0; % is the mRNA terminated?.
        mRNA.DN=0; % is the mRNA being degraded by RNAse at N?
        mRNA.DcIII=0; % is the mRNA being degraded by RNAse at cIII?
        mRNA.km=zeros(1,7);
        mRNA.kr=[];
    otherwise
        disp('error, promoter not found')
end
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





