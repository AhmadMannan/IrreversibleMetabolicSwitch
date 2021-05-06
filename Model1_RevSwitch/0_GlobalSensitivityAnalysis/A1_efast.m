function s_KM = A1_efast(arch,metric,fig)

% First order and total effect indices for a given
% model computed with Extended Fourier Amplitude
% Sensitivity Test (EFAST).
% Andrea Saltelli, Stefano Tarantola and Karen Chan.
% 1999. % "A quantitative model-independent method for global
% sensitivity analysis of model output". % Technometrics 41:39-56
% clear all;
% clc
% close all;
%% INPUT
NR = 3; %: 6; no. of search curves - RESAMPLING
k = 7 + 1; % # of input factors (parameters varied) + dummy parameter
NS = 100; % 257; # of samples per search curve
wantedN = NS*k*NR; % wanted no. of sample points
VarOI = 1; % variable of interest, on which eFAST sensitivity index is based on
% f = @ODE_LotkaVolterra; % Model file name, with function handle

% OUTPUT
% SI[] : first order sensitivity indices
% STI[] : total effect sensitivity indices
% Other used variables/constants:
% OM[] : vector of k frequencies
% OMi : frequency for the group of interest
% OMCI[] : set of freq. used for the compl. group
% X[] : parameter combination rank matrix
% AC[],BC[]: fourier coefficients
% FI[] : random phase shift
% V : total output variance (for each curve)
% VI : partial var. of par. i (for each curve)
% VCI : part. var. of the compl. set of par...
% AV : total variance in the time domain
% AVI : partial variance of par. i
% AVCI : part. var. of the compl. set of par.
% Y[] : model output

MI = 4; %: maximum number of fourier coefficients
% that may be retained in calculating the partial
% variances without interferences between the
% assigned frequencies

%% PARAMETERS AND ODE SETTINGS (they are included in the following file)
[p0,p_var_pos,efast_var,pmin,pmax,time_points,y0,y_var_label] = m_efast_ParamSettings(arch,metric);

% Computation of the frequency for the group
% of interest OMi and the # of sample points NS (here N = NS)
OMi = floor(((wantedN/NR)-1)/(2*MI)/k);
NS = 2*MI*OMi+1;
if(NS*NR < 65)
    fprintf(['Error: sample size must be >= ' ...
    '65 per factor.\n']);
    return;
end


%% Pre-allocation of the output matrix Y
%% Y will save only the points of interest specified in
%% the vector time_points
% Y = [];
Y(NS,length(time_points),length(y0),length(pmin),NR)=0;  % pre-allocation

% Loop over k parameters (input factors)

tic
% OM = [];
% Y = [];
for i = 1:k % i = # of replications (or blocks)
    
    % Algorithm for selecting the set of frequencies.
    % OMci(i), i = 1:k-1, contains the set of frequencies
    % to be used by the complementary group.
    OMci = SETFREQ(k,OMi/2/MI,i);   
    
    % Loop over the NR search curves.
    for L = 1:NR
        
        % Setting the vector of frequencies OM for the k parameters
        cj = 1;
        for j=1:k
            if(j==i)
                % For the parameter (factor) of interest
                OM(i) = OMi;
            else
                % For the complementary group.
                OM(j) = OMci(cj);
                cj = cj+1;
            end
        end
        
        % Setting the relation between the scalar
        % variable S and the coordinates
        % {X(1),X(2),...X(k)} of each sample point.
        FI = rand(1,k)*2*pi; % random phase shift
        S_VEC = pi*(2*(1:NS)-NS-1)/NS;
        OM_VEC = OM(1:k);
        FI_MAT = FI(ones(NS,1),1:k)';
        ANGLE = OM_VEC'*S_VEC+FI_MAT;
        
        X(:,:,i,L) = 0.5+asin(sin(ANGLE'))/pi; % between 0 and 1
        
        % Transform distributions from standard
        % uniform to general.
        X(:,:,i,L) = parameterdist(X(:,:,i,L),pmax,pmin,0,1,NS,'unif'); %%this is what assigns 'our' values rather than 0:1 dist
        
        % Do NS model evaluations.
        parfor run_num = 1:NS
            [i run_num L] % keeps track of [parameter run NR]
            % ODE system file
            % f = @ODE_LotkaVolterra;
            % ODE solver call
            % [t,y] = ode15s(@(t,y)f(t,y,X(:,:,i,L),run_num),tspan,y0,[]);
            
            % Redefining parameters as needed in model:
            p_t = X(:,:,i,L);
            p_t = p_t(run_num,:);
            params_t = p0;
            params_t(p_var_pos) = p_t(1:end);
            
            % Generating dose-response curve and measuring BR, IT, RT:
            [BR,IT,RT] = m_QuantDoseResp_BR_IT_RT(params_t,arch,[]);
            
            % Defining which metric we are measuring sensitivity of:
            res = [];
            if metric == 'BR' %#ok<*BDSCA>
                res = BR;
            elseif metric == 'IT'
                res = IT;
            elseif metric == 'RT'
                res = RT;
            else
                error('State the dose-response feature we are studying the sensitivity of!')
            end
            
            % It saves only the output at the time points of interest
            % Y(run_num,:,:,i,L) = y(time_points+1,:);
            Y(run_num,:,:,i,L) = res(time_points,:);
        end %run_num=1:NS
    end % L=1:NR
end % i=1:k
% save Model_efast.mat;
toc

%
% CALCULATE Si AND STi for each resample (1,2,...,NR) [ranges]
[Si,Sti,rangeSi,rangeSti] = efast_sd(Y,OMi,MI,time_points,1:length(y0));

% Calculate Coeff. of Var. for Si and STi for Viral load (variable 4). See
% online Supplement A.5 for details.
[CVsi CVsti] = CVmethod(Si,rangeSi,Sti,rangeSti,1);

% T-test on Si and STi for Viral load (variable 4)
s_KM = efast_ttest(Si,rangeSi,Sti,rangeSti,1:length(time_points),efast_var,VarOI,y_var_label,0.05);


% PLOTTING RESULTS
figure(fig); clf
m_BarPlot_Sensitivities
title({['eFAST (GSA) to find sensitivity of: ',char(y_var_label)]; ['Circuit: ', char(arch)]})

