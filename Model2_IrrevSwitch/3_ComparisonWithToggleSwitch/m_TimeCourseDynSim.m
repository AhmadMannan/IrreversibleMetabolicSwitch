% ------------------------------------------------------------------------
% RUNNING TIME COURSE SIMULATIONS FOR A FIXED INDUCTION PROCESS
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 22/03/2021
% ------------------------------------------------------------------------

function [t,x] = m_TimeCourseDynSim(params,hoursCult)


% --- [1] Finding steady state w/o induction = init cond -----------------

I0 = 0; % absense of inducer
[SS_0I,~] = m_DoseResp(params,I0);
SS_0I = abs(real(SS_0I));

% There may be multiple stable steady states of LacI in the absence of
% inducer, finding steady state when LacI conc is highest:
stbl_p = find(SS_0I(:,1+3) == max(SS_0I(:,1+3)));
SS_0I = SS_0I(stbl_p,2:end); %#ok<*FNDSB>


% --- [2] Generate simulation time course --------------------------------

% Define simulation timespan and initial conditions:
timespan = linspace(0,hoursCult,20000);
% timespan = linspace(0,hoursCult,3000);
initCond = SS_0I;

% Simulate time course simulation:
[t,x] = ode15s(@(t,x)Model3_TS_Dyn(t,x,params),timespan,initCond);
% [t,x] = odeSun(@Model2_Dyn,timespan,initCond,params);
