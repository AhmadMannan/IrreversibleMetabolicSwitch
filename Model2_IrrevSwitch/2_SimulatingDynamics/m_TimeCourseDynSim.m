% ------------------------------------------------------------------------
% RUNNING TIME COURSE SIMULATIONS FOR A FIXED INDUCTION PROCESS
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 24/11/2020
% ------------------------------------------------------------------------

function [t,x] = m_TimeCourseDynSim(params,comp,hoursCult)


% --- [1] Finding steady state w/o induction = init cond -----------------

OA0 = 0; % absense of inducer
SS_0OA = m_DoseResp(params,OA0,comp,[]);
SS_0OA = abs(real(SS_0OA));

% There may be multiple stable steady states of FadR in the absence of
% inducer. Find stable steady state for when FadR level is highest:
stbl_p = find(SS_0OA(:,end) == 1); % find position of all stable states
SS_0OA = SS_0OA(stbl_p,2:end-1); %#ok<*FNDSB>
p_sshR = find(SS_0OA(:,1) == max(SS_0OA(:,1))); % select stable state with highest FadR level
SS_0OA = SS_0OA(p_sshR,:);


% --- [2] Generate simulation time course --------------------------------

% Define simulation timespan and initial conditions:
timespan = linspace(0,hoursCult,20000);
% timespan = linspace(0,hoursCult,3000);
initCond = [SS_0OA,0];

% Simulate time course simulation:
[t,x] = ode15s(@(t,x)Model2_Dyn(t,x,params,comp),timespan,initCond);
% [t,x] = odeSun(@Model2_Dyn,timespan,initCond,params);
