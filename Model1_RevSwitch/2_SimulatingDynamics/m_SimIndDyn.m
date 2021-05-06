% ------------------------------------------------------------------------
%   SIMULATE DYNAMICS OF MODEL WITH FADR REGULATING SWITCH TO PRODUCTION
%                  Under dynamic availability of inducer
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 12/11/2020
% ------------------------------------------------------------------------

function [t,x] = m_SimIndDyn(arch,fig,params,hoursCult)


%%
% Define model parameters:
% arch = 'PAR';
% hoursCult = 100;
% p = Model1_Params(arch);
p = params;


% --- [1] Calculate system steady state without inducer ------------------

% SS_0OA = Model_FAU_DetSteadyStates(p,arch);
SS0 = m_DoseResp(p,arch,0,[],[]);
SS0 = real(round(10000*SS0(2:end-1))/10000);


% --- [2] Simulate induction dynamics ------------------------------------

% Define simulation timespan and initial conditions:
timespan = [0,hoursCult];
OA0 = 0;
initCond = [SS0,OA0];

% Simulate time course simulation:
[t,x] = ode15s(@(t,x)Model1_Dyn(t,x,p,arch),timespan,initCond);

% ### ADDITIONAL FOR PLOTTING ###
% Show steady state for 10 hours before induction:
% t = [-10; 0; t];
% x = [[SS_0OA(1:end-1),0];[SS_0OA(1:end-1),0]; x];


% --- [3] Plotting time courses ------------------------------------------

if isempty(fig) == 0
    figure(fig); %clf

    if strcmp(arch,'NAR') == 1
        colr = 'r';
    elseif strcmp(arch,'PAR') == 1
        colr = 'b';
    end
    
    subplot(1,2,1)
    % Plot of OA dynamics vs time:
    plot(t,x(:,6),'-','LineWidth',1,'Color',[0.3,0.3,0.3])
    % Plot of free FadR vs time:
    hold on
    plot(t,x(:,1),'-','LineWidth',2,'Color',colr)
    hold off
    xlabel('Time (h)')
    ylabel('Concentration (\mu M)')
    xlim([min(t),max(t)])
    ylim([0,0.1])
    hold on
    plot([-10,100],[1/p(18),1/p(18)],'k:')
    hold off
    legend('OA','FadR','1/K_p','Location','NorthEast')
    
    subplot(1,2,2)
    % Plot of acyl-CoA vs time:
    plot(t,x(:,3),'k-','LineWidth',1)
    ylabel('Acyl-CoA (\mu M)')
    xlabel('Time (h)')
end

