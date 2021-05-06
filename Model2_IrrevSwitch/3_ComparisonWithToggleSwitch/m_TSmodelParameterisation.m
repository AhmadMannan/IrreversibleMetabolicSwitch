% ------------------------------------------------------------------------
% FITTING MODEL TO DATA
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 23/03/2021
% ------------------------------------------------------------------------

function [dataDR,dataK,DR,K] = m_TSmodelParameterisation(p_2opt)


%%
% --- Load and plot data for fitting params ------------------------------

% Data of the dose-response of IPTG vs monomers of LacI:
dataDR = xlsread('Data_LacIvsIPTG_Xu2009.xlsx','B2:C46');
figure(1); clf
hold on
plot(dataDR(:,1),dataDR(:,2),'x','Color',[0.6,0.6,0.6])
hold off
set(gca,'XScale','log')
xlabel('IPTG (\mu M)')
ylabel('Relative GFP prop. to LacI (a.u.)')

% Data of the kinetics of IPTG binding to monomers of LacI:
dataK = xlsread('Data_LacIIPTGbindkinetics_Xu2009.xlsx','A1:B80');
figure(2); clf
hold on
plot(dataK(:,1)/3600,dataK(:,2),'x','Color',[0.6,0.6,0.6])
hold off
xlabel('Time (h)')
ylabel('Relative GFP prop. to LacI (a.u.)')


% --- Simulate model of LacI-IPTG binding kinetics and dose-response -----

% Defining parameters:
params = Model3_TS_Dyn_Params;
params(17) = 1; % specifying 1 IPTG binds to 1 lacI monomer
% ################################
  params([6,7]) = [p_2opt(1),p_2opt(2)];
%   params([6,7]) = [750,900]*1.3;
% ################################

% Generating dose-response:
% rngInd = logspace(-2,2,50);
rngInd = dataDR(:,1);
dr = [];
for i = 1:length(rngInd)
    ic = [rngInd(i),0.15,0]; % in conc units of micoM
    tspan = [0,10]; % in time units of hours
    [t,x] = ode15s(@(t,x)model_LacIIPTGparams(t,x,params),tspan,ic);
    
    % check if steady state reached:
    dx = model_LacIIPTGparams(t(end),x(end,:),params);
    if max(abs(dx)) > 1e-4
        disp('Steady state not reached!')
    end
    
    % saving data for dose-response:
    dr = [dr; rngInd(i), x(end,2)]; %#ok<*AGROW>
end
DR = [dr(:,1),dr(:,2)./max(dr(:,2))];

% Plotting dose-response:
figure(1)
hold on
plot(DR(:,1),DR(:,2),'k-','LineWidth',2)
hold off


% Simulating kinetic response of LacI-IPTG binding:
ic2 = [50,1,0]; % in conc units of microM
tspan2 = dataK(:,1)/3600;
[t2,x2] = ode15s(@(t,x)model_LacIIPTGparams(t,x,params),tspan2,ic2);
K = [t2,x2(:,2)./max(x2(:,2))];

% Plotting time course:
figure(2)
hold on
plot(K(:,1),K(:,2),'k-','LineWidth',2)
hold off
