% ------------------------------------------------------------------------
% Solve ODE model with Sundials solver CVode
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% ------------------------------------------------------------------------

function [t,y] = odeSun(model,timePts,initConds,params)


% Define parameters as model input:
options = CVodeSetOptions('UserData',params,'RelTol',1e-4,'AbsTol',1e-6);

% Define time points over which to integrate:
tStart = timePts(1);
tEnd = timePts(end);
numtimePts = length(timePts);
dt = (tEnd-tStart)/(numtimePts-1);
tt = linspace(tStart+dt,tEnd,numtimePts-1);

% Define initial conditions:
s = size(initConds);
if s(1) == 1
    initConds = initConds';
end

% Define solver method and integrate:
CVodeInit(model, 'BDF', 'Newton', tStart, initConds, options); 
[~,t,y] = CVode(tt,'Normal');

% Add initial point to full time course:
t = [tStart,t]';
y = [initConds, y]';
