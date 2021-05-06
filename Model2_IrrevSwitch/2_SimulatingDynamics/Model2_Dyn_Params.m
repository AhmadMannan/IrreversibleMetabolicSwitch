% ------------------------------------------------------------------------
% MODEL 2 PARAMETERS FOR SIMULATING INDUCTION REGIME
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 24/11/2020
% ------------------------------------------------------------------------

function params = Model2_Dyn_Params

% Load parameters found from previous analysis that achieves an
% irreversible switch:
load('Params_Eng_IrrevSw')
p0 = Params_Eng_IrrevSw;

% Define 1/Kp as FadR conc when system is said to be in induced state:
K_p = p0(6); % same as affinity of FadR to bind to operator on fadD prom., assuming we use fadD promoter to control expression of synthesis genes

% Define induction process parameters:
OAin  = 0; % Feed-in rate of inducing oleic acid
OAxt  = 0; % Time length of induction exposure

params = [p0; K_p; OAin; OAxt];