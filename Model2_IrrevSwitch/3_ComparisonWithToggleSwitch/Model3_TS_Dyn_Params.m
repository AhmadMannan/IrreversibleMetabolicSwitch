% ------------------------------------------------------------------------
% MODEL 3 PARAMETERS
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 12/03/2021
% ------------------------------------------------------------------------

function params = Model3_TS_Dyn_Params


% Growth rate:
  muMax = 0.1818;  % Factor converting number of enzymes to growth rate, assuming max growth rate (mu0) achieved for max enzyme (Eg = 1)
  sT    = 0.7;  % Minimum growth rate when no Ep expressed, based on Usui 2012 report that mu is reduced to 30% of its value when pgi is knocked out (Usui 2013), and growth in our cell is at 0.1818 1/h

% Induction parameters:
  Ifeed = 0;  % Feed-in rate of inducing oleic acid
  Texp  = 0;  % Time length of induction exposure

% Fwd and rev rates of passive diffusion of IPTG from out to in cell:
  ki    = 3600;  % (1/h) rate of diffusion - set very high, assuming rapid diffusion, with uptake working close to chemical equilibrium

% Fwd and rev rates of IPTG sequestration of LacI:
  load('kOpts')
  kf    = kOpts(1); % 800 (1/h)
  kr    = kOpts(2); % 1000 (1/h)

% Expression and regulation of LacI:
  a_L   = 0.1409*3.1623*1.2; %156.25*muMax;  % based on Gardner 2000  
  b_L   = 0.0005;  % arbitrary
  K_L   = 100; % affinity of TetR to inhibit LacI synthesis

% Expression and regulation of TetR:
  a_T   = 0.0517*5.0119; %15.6*muMax; % maximum synthesis rate
  b_T   = 0.1*0.0108; % leaky expression
  K_T   = 305.95; % affinity of LacI to inhibit TetR synthesis
  
% Eg expression and regulation:
  a_g   = muMax; % Maximum synthesis rate of enzyme Eg
  K_g   = 100; % affinity of TetR to inhibit Eg synthesis
  
  
% Indicator of whether we are considering dynamics with constant inducer
% presence:
  constInd = 0;
  
  
% Number of IPTG molecules binding LacI:
  n = 2;
  
  
params = [muMax; sT; Ifeed; Texp; ki; kf; kr; b_L; a_L; K_L; b_T; a_T; K_T; a_g; K_g; constInd; n];
