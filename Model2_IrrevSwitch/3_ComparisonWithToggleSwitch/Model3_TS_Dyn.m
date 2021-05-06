% ------------------------------------------------------------------------
% MODEL 3 - KINETIC MODEL OF BISTABLE GENETIC TOGGLE SWITCH
% - Gardner et al (2000)
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 11/03/2021
% ------------------------------------------------------------------------

function [dx,flag,new_data,Ifeed_t,mu] = Model3_TS_Dyn(t,xvals,params)


% --- [0] MODEL DETAILS --------------------------------------------------
%
% SPECIES
% ----------------------------------------
% METABOLISM AND METABOLIC CONTROL
% ----------------------------------------
%    1    | Ix = IPTG in external media (inducer)
%    2    | I  = IPTG internal
%    4    | C  = LacI.IPTG seq. complex
% ----------------------------------------
% TRANSCRIPTIONAL CONTROL AND EXPRESSION
% ----------------------------------------
%    3    | L  = Regulator LacI
%    5    | T  = TetR
%    6    | Eg = Enzyme governing growth
% ----------------------------------------
%
% UNITS
%  - Concentrations = micro-Molar (uM)
%  - Time = hours
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% [1] PRELIMINARY
% ------------------------------------------------------------------------
% Output vector:
dx = zeros(6,1);

% Defining model variables:
Ix = xvals(1); % external IPTG
I  = xvals(2); % internalized IPTG
L  = xvals(3); % LacI
C  = xvals(4); % LacI.IPTG Complex
T  = xvals(5); % TetR
Eg = xvals(6); % Eg - enzyme governing growth


% ------------------------------------------------------------------------
% [2] MODEL - PARAMETERS
% ------------------------------------------------------------------------
p = params;

% Growth rate:
  muMax = p(1);  % Factor converting number of enzymes to growth rate, assuming max growth rate (mu0) achieved for max enzyme (Eg = 1).  
  sT    = p(2);  % Minimum growth rate when no Ep expressed, based on Usui 2012 report that mu is reduced to 30% of its value when pgi is knocked out (Usui 2013), and growth in our cell is at 0.1818 1/h.

% Induction parameters:
  Ifeed = p(3);  % Feed-in rate of inducing oleic acid
  Texp = p(4);  % Time length of induction exposure

% Fwd and rev rates of passive diffusion of IPTG from out to in cell:
  ki    = p(5);  % rate of diffusion - set very high, assuming rapid diffusion, with uptake working close to chemical equilibrium

% Fwd and rev rates of IPTG sequestration of LacI:
  kf    = p(6);  % Forward rate of sequestration
  kr    = p(7);  % Reverse rate of sequestration

% Expression and regulation of LacI:
  b_L   = p(8);  % Basal leaky expression of LacI synthesis.
  a_L   = p(9);  % Maximum rate of LacI synthesis.
  K_L   = p(10); % Affinity of TetR to inhibit LacI synthesis.

% Expression and regulation of TetR:
  b_T   = p(11); % Basal leaky expression of T synthesis.
  a_T   = p(12); % Maximum rate of T synthesis.
  K_T   = p(13); % Affinity of LacI to inhibit TetR synthesis.
  
% Eg expression and regulation:
  a_g   = p(14); % Maximum synthesis rate of enzyme Eg.
  K_g   = p(15); % Affinity of FadR to activate Eg synthesis.

% Indicator for when we consider dynamics with fixed inducer conc in media:
  conI  = p(16);
  
% Number of IPTG binding to LacI:
  n = p(17);
  

% ------------------------------------------------------------------------
% [3] MODEL - SYSTEM OF ODES
% ------------------------------------------------------------------------

% Growth Rate - As a function of Eg:
mu = muMax * ((Eg*sT) - sT + 1);

% Ensuring concs do not become negative:
  if Ix <= 0
      Ix = 0;
  end
  if I <= 0
      I = 0;
  end
  if C <= 0
      C = 0;
  end
  
% Define ODE for external IPTG (Ix):
  uptake = ki * (Ix - I);
  if conI == 1
      dx(1) = 0;
  elseif conI ~= 1
      % Define flux of IPTG feed into media:
      if t <= Texp
          Ifeed_t = Ifeed; % parameter: feed-in rate
      else
          Ifeed_t = 0;
      end
      dx(1) = Ifeed_t - uptake;
  end
  
% IPTG in cell (I):
%   n = 2;
  sequestration = kf*(I^n)*L - kr*C;
  dx(2) = uptake - (n*sequestration) - (mu*I);
  
% LacI expression and sequestration (dimer - L):
  dx(3) = b_L + a_L./(1 + (K_L*T).^2) - sequestration - (mu*L);
  
% LacI-IPTG-IPTG complex (C):
  dx(4) = sequestration - (mu*C);
  
% TetR expression (T):
  dx(5) = b_T + a_T/(1 + (K_T*L).^1) - (mu*T);
  
% Eg synthesis (Eg):
  dx(6) = (a_g ./ (1 + (K_g*T).^2)) - (mu*Eg);

% Outputs required for Sundials solver:
  flag = 0;
  new_data = [];
