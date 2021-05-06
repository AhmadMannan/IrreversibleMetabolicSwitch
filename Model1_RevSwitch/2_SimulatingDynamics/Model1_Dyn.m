% ------------------------------------------------------------------------
%     KINETIC MODEL - FA REG APPLIED TO CONTROL GROWTH AND PRODUCTION
%                   Modelling dynamic feed-in of inducer
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 12/11/2020
% ------------------------------------------------------------------------


function [dx,OAin] = Model1_Dyn(t,xvals,params,arch)


% --- [0] MODEL DETAILS --------------------------------------------------
%
% SPECIES
% -----------------------------------------
% TRANSCRIPTIONAL CONTROL AND EXPRESSION
% -----------------------------------------
%    1    | R  = Regulator FadR
%    2    | D  = Uptake enzyme FadD
%         | B  = PlsB - Assumed at quasi-steady state (set as param)
%    5    | Eg = Enzyme governing growth
% -----------------------------------------
% METABOLISM AND METABOLIC CONTROL
% -----------------------------------------
%    3    | A  = Acyl-CoA
%    4    | C  = FadR.Acyl-CoA seq. complex
%    6    | OA = Oleic acid availability
% -----------------------------------------
%
% UNITS
%  - Concentrations = micro-Molar (uM)
%  - Time = Hours
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% [1] PRELIMINARY
% ------------------------------------------------------------------------
% Output vector:
dx = zeros(6,1);

% Defining model variables:
R  = xvals(1); % FadR
D  = xvals(2); % FadD
A  = xvals(3); % Acyl-CoA
C  = xvals(4); % FadR.Acyl-CoA Complex
Eg = xvals(5); % Eg - enzyme governing growth
OA = xvals(6); % Oleic acid


% ------------------------------------------------------------------------
% [2] MODEL - PARAMETERS
% ------------------------------------------------------------------------
p = params;

% FadR expression and regulation:
  b_R   = p(1);  % Basal leaky expression of FadR.
  a_R   = p(2);  % Maximum rate of expression of FadR.
  K_R   = p(3);  % Affinity (ratio of rate of binding to unbinding) of FadR to its promoter.

% FadD expression and regulation:
  b_D   = p(4);  % Basal leaky expression of FadD.
  a_D   = p(5);  % Maximum rate of expression of FadD.
  K_D   = p(6);  % Affinity (ratio of rate of binding to unbinding) of FadR to FadD promoter.

% Uptake Reaction (single pseudo-rxn representing FadL and FadD):
  kcatD = p(7); % Turnover rate of enzyme of FadD.
  Km_D  = p(8); % Michaelis-constant of enzyme of FadD for Fatty Acid (OA).
  
% PlsB enzyme kinetics - Glycerol-3-Phosphate Acyltransferase:
  kcatB = p(9); % Turnover rate of enzyme of PlsB.
  Km_B  = p(10); % Michaelis-constant of enzyme of PlsB for Acyl-CoA.
  B     = p(11); % Assumed constant levels of PlsB enzyme.

% Sequestration of FadR by Acyl-CoA:
  kf    = p(12); % Effective forward rate of FadR sequestration by Acyl-CoA.
  kr    = p(13); % Effective reverse rate of sequestering.

% Eg expression and regulation:
  a_g   = p(14); % Maximum synthesis rate of enzyme Eg.
  K_g   = p(15); % Affinity of FadR to activate Eg synthesis.
  
% Growth rate:
  muMax = p(16); % Factor converting number of enzymes to growth rate, assuming max growth rate (mu0) achieved for max enzyme (Eg = 1).  
  sT    = p(17); % Minimum growth rate when no Ep expressed, based on Usui 2012 report that mu is reduced to 30% of its value when pgi is knocked out (Usui 2013), and growth in our cell is at 0.1818 1/h.

% Affinity of FadR to inhibit Ep synthesis.
  K_p   = p(18);

% Feed-in rate of inducing OA:
  OAin  = p(19);
  SigOS = p(20); % Signal overshoot (i.e. space between low FadR and even lower FadR states).


% ------------------------------------------------------------------------
% [3] MODEL - SYSTEM OF ODES
% ------------------------------------------------------------------------

% Growth Rate - As a function of Eg:
mu = muMax * ((Eg*sT) - sT + 1);


% --- TRANSCRIPTIONAL CONTROL --------------------------------------------

% FadR expression:
  if strcmp(arch,'NAR') == 1     % Negative autoregulation
      P_R = a_R / (1 + (K_R*R));
  elseif strcmp(arch,'PAR') == 1 % Positive autoregulation
      P_R = a_R * (K_R*R) / (1 + (K_R*R));
  end
  dx(1) = b_R + P_R - (kf*R*(A^2) - kr*C) - (mu*R);

% FadD expression:
  dx(2) = b_D + (a_D/(1 + ((K_D*R)^2))) - (mu*D);


% --- METABOLISM ---------------------------------------------------------

% Oleic Acid (OA) concentration in media:
% OA = p(26);
  if OA <= 0
      OA = 0;
  end
  if A <= 0
      A = 0;
  end
  if C <= 0
      C = 0;
  end
  
% Acyl-CoA concentration:
  uptake = kcatD * OA * D / (Km_D + OA);
  PlsBdrain = kcatB * B * A / (Km_B + A);
  dx(3) = uptake - PlsBdrain - (2*(kf*R*(A^2) - kr*C)) - (mu*A);

% FadR.Acyl-CoA complex concentration:
  dx(4) = (kf*R*(A^2) - kr*C) - (mu*C);


% --- CONTROL ON GROWTH --------------------------------------------------

% Eg synthesis:
  dx(5) = (a_g*K_g*R / (1 + K_g*R)) - (mu*Eg);


% --- OLEIC ACID INDUCTIONS AND CONSUMPTION ------------------------------

% Oleic acid:
  if R > 1/K_p
      OAin_t = OAin; % parameter: feed-in rate
  elseif R <= SigOS/K_p
      OAin_t = 0;
  elseif R <= 1/K_p && R > SigOS/K_p
      if dx(1) < 0
          OAin_t = OAin;
      else
          OAin_t = 0;
      end
  end
  dx(6) = OAin_t - uptake;
