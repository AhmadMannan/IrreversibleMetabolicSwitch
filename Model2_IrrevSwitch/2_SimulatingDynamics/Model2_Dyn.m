% ------------------------------------------------------------------------
%             MODEL 2 -  KINETIC MODEL OF IRREVERSIBLE SWITCH
%                       Model of irreversible switch
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 24/11/2020
% ------------------------------------------------------------------------


function [dx, OAin_t] = Model2_Dyn(t,xvals,params,comp)


% --- [0] MODEL DETAILS --------------------------------------------------
%
% SPECIES
% ----------------------------------------
% TRANSCRIPTIONAL CONTROL AND EXPRESSION
% ----------------------------------------
%    1    | R  = Regulator FadR
%    2    | D  = Uptake enzyme FadD
%         | B  = PlsB - Assumed at quasi-steady state (set as param)
%    5    | Eg = Enzyme governing growth
%    6    | T  = TetR
% ----------------------------------------
% METABOLISM AND METABOLIC CONTROL
% ----------------------------------------
%    3    | A  = Acyl-CoA
%    4    | C  = FadR.Acyl-CoA seq. complex
%    7    | OA = Oleic acid (inducer)
% ----------------------------------------
%
% UNITS
%  - Concentrations = micro-Molar (uM)
%  - Time = Hours
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% [1] PRELIMINARY
% ------------------------------------------------------------------------
% Output vector:
dx = zeros(7,1);

% Defining model variables:
R  = xvals(1); % FadR
D  = xvals(2); % FadD
A  = xvals(3); % Acyl-CoA
C  = xvals(4); % FadR.Acyl-CoA Complex
Eg = xvals(5); % Eg - enzyme governing growth
T  = xvals(6); % TetR
OA = xvals(7); % Oleic acid 


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

% TetR expression reg. by FadR, and FadR expression reg. by TetR:
  b_T   = p(18); % Basal leaky expression of T synthesis.
  a_T   = p(19); % Maximum rate of T synthesis.
  K_Ri  = p(20); % Affinity of FadR to inhibit TetR synthesis.
  K_T   = p(21); % Affinity of TetR to inhibit FadR synthesis.

% Affinity of FadR to inhibit production enzymes:
  K_p   = p(22);
  
% Induction parameters:
  OAin  = p(23); % Feed-in rate of inducing oleic acid
  OAxt  = p(24); % Time length of induction exposure


% ------------------------------------------------------------------------
% [3] MODEL - SYSTEM OF ODES
% ------------------------------------------------------------------------

% Growth Rate - As a function of Eg:
mu = muMax * ((Eg*sT) - sT + 1);


% --- TRANSCRIPTIONAL CONTROL --------------------------------------------

% FadR expression:
  if strcmp(comp,'COMP') == 1        % FadR-TetR competitive binding
      P_R = a_R * K_R * R / (1 + K_R*R + (K_T*T).^2);
  elseif strcmp(comp,'NONCOMP') == 1 % FadR-TetR non-competitive binding
      P_R = a_R * K_R * R / ((1 + K_R*R)*(1 + (K_T*T).^2));
  elseif strcmp(comp,'NOPAR') == 1
      P_R = a_R ./ (1 + (K_T*T).^2);
  end
  dx(1) = b_R + P_R - (kf*R*(A^2) - kr*C) - (mu * R);

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


% --- EXTRA CONTROL ON FADR EXPRESSION BY TETR ---------------------------

% TetR synthesis:
  dx(6) = b_T + a_T/(1 + (K_Ri*R).^2) - (mu*T);

  
% --- OLEIC ACID CONSUMPTION ---------------------------------------------

% Oleic acid:
  if t <= OAxt
      OAin_t = OAin; % parameter: feed-in rate
  else
      OAin_t = 0;
  end
  dx(7) = OAin_t - uptake;
