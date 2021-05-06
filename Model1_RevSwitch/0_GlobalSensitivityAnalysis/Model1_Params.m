% ------------------------------------------------------------------------
% DEFINING PARAMETERS OF MODEL 1
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 03/11/2020
% ------------------------------------------------------------------------

function params = Model1_Params(arch)


% Load optimally fitted parameter values:
load('Params_OptFit_NAR.mat') %#ok<*LOAD>
p = Params_OptFit_NAR;

% FadR expression and regulation:
  b_R   = p(2);  % Basal leaky expression of FadR.
  a_R   = p(3);  % Maximum rate of expression of FadR.
  K_R   = p(4);  % Affinity (ratio of rate of binding to unbinding) of FadR to its promoter.

% FadD expression and regulation:
  b_D   = p(6);  % Basal leaky expression of FadD.
  a_D   = p(7);  % Maximum rate of expression of FadD.
  K_D   = p(8);  % Affinity (ratio of rate of binding to unbinding) of FadR to FadD promoter.

% Uptake Reaction (single pseudo-rxn representing FadL and FadD):
  kcatD = p(10); % Turnover rate of enzyme of FadD.
  Km_D  = p(11); % Michaelis-constant of enzyme of FadD for Fatty Acid (OA).
  
% PlsB enzyme kinetics - Glycerol-3-Phosphate Acyltransferase:
  kcatB = p(12); % Turnover rate of enzyme of PlsB.
  Km_B  = p(13); % Michaelis-constant of enzyme of PlsB for Acyl-CoA.
  B     = p(14); % Assumed constant levels of PlsB enzyme.

% Sequestration of FadR by Acyl-CoA:
  kf    = p(15); % Effective forward rate of FadR sequestration by Acyl-CoA.
  kr    = p(16); % Effective reverse rate of sequestering.

% Eg expression and regulation:
  a_g   = p(1); % Maximum synthesis rate of enzyme Eg.
  K_g   = 4.9114; % Affinity of FadR to activate Eg synthesis, based on FadR PAR promoter.
  
% Growth rate:
  muMax = p(1); % Assumed constant, as seen in experiments of Di Liu.
  sT    = 0.7; % Minimum growth rate when no Ep expressed, based on Usui 2012 report that mu is reduced to 30% of its value when pgi is knocked out (Usui 2013), and growth in our cell is at 0.1818 1/h.

% Defining the vector of parameters:
  p2 = [b_R;   % 1
        a_R;   % 2
        K_R;   % 3
        b_D;   % 4
        a_D;   % 5
        K_D;   % 6
        kcatD; % 7
        Km_D;  % 8
        kcatB; % 9
        Km_B;  % 10
        B;     % 11
        kf;    % 12
        kr;    % 13
        a_g;   % 14
        K_g;   % 15
        muMax; % 16
        sT];   % 17

% Redefine parameters for circuit topologies with NAR and PAR, so that all
% three synthesis same steady state conc of FadR in absence of inducer:
    params = p2;
if strcmp(arch,'NAR') == 1
    % do nothing
elseif strcmp(arch,'PAR') == 1
    params(3) = p2(3) * 7;
end
