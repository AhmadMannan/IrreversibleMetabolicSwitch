% ------------------------------------------------------------------------
% DEFINING PARAMETERS OF MODEL 2
% Model of irreversible switch
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 16/11/2020
% ------------------------------------------------------------------------

function params = Model2_Params


% Load optimally fitted parameter values:
load('Params_OptFit_NAR.mat') %#ok<*LOAD>
p = Params_OptFit_NAR;

% FadR expression and regulation:
%   b_R   = p(2);  % Basal leaky expression of FadR.
%   a_R   = p(3);  % Maximum rate of expression of FadR.
%   K_R   = p(4);  % Affinity (ratio of rate of binding to unbinding) of FadR to its promoter.
  b_R   = 0.0005;  % Basal leaky expression of FadR.
  a_R   = 0.1409;  % Maximum rate of expression of FadR.
  K_R   = 4.9114;  % Affinity (ratio of rate of binding to unbinding) of FadR to its promoter.

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

% TetR exp and reg by FadR, and affinity of TetR to inhibit FadR exp:
  b_T   = p(6); % p(2) TetR expression leakiness - based on FadD promoter.
  a_T   = p(7); % p(3) TetR promoter strength - based on FadD promoter.
  K_Ri  = p(8); % Affinity of FadR to inhibit TetR exp - based on FadD promoter.
  K_T   = (1/(1e6/(5.6e9))); % Based on an affinity of 5.6e9 1/M affinity of TetR to tetO, from Kamionka2004.

% Defining the vector of parameters:
  p2 = [b_R;
        a_R;
        K_R;
        b_D;
        a_D;
        K_D;
        kcatD;
        Km_D;
        kcatB;
        Km_B;
        B;
        kf;
        kr;
        a_g;
        K_g;
        muMax;
        sT;
        b_T;
        a_T;
        K_Ri;
        K_T];

% Redefine parameters:
    params = p2;
% b_R:
% params(2) = p(2) * 0.6981; % * 0.0211;% * 0.0725;
% % a_R:
% params(3) = p(7) * 2.7271; % * 12.7974; % p(3) * 50.64; * 15.75; % Tuning so FadR conc before induction is same as system w/o toggle-switch
% % K_R
% params(4) = p(4) * 1.1363; % * 0.4866; %p(8) * 0.6; % * 0.4; % * 14 * 1/(14*4);
