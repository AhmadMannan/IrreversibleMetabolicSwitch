% ------------------------------------------------------------------------
% PARAMETER SETTING OF FA MODEL, FOR GLOBAL SENSITIVITY ANALYSIS
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% ------------------------------------------------------------------------

function [p,p_var_pos,efast_var,pmin,pmax,time_points,y0,y_var_label] = m_efast_ParamSettings(arch,metric)


% --- PARAMETERS ---------------------------------------------------------

% Baseline parameter values:
p = Model1_Params(arch);

% Expression of FadR:
  b_R   = p(1);  % Basal leaky expression of FadR.
  a_R   = p(2);  % Maximum rate of expression of FadR.
  K_R   = p(3);  % Ratio of rate of binding to unbinding of FadR to its promoter.
% Expression of FadD:
  b_D   = p(4);  % Basal leaky expression of FadD.
  a_D   = p(5);  % Maximum rate of expression of FadD.
  K_D   = p(6);  % Ratio of rate of binding to unbinding of FadR to FadD promoter.
% Uptake Reaction (single pseudo-rxn representing FadL and FadD):
  % kcatD = p(7); % Turnover rate of enzyme of FadD.
  Km_D  = p(8); % Michaelis-constant of enzyme of FadD for Fatty Acid (OA).
% PlsB enzyme kinetics - Glycerol-3-Phosphate Acyltransferase:
  % kcatB = p(9); % Turnover rate of enzyme of PlsB.
  Km_B  = p(10); % Michaelis-constant of enzyme of PlsB for Acyl-CoA.
  % B     = p(11); % Assumed constant levels of PlsB enzyme.
% Sequestration of FadR by Acyl-CoA:
  kf    = p(12); % Effective forward rate of FadR sequestration by Acyl-CoA.
  kr    = p(13); % Effective reverse rate of sequestering.
% Severity of growth attenuation during induction:
  sT    = p(17);
% dummy parameter:
  dummy = 1;

% Define parameters to vary for sensitivity analysis:
% params = [b_R, a_R, K_R, b_D, a_D, K_D, Km_D, Km_B, kf, kr, sT, dummy];
% p_var_pos = [1,2,3,4,5,6,8,10,12,13,17,18];
params = [b_R, a_R, K_R, b_D, a_D, K_D, sT, dummy];
p_var_pos = [1,2,3,4,5,6,17,18];

% Parameter labels:
% efast_var = {'b_R','a_R','K_R','b_D','a_D','K_D','Km_D','Km_B','k_f','k_r','s_T','dummy'};
efast_var = {'b_R','a_R','K_R','b_D','a_D','K_D','s_T','dummy'};

% Define min and max vectors of parameters:
pmin = params * 0.1; % 0.1
pmax = params * 5; % 10
pmin(end-1) = 0.1; % for param sT
pmax(end-1) = 0.9; % for param sT


% --- INPUTS FOR SOLVING MODEL -------------------------------------------

% Simulation time span:
time_points = 1;

% Variables Labels
y0 = rand(1); % dummy initial condition
if metric == 'BR' %#ok<*BDSCA>
    y_var_label = {'Bistable range (BR)'};
elseif metric == 'IT'
    y_var_label = {'Induction threshold (IT)'};
elseif metric == 'RT'
    y_var_label = {'Reversion threshold (RT)'};
end
