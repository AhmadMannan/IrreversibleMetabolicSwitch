% ------------------------------------------------------------------------
%           PLOTTING CHANGE IN STEADY STATE DOSE-RESPONSE CURVE
%         FOR CHANGES IN PARAMETERS FOR CIRCUITS WITH NAR AND PAR
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 04/11/2020
% ------------------------------------------------------------------------

function A1_RunLocalSensitivityAnalysis


%% Run analysis and plots dose-response curves for varying parameters for
% circuit with ...
% ... NAR:
A1_DoseResp_VaryParams('NAR',1)

% ... PAR:
A1_DoseResp_VaryParams('PAR',2)
