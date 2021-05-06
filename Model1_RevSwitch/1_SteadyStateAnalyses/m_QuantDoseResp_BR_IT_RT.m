% ------------------------------------------------------------------------
% GENERATE DOSE RESPONSE AND QUANTIFY BISTABLE RANGE (BR), INDUCTION
% THRESHOLD (IT) AND REVERSION THRESHOLD (RT), FOR VARIATIONS IN PARAMS
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 05/11/2020
% ------------------------------------------------------------------------

function [BR,IT,RT,ssR] = m_QuantDoseResp_BR_IT_RT(params,arch,fig)


%%
% --- Defining the key inputs --------------------------------------------

% arch = 'NAR';
% params = Model1_Params(arch);
rangeOA = logspace(-5,-1,2000);


% --- Generate the dose response for given params ------------------------
SSsVsOA = m_DoseResp(params,arch,rangeOA,[],[]);
ssR = [SSsVsOA(:,[1,2]),SSsVsOA(:,end)];

% Reorder results (in descending order of steady state FadR value):
[~,idx] = sort(ssR(:,2),'descend');
ssR = ssR(idx,:);

% Finding unstable steady states:
unst_p = find(ssR(:,end) == 0);


% --- Plotting -----------------------------------------------------------

if isempty(fig) == 0
    figure(fig); clf
    plot(ssR(:,1),ssR(:,2),'k-','LineWidth',3) % plotting all FadR steady states
    hold on
    plot(ssR(unst_p,1),ssR(unst_p,2),'r-','LineWidth',3) % highlighting unstable steady states
    hold off
    set(gca,'XScale','log')
    xlabel('OA (\mu M)')
    ylabel('FadR (\mu M)')
end


% --- Calculating the dose-response characteristics BR, IT, RT -----------

if isempty(unst_p) == 0
    % Defining the steady state FadR and OA inducing levels where the
    % system is unstable:
    ssR_unst = ssR(unst_p,:);
    
    % Measuring the approximate ...
    % ... induction threshold:
    IT = max(ssR_unst(:,1));
    % ... reversion threshold:
    RT = min(ssR_unst(:,1));
    % ... bistable range:
    % BR = max(ssR_unst(:,1)) - min(ssR_unst(:,1));
    BR = IT - RT;
else
    pos_TH = min(find(ssR(:,2) <= (max(ssR(:,2)) + min(ssR(:,2)))/2)); %#ok<*MXFND> % finding position of data point of dose-response where expression is half-way between max and min, i.e. threshold value.
    IT = ssR(pos_TH,1);
    RT = IT;
    BR = IT - RT;
end
