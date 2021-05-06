% ------------------------------------------------------------------------
% DETERMINE THE NUMBER OF STEADY STATES AT EACH INDUCER CONCENTRATION
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 19/03/2021
% ------------------------------------------------------------------------

function numSS = m_Det_MonoBistIrrev(params)


% Determine two dose response curves: for increasing and decreasing doses
% of inducer:
rngInd = [0,logspace(-4,3,10)];
[dr_inc,dr_dec] = m_DoseResp(params,rngInd);

% Determine difference between steady state levels of LacI in both dose
% response curves:
% [~,idx_i] = sort(dr_inc(:,1),'ascend');
% [~,idx_d] = sort(dr_dec(:,1),'ascend');
% dr_inc = dr_inc(idx_i,:);
% dr_dec = dr_dec(idx_d,:);

diff = (dr_inc(:,1+3) - dr_dec(end:-1:1,1+3))./dr_inc(:,1+3);
% diff = (dr_inc(:,1+3) - dr_dec(:,1+3))./dr_inc(:,1+3);

if max(abs(diff)) > 0.1 % more than 10% difference
    if abs(diff(1)) > 0.1
        numSS = 3;
    else
        numSS = 2;
    end
else
    numSS = 1;
end
