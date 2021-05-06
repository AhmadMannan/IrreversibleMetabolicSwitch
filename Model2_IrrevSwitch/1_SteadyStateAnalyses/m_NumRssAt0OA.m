% ------------------------------------------------------------------------
%  DETERMINING NUMBER OF STEADY STATES OF FADR WHEN OA = 0uM.
%  Model of irreversible switch
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 17/11/2020
% ------------------------------------------------------------------------

function NumSS = m_NumRssAt0OA(params,comp)


% Model parameters:
p     = params;
p     = params;
b_R   = p(1);  % Basal leaky expression of FadR.
a_R   = p(2);  % Maximum rate of expression of FadR.
K_R   = p(3);  % Ratio of rate of binding to unbinding of FadR to its promoter.
b_D   = p(4);  % Basal leaky expression of FadD.
a_D   = p(5);  % Maximum rate of expression of FadD.
K_D   = p(6);  % Ratio of rate of binding to unbinding of FadR to FadD promoter.
kcatD = p(7); % Turnover rate of enzyme of FadD.
Km_D  = p(8); % Michaelis-constant of enzyme of FadD for Fatty Acid (OA).
kcatB = p(9); % Turnover rate of enzyme of PlsB.
Km_B  = p(10); % Michaelis-constant of enzyme of PlsB for Acyl-CoA.
B     = p(11); % Assumed constant levels of PlsB enzyme.
kf    = p(12); % Effective forward rate of FadR sequestration by Acyl-CoA.
kr    = p(13); % Effective reverse rate of sequestering.
a_g   = p(14); % Maximum synthesis rate of Eg (enzyme for growth).
K_g   = p(15); % Affinity of FadR to activating expression of Eg.
muMax = p(16); % Maximum achievable growth rate (as set).
sT    = p(17); % Severity of growth-production enzyme concentration trade-off, ranging from 0 (no trade-off) to 1 (most severe trade-off).
b_T   = p(18); % TetR expression leakiness
a_T   = p(19); % TetR promoter strength
K_Ri  = p(20); % Affinity of FadR to inhibit TetR expression
K_T   = p(21); % Affinity of TetR to inhibit FadR expression

% Range of FadR values to generate curve over:
R = logspace(-10,0,200);

% Calculate a curve over many FadR (R) values:
mu = (muMax*(1 - sT) + sqrt(((muMax*(1 - sT))^2) + (4*muMax*sT*a_g*K_g*R./(1 + (K_g*R)))))./2;
T_R = (1./mu) .* (b_T + a_T./(1 + (K_Ri*R).^2));
if strcmp(comp,'COMP') == 1
    P_R = a_R * K_R .* R ./ (1 + K_R.*R + ((K_T*T_R).^2));
elseif strcmp(comp,'NONCOMP') == 1
    P_R = a_R * K_R .* R ./ ((1 + K_R.*R).*(1 + (K_T*T_R).^2));
elseif strcmp(comp,'NOPAR') == 1
    P_R = a_R ./ (1 + (K_T*T_R).^2);
end
curve = b_R + P_R - (mu .* R);

% Finding intervals where curve crosses zero, to determine number of
% steady states:
cs = sign(curve); % Convert all points on curve to -1 or +1.
signChgs = 0; % Sign change counter
s_t = cs(1);
for j = 1:length(cs)
    if cs(j) ~= s_t
        signChgs = signChgs + 1;
        s_t = cs(j);
    end
end
NumSS = signChgs;
