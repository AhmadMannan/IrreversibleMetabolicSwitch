% ------------------------------------------------------------------------
%     DETERMINE MAXIMUM ACHIEVABLE FADR IN THE ABSENCE OF INDUCING OA
%             Continuous culture of system controlling growth
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 04/11/2020
% ------------------------------------------------------------------------

function [Rmax,fval,exitflag] = m_MaxSSFadR(params,arch)

[Rmax,fval,exitflag] = fminbnd(@(R)DetMaxSSFadR(R),0,1000);


function Obj = DetMaxSSFadR(R)


% ---- Define parameters -------------------------------------------------

p     = params;
b_R   = p(1);  %#ok<*NASGU> % Basal leaky expression of FadR.
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
muMax = p(16); % Factor converting number of Eg into growth rate.
sT    = p(17); % Growth rate reduced to 30% of its value when Eg = 0.


% ---- Define growth as a function of FadR -------------------------------

mu = (muMax*(1-sT) + sqrt(((muMax*(sT-1))^2) + (4*muMax*sT*a_g*K_g*R/(1 + (K_g*R)))))/2;


% ---| Define objective whose minimum will determine |--------------------
% ---| the max steady state FadR for given params.   |--------------------

if strcmp(arch,'NAR') == 1
    P_R = a_R/(1 + K_R*R);
elseif strcmp(arch,'PAR') == 1
    P_R = a_R*K_R*R/(1 + K_R*R);
end

Obj = abs(b_R + P_R - R*mu);

end

end