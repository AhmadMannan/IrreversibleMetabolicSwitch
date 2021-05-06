% ------------------------------------------------------------------------
%              ANALYTICAL FORM OF JACOBIAN MATRIX FOR MODEL 1
%                 Continuous culture of control on growth
% ------------------------------------------------------------------------
% We write down the analytical forms of the elements of the Jacobian
% matrix, and use this to determine steady state stability.
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 04/11/2020
% ------------------------------------------------------------------------

function J = m_Jacobian(params,arch,varSS,OA)


% --- Model parameters ---------------------------------------------------

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


% --- Steady state values of variables -----------------------------------
R  = varSS(1);
D  = varSS(2);
A  = varSS(3);
C  = varSS(4);
Eg = varSS(5);


% --- Expression of Jacobian elements ------------------------------------

% f1 = dR/dt, f2 = dD/dt, f3 = dA/dt, f4 = dC/dt, f5 = dEg/dt

% Initialize 6x6 Jacobian:
J = zeros(5);

% Define growth rate:
mu = muMax * (Eg*sT - sT + 1);


% ----- f1 -----
% df1/dR:
if strcmp(arch,'NAR') == 1
    s = -1;
elseif strcmp(arch,'PAR') == 1
    s = 1;
end
J(1,1) = (s * (a_R*K_R/(1 + (K_R*R))^2)) - (kf*(A^2)) - mu;

% df1/dD:
J(1,2) = 0;

% df1/dA:
J(1,3) = - 2 * kf * A * R;

% df1/dC:
J(1,4) = kr;

% df1/dEg:
J(1,5) = - muMax * R * sT;


% ----- f2 -----
% df2/dR:
J(2,1) = - 2 * a_D * (K_D^2) * R / ((1 + (K_D*R)^2)^2);

% df2/dD:
J(2,2) = - mu;

% df2/dA:
J(2,3) = 0;

% df2/dC:
J(2,4) = 0;

% df2/dEg:
J(2,5) = - muMax * D * sT;


% ----- f3 -----
% df3/dR:
J(3,1) = - 2 * kf * (A^2);

% df3/dD:
J(3,2) = kcatD * OA / (Km_D + OA);

% df3/dA:
J(3,3) =  - (kcatB*B*Km_B/((Km_B + A)^2)) - (4*kf*A*R) - mu;

% df3/dC:
J(3,4) = 2 * kr;

% df3/dEg:
J(3,5) = - muMax * A * sT;


% ----- f4 -----
% df4/dR:
J(4,1) = kf * (A^2);

% df4/dD:
J(4,2) = 0;

% df4/dA:
J(4,3) = 2 * kf * A * R;

% df4/dC:
J(4,4) = - kr - mu;

% df4/dEg:
J(4,5) = - muMax * C * sT;


% ----- f5 -----
% df5/dR:
J(5,1) = a_g * K_g / ((1 + K_g*R)^2);

% df5/dD:
J(5,2) = 0;

% df5/dA:
J(5,3) = 0;

% df5/dC:
J(5,4) = 0;

% df5/dEg:
J(5,5) = - (2*muMax*Eg*sT) + (muMax*(sT - 1));
