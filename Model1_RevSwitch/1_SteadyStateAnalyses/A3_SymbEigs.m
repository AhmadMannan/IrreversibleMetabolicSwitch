
function Eigs = A3_SymbEigs(arch)


% Define symbolic terms:
syms b_R a_R K_R b_D a_D K_D kcatD Km_D kcatB Km_B B kf kr a_g K_g sT muMax R D A C Eg OA mu0


% --- Expression of Jacobian elements ------------------------------------

% f1 = dR/dt, f2 = dD/dt, f3 = dA/dt, f4 = dC/dt, f5 = dEg/dt


% ----- f1 -----
% df1/dR:
if strcmp(arch,'NAR') == 1
    s = -1;
elseif strcmp(arch,'PAR') == 1
    s = 1;
end
J(1,1) = (s*(a_R*K_R/(1 + K_R*R)^2)) - (kf*(A^2)) - (muMax*(Eg*sT - sT + 1));

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
J(2,1) = - a_D * 2 * (K_D^2) * R / ((1 + (K_D*R)^2)^2);

% df2/dD:
J(2,2) = - muMax * (Eg*sT - sT + 1);

% df2/dA:
J(2,3) = 0;

% df2/dC:
J(2,4) = 0;

% df2/dEg:
J(2,5) = - D * muMax * sT;


% ----- f3 -----
% df3/dR:
J(3,1) = - 2 * kf * (A^2);

% df3/dD:
J(3,2) = kcatD * OA / (Km_D + OA);

% df3/dA:
J(3,3) = - (B*Km_B*kcatB/((Km_B + A)^2)) - (4 * kf * A * R) - (muMax*(Eg*sT - sT + 1));

% df3/dC:
J(3,4) = 2 * kr;

% df3/dEg:
J(3,5) = - A * muMax * sT;


% ----- f4 -----
% df4/dR:
J(4,1) = kf * (A^2);

% df4/dD:
J(4,2) = 0;

% df4/dA:
J(4,3) = 2 * kf * A * R;

% df4/dC:
J(4,4) = - kr - (muMax*(Eg*sT - sT + 1));

% df4/dEg:
J(4,5) = - C * muMax * sT;


% ----- f5 -----
% df5/dR:
J(5,1) = a_g*K_g/((1 + K_g*R)^2);

% df5/dD:
J(5,2) = 0;

% df5/dA:
J(5,3) = 0;

% df5/dC:
J(5,4) = 0;

% df5/dEg:
J(5,5) = - (2*muMax*Eg*sT) + (muMax*sT) - muMax;


% --- Substituting steady state w/o inducer ------------------------------
% To calculate the convergence rate close to the uninduced steady state of
% the system, we (i) evaluate the Jacobian at the uninduced steady state,
% defined as (R,D,A,C,Eg)=(R,D,0,0,1), and (ii) determine its eignevalues:

% (i) evaluate Jacobian at induced steady state:
J = subs(J,OA,0);
J = subs(J,A,0);
J = subs(J,C,0);
J = subs(J,Eg,1);
J = subs(J,sT,0);
J = subs(J,muMax,mu0);


% (ii) determine its eigenvalues:
Eigs = eig(J);
