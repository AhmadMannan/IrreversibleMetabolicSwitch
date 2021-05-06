% ------------------------------------------------------------------------
% DEFINE THE OBJECTIVES FOR THE MULTIOBJECTIVE OPTIMIZATION PROBLEM
%  - Obj1 = Total inducer used
%  - Obj2 = Switch time
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% ------------------------------------------------------------------------

function [Obj,t,x] = m_MultiObjs(p_IndParams,p0,comp,hoursCult)


% --- Defining params ----------------------------------------------------

params_all = p0;
params_all([23,24]) = p_IndParams;


% --- Simulate dynamic time course ---------------------------------------

[t,x] = m_TimeCourseDynSim(params_all,comp,hoursCult);


% --- Calculate obj1 = total oleic acid used -----------------------------

% Measure inducer feed-in rate at each time point of process:
Iin_t = zeros(size(t));
for j = 1:length(t)
    [~,Iin_r] = Model2_Dyn(t(j),x(j,:),params_all,comp);
    Iin_t(j) = Iin_r;
end

Obj1 = trapz(t,Iin_t);


% --- Calculate obj2 = switch time ---------------------------------------

% Low state of TF we define as the system having achieved induced state:
IndState = 0.05;

if x(end,1) <= IndState
    % ... measure approx. time taken to achieve low state:
    pos_t = min(find(x(:,1) <= IndState)); %#ok<*MXFND>
    Obj2 = t(pos_t);
else
    Obj2 = t(end);
end

% Define objective:
Obj = [Obj1, Obj2];
