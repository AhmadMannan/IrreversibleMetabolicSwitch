% ------------------------------------------------------------------------
% CALCULATING THE DOSE-RESPONSE
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 15/03/2021
% ------------------------------------------------------------------------

function [ss_inc,ss_dec] = m_DoseResp(params,rngInd)


% Global parameters:
tol = 1e-4;
maxiter = 100;
tpoints = 0:0.01:100;

% Define model parameters:
p = params;
p(16) = 1; % indicating constant inducer conc during dynamics, as required


% --- Find steady states for increasing inducer levels -------------------

ss_inc = [];
ic0 = [0,0,1,0,0,1]; % [Ix,I,L,C,T,Eg]
for i = 1:length(rngInd)
    
    % Define inducer conc:
    ic0(1) = rngInd(i);
    
    % Integrate till steady state (within tol):
    iter = 1;
    while iter < maxiter
        % integrate over fixed time span:
        [t_t,x_t] = ode15s(@(t,x)Model3_TS_Dyn(t,x,p),[0,max(tpoints)],ic0);
%         [t_t,x_t] = odeSun(@Model3_TS_Dyn,tpoints,ic0,p);
        
        % check max derivative value:
        dx = Model3_TS_Dyn(t_t(end),x_t(end,:),p);
        if max(abs(dx)) > tol
            ic0 = x_t(end,:);
            iter = iter + 1;
        else
            break
        end
    end
    if iter == maxiter
        disp('Did not reach steady state!')
%     elseif iter < maxiter
%         disp('Reached steady state, within tolerance!')
    end
    
    % Save steady state levels:
    ss_inc = [ss_inc; rngInd(i), x_t(end,:)]; %#ok<*AGROW>
    
    % Redefine initial conditions for next inducer level:
    ic0 = x_t(end,:);
end


% --- Find steady states for decreasing inducer levels -------------------

ss_dec = [];
rngInd = rngInd(end:-1:1);
for i = 1:length(rngInd)
    
    % Define inducer conc:
    ic0(1) = rngInd(i);
    
    % Integrate till steady state (within tol):
    iter = 1;
    while iter < maxiter
        % integrate over fixed time span:
        [t_t,x_t] = ode15s(@(t,x)Model3_TS_Dyn(t,x,p),[0,max(tpoints)],ic0);
%         [t_t,x_t] = odeSun(@Model3_TS_Dyn,tpoints,ic0,p);
        
        % check max derivative value:
        dx = Model3_TS_Dyn(t_t(end),x_t(end,:),p);
        if max(abs(dx)) > tol
            ic0 = x_t(end,:);
            iter = iter + 1;
        else
            break
        end
    end
    if iter == maxiter
        disp('Did not reach steady state!')
%     elseif iter < maxiter
%         disp('Reached steady state, within tolerance!')
    end
    
    % Save steady state levels:
    ss_dec = [ss_dec; rngInd(i), x_t(end,:)]; %#ok<*AGROW>
    
    % Redefine initial conditions for next inducer level:
    ic0 = x_t(end,:);
end


% --- Output -------------------------------------------------------------

% DR = [ss_inc; ss_dec];

