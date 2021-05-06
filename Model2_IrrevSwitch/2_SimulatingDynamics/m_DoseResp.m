% ------------------------------------------------------------------------
%    DOSE-RESPONSE CURVES - FADR STEADY STATES FOR VARYING OA INDUCTION
% ------------------------------------------------------------------------
% Generate dose-response curves for titrations of OA, and identify stable
% and unstable steady states of system with integrated toggle-switch.
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 17/11/2020
% ------------------------------------------------------------------------

function [SSsVsOA,RVsOA] = m_DoseResp(p2,rangeOA,comp,plotID)


%% -----------------------------------------------------------------------
% Define global settings
% ------------------------------------------------------------------------
% rangeOA = logspace(-3,-1,100);
OA = rangeOA; % Titrations of inducing OA.
SSsVsOA = []; % Matrix of steady states of all species vs OA
RVsOA = [];   % Steady states of only free FadR vs OA


% ------------------------------------------------------------------------
% Find steady states for each value of inducing OA
% ------------------------------------------------------------------------

% w = waitbar(0,'Please wait ...');
for i = 1:length(OA)
    % --- Determining range of 'R' to generate curve over ----------------
    
    % Defining range of 'R' to generate nullclines:
    R = logspace(-20,3,500);
    
    % --- Generating master nullcline (MN) -------------------------------
    MN = [];
    for k = 1:length(R)
        [~,mn] = m_Obj_NullclinesAndSS(R(k),p2,OA(i),comp);
        MN = [MN; mn];
    end
    % [~,MN] = m_NullclinesAndSS(R,params,OA(i),comp);
    MN = real(MN);
%     figure(10)
%     if i == 1
%         clf
%     end
%     hold on
%     plot(R,MN,'o-','MarkerSize',3)
%     plot(R,zeros(size(R)),'k-')
%     xlabel('FadR (\mu M)'); ylabel('Master Nullcline')
%     ylim([-10,10])
%     set(gca,'XScale','log')
%     hold off
    
    % --- Determining FadR steady states ---------------------------------
    % Note: FadR steady states are found as the intersection between y = 0
    %       and the MN curve.
    
    % Finding intervals where MN and y = 0 intersect:
    MNsign = sign(MN);
    c = MNsign(1);
    pos_chg = [];
    for j = 1:length(MNsign)
        if MNsign(j) ~= c
            pos_chg = [pos_chg; j-1 j]; %#ok<*AGROW>
            c = MNsign(j);
        end
    end
    rng_R = [];
    length(pos_chg(:,1));
    for j = 1:length(pos_chg(:,1))
        rng_R = [rng_R; sort(R(pos_chg(j,:)),'ascend')];
    end
    
    % Steady states of FadR (R):
    rss_t = [];
    for j = 1:length(pos_chg(:,1))
        lb = rng_R(j,1);
        ub = rng_R(j,2);
        [R_ss_t,fval] = m_Bisection(@(rss)m_Obj_NullclinesAndSS(rss,p2,OA(i),comp),lb,ub);
        if abs(fval) <= 1e-6
            rss_t = [rss_t; R_ss_t];
        end
    end
    
    % Steady states and stabilities:
    for j = 1:length(rss_t)
        % Vector of steady state:
        [~,~,SS_t] = m_Obj_NullclinesAndSS(rss_t(j),p2,OA(i),comp);
        
        % Checking stability of steady state (eigenvalues of Jacobian):
        J_t = m_Jacobian(p2,comp,SS_t,OA(i));
        ev_t = eig(J_t);
        if max(real(ev_t)) <= 0
            stability = 1;
        else
            stability = 0;
        end
        
        % Saving steady states of all variables:
        SSsVsOA = [SSsVsOA; OA(i),SS_t,stability];
        RVsOA = [RVsOA; OA(i),SS_t(1),stability];
    end
    
%     waitbar(i/length(OA))
end
% close(w)


% ------------------------------------------------------------------------
% Plotting steady state of FadR and FadD vs OA
% ------------------------------------------------------------------------

if isempty(plotID) == 0
    fig = plotID;
    
    % Sort FadR steady states in descending order:
    [~,idx] = sort(SSsVsOA(:,2),'descend');
    FadRss = SSsVsOA(idx,[1,2,end]);
    
    % Plot dose-response curve:
    figure(fig); clf
    
    hold on; plot(FadRss(:,1),FadRss(:,2),'ko','LineWidth',2); hold off
    stSS = find(FadRss(:,3) == 1);
    ustSS = find(FadRss(:,3) == 0);
    hold on; plot(FadRss(ustSS,1),FadRss(ustSS,2),'ro','LineWidth',2); hold off
    xlabel('Inducing OA (\mu M)','FontSize',14)
    ylabel('Free FadR (\mu M)','FontSize',14)
%     set(gca,'XScale','log')
end

