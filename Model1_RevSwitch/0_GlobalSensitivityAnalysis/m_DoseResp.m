% ------------------------------------------------------------------------
%           KINETIC MODEL STEADY STATES FOR VARYING OA INDUCTION
%             Continuous culture of system controlling growth
% ------------------------------------------------------------------------
% Generate dose-response curves for titrations of OA, and identify stable,
% unstable steady states for bistability (hysteresis).
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 04/11/2020
% ------------------------------------------------------------------------

function SSsVsOA = m_DoseResp(params,arch,rangeOA,checkBistability,plotID)


%% ------------------------------------------------------------------------
% Define global settings
% ------------------------------------------------------------------------
% rangeOA = logspace(-3,-1,100);
OA = rangeOA; % Titrations of inducing OA.
fig = 1;      % Figure and plotting settings.
SSsVsOA = []; % Matrix of steady states of all species vs OA
R_ss = [];    % Steady states of only free FadR vs OA


% ------------------------------------------------------------------------
% Find steady states for each value of inducing OA
% ------------------------------------------------------------------------

% w = waitbar(0,'Please wait ...');
for i = 1:length(OA)
    % --- Determining range of 'R' to generate curve over ----------------
    
    % Max value of R steady state (i.e. where nullcline1 intercepts A = 0):
    Rmax = m_MaxSSFadR(params,arch);
    
    % Defining range of 'R' to generate nullclines:
    R = linspace(1e-20,Rmax*1.1,500);
    
    
    % --- Generating master nullcline (MN) -------------------------------
    
    [~,MN] = Model1_Obj_NullclinesAndSS(R,params,OA(i),arch);
%     figure(7)
%     if i == 1
%         clf
%     end
%     hold on
%     plot(R,real(MN),'o','MarkerSize',3)
%     plot(R,zeros(size(R)),'k-')
%     ylim([-2,2])
%     hold off
    
    % --- Determining steady state points --------------------------------
    % Note: The steady states are determined as the intersection between
    % y = 0 and the MN.
    
    % Finding intervals where MN and y = 0 intersect:
    MNsign = sign(real(MN));
    c = MNsign(1);
    pos_chg = [];
    for j = 1:length(MNsign)
        if MNsign(j) ~= c
            pos_chg = [pos_chg; j-1 j]; %#ok<*AGROW>
            c = MNsign(j);
        end
    end
    rng_R = [];
    for j = 1:length(pos_chg(:,1))
        rng_R = [rng_R; sort(R(pos_chg(j,:)),'ascend')];
    end
    
    % Break loop if we find system is bistable:
    % Note: We only want to do this when determining region of parameter
    %       space where system is bistable and monostable.
    if isempty(checkBistability) == 0
        if length(pos_chg(:,1)) >= 2
            break
        end
    end
    
%     % Steady states of FadR (R):
%     options = optimset('Display','final','TolFun',1e-10,'TolX',1e-100,'MaxFunEval',1e4,'MaxIter',1e4);
%     rss_t = [];
%     for j = 1:length(pos_chg(:,1))
%         lb = rng_R(j,1);
%         ub = rng_R(j,2);
%         [R_ss_t,fval] = fminbnd(@(rss)Model_FAU_MetEng_Obj_NullclinesAndSS(rss,params,OA(i),arch),lb,ub,options);
%         % [R_ss_t,fval] = fmincon(@(rss)Model_FAU_MetEng_Obj_NullclinesAndSS(rss,params,OA(i),arch),lb,[],[],[],[],lb,ub,[],options);
%         if abs(fval) < 1e-1
%             % R_ss = [R_ss; OA(i) R_ss_t];
%             rss_t = [rss_t; R_ss_t];
%         end
%     end
    
    % Steady states of FadR (R):
    rss_t = [];
    for j = 1:length(pos_chg(:,1))
        lb = rng_R(j,1);
        ub = rng_R(j,2);
        [R_ss_t,fval] = m_Bisection(@(rss)Model1_Obj_NullclinesAndSS(rss,params,OA(i),arch),lb,ub);
        if abs(fval) <= 1e-4
            % R_ss = [R_ss; OA(i) R_ss_t];
            rss_t = [rss_t; R_ss_t];
        end
    end
    
    % Steady states and stabilities:
    for j = 1:length(rss_t)
        % Vector of steady state:
        [~,~,SS_t] = Model1_Obj_NullclinesAndSS(rss_t(j),params,OA(i),arch);
        
        % Checking stability of steady state (eigenvalues of Jacobian):
        J_t = Model1_Jacobian(params,arch,SS_t,OA(i));
        ev_t = eig(J_t);
        if max(real(ev_t)) <= 0
            stable = 1;
        else
            stable = 0;
        end
        
        % Saving steady states of all variables:
        SSsVsOA = [SSsVsOA; OA(i),SS_t,stable];
        R_ss = [R_ss; OA(i),SS_t(1),stable];
    end
    
    % Condition to terminate iterations earlier, once we have reached the
    % FadR off-state:
    DiffVs1st = R_ss(end,2)/R_ss(1,2);
    if isempty(checkBistability) == 0 && DiffVs1st < 0.5
        diff = abs(R_ss(end-1,2) - R_ss(end,2));
        if diff < 1e-4
            break
        end
    end
    
%     waitbar(i/length(OA))
end
% close(w)


% ------------------------------------------------------------------------
% Plotting steady state of FadR and FadD vs OA
% ------------------------------------------------------------------------

if isempty(plotID) == 0
    % Sort FadR steady states in descending order:
    [~,idx] = sort(SSsVsOA(:,2),'descend');
    FadRss = SSsVsOA(idx,[1,2,end]);
    FadDss = SSsVsOA(idx,[1,3,end]);
    
    % Plot dose-response curve:
    figure(fig); clf

    subplot(1,2,1)
    plot(FadRss(:,1),FadRss(:,2),'k-','LineWidth',2)
    stSS = find(FadRss(:,3) == 1);
    ustSS = find(FadRss(:,3) == 0);
    hold on
    plot(FadRss(stSS,1),FadRss(stSS,2),'ko','LineWidth',1,'MarkerSize',8,'MarkerFaceColor','k')
    plot(FadRss(ustSS,1),FadRss(ustSS,2),'ro','LineWidth',1,'MarkerSize',8,'MarkerFaceColor','r')
    hold off
    xlabel('Inducing OA (\mu M)','FontSize',14)
    ylabel('Free FadR (\mu M)','FontSize',14)
    set(gca,'XScale','log')

    subplot(1,2,2)
    plot(FadDss(:,1),FadDss(:,2),'k-','LineWidth',2)
    stSS = find(FadDss(:,3) == 1);
    ustSS = find(FadDss(:,3) == 0);
    hold on
    plot(FadDss(stSS,1),FadDss(stSS,2),'ko','LineWidth',1,'MarkerSize',8,'MarkerFaceColor','k')
    plot(FadDss(ustSS,1),FadDss(ustSS,2),'ro','LineWidth',1,'MarkerSize',8,'MarkerFaceColor','r')
    hold off
    xlabel('Inducing OA (\mu M)','FontSize',14)
    ylabel('FadD (\mu M)','FontSize',14)
    set(gca,'XScale','log')
end


% --- Checking for bistability -------------------------------------------
% if isempty(checkBistability) == 0
%     if length(pos_chg(:,1)) >= 2
%         BistableInd = 1;
%     else
%         BistableInd = 0;
%     end
% else
%     BistableInd = [];
% end

