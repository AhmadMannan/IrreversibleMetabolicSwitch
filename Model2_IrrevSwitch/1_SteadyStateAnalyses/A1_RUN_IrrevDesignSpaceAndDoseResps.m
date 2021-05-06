% ------------------------------------------------------------------------
% SHOW DESIGN SPACE IN WHICH THREE DIFFERENT SYSTEMS CAN BE IRREVERSIBLE,
% where the systems include (i) w/o PAR (NOPAR), (ii) PAR +toggle with
% competitive binding of FadR & TetR for fadR promoter (COMP), and (iii)
% with non-competitive binding to fadR promoter (NONCOMP).
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 19/11/2020
% ------------------------------------------------------------------------

function p_IS = A1_RUN_IrrevDesignSpaceAndDoseResps


%%
% --- [A] Irreversible design space --------------------------------------

% Determine a select parameter regime in which the following systems are
% irreversible:
circuits = {'COMP','NONCOMP','NOPAR'};
P_mss = []; % parameter regimes with multiple steady states
Rng = [];
for i = 1:3
    tic
    [P_mss{i},Rng{i}] = m_FindParams4IrrevSwitch(circuits{i},i); %#ok<*AGROW>
    toc
end

% save('P_mss','P_mss');

%
% Selecting one of the parameter regimes for irreversible switch (IS):
params = Model2_Params; % set of nominal param values
SelParTuning = P_mss{3}(round(length(P_mss{3}(:,1))*0.23),:);
disp(SelParTuning)
SelParTuning

% Plotting our selected parameter regime for irreversible switch on figs:
numPs = length(P_mss{1}(1,:));
sppos = zeros(numPs); sppos(1:numPs^2) = 1:numPs^2; sppos = sppos'; % subplot positions
for c = 1:3 % the 3 diff circuit topologies
    figure(c)
    for i = 1:numPs-1
        for j = i+1:numPs
            % Go to correct subplot:
            subplot(numPs,numPs,sppos(i,j))
            hold on
            plot(find(SelParTuning(j) == Rng{c}{j}),find(SelParTuning(i) == Rng{c}{i}),'ws','MarkerSize',6,'MarkerFaceColor','g');
            hold off
        end
    end
end


%%
% --- [B] Plot of dose-response of four systems with same param regime ---

% [ii] Select a parameter regime for irreversibility and plot dose-resp of
%      the four systems with only PAR, that w/o PAR and those with PAR and
%      COMP and NONCOMP binding.
% ------------------------------------------------------------------------

% Select one of the parameter regimes for irreversible switch (IS):
psc = ones(size(params)); % parameter scaling
psc([2,19,20,21]) = SelParTuning;
p_IS = params .* psc; % params for irreversible switch

% Save parameters of engineered, irrversible switch with Comp topology:
Params_Eng_IrrevSw = p_IS;
% save('Params_Eng_IrrevSw','Params_Eng_IrrevSw')

p_IS = Params_Eng_IrrevSw;
% Generate dose-response for systems with ...
figure(4); clf
% rangeOA = logspace(-4.5,-1.5,1000);
rangeOA = logspace(-4.5,0,1000);
set(gca,'XScale','log','YScale','log')
xlabel('Oleic acid (\mu M)')
ylabel('FadR (\mu M)')
  % ... only PAR (stable states):
    p_t = p_IS; p_t([18,19,21]) = 0;
    SSvsOA = m_DoseResp(p_t,rangeOA','COMP',[]);
    [~,idx] = sort(SSvsOA(:,2),'descend');
    FadRss_PAR = SSvsOA(idx,[1,2,end]);
    unstPts_PAR = find(FadRss_PAR(:,end) == 0);
    hold on
    plot(FadRss_PAR(:,1),FadRss_PAR(:,2),'-','Color',[0.1,0.4,0.8],'LineWidth',3)
    hold off

  % ... w/o PAR but with only toggle (stable states):
    p_t = p_IS;
    SSvsOA = m_DoseResp(p_t,rangeOA','NOPAR',[]);
    [~,idx] = sort(SSvsOA(:,2),'descend');
    FadRss_NOPAR = SSvsOA(idx,[1,2,end]);
    unstPts_NOPAR = find(FadRss_NOPAR(:,end) == 0);
    hold on
    plot(FadRss_NOPAR(:,1),FadRss_NOPAR(:,2),'-','Color',[0.7,0.7,0.7],'LineWidth',3)
    hold off
    
  % ... PAR with COMP binding of FadR and TetR (stable states):
    p_t = p_IS;
    SSvsOA = m_DoseResp(p_t,rangeOA','COMP',[]);
    [~,idx] = sort(SSvsOA(:,2),'descend');
    FadRss_COMP = SSvsOA(idx,[1,2,end]);
    unstPts_COMP = find(FadRss_COMP(:,end) == 0);
    hold on
    plot(FadRss_COMP(:,1),FadRss_COMP(:,2),'-','Color',[0,0,0],'LineWidth',3)
    hold off
    
  % ... PAR with NONCOMP binding of FadR and TetR (stable states):
    p_t = p_IS;
    SSvsOA = m_DoseResp(p_t,rangeOA','NONCOMP',[]);
    [~,idx] = sort(SSvsOA(:,2),'descend');
    FadRss_NONCOMP = SSvsOA(idx,[1,2,end]);
    unstPts_NONCOMP = find(FadRss_NONCOMP(:,end) == 0);
    hold on
    plot(FadRss_NONCOMP(:,1),FadRss_NONCOMP(:,2),'-','Color',[0.5,0.5,0.5],'LineWidth',3)
    hold off
    
  % Plotting unstable states of ...
    hold on
    plot(FadRss_PAR(unstPts_PAR,1),FadRss_PAR(unstPts_PAR,2),'w:','LineWidth',3.5)                 % .. PAR
    plot(FadRss_NOPAR(unstPts_NOPAR,1),FadRss_NOPAR(unstPts_NOPAR,2),'w:','LineWidth',3.5)         % .. NoPAR
    plot(FadRss_COMP(unstPts_COMP,1),FadRss_COMP(unstPts_COMP,2),'w:','LineWidth',3.5)             % .. COMP
    plot(FadRss_NONCOMP(unstPts_NONCOMP,1),FadRss_NONCOMP(unstPts_NONCOMP,2),'w:','LineWidth',3.5) % .. NonCOMP
    hold off
    
    axis([min(rangeOA), max(rangeOA), 1e-6, 10])
    legend('PAR-stable','NoPAR-stable','COMP-stable','NonCOMP-stable','PAR-unstable','NoPAR-unstable','COMP-unstable','NonCOMP-unstable','Location','SouthEast')
    