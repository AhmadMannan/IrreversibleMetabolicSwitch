% ------------------------------------------------------------------------
% STUDY HOW TUNING OA FEED-IN AFFECTS TOTAL INDUCER USED AND SWITCH TIME
% FOR CIRCUITS WITH NAR AND PAR
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 12/11/2020
% ------------------------------------------------------------------------


function A1_Run


%% -----------------------------------------------------------------------
% Generate Figure 3(b), i.e. plot of total inducer used vs switch time ...
% ------------------------------------------------------------------------

% Defining induction regimes:
rng_OAin = logspace(log10(0.25),log10(5),50); % ... range of OA feed-in flux
SigOS = 0.3; % Define fraction of 1/Kp for which induction feed-in should stop:

% ... for circuit with NAR:
p_NAR = Model1_Params('NAR');
[tOA_NAR,st_NAR] = A1_VarOAFeedin('NAR',p_NAR',rng_OAin,SigOS,[],[]);

% ... for circuit with PAR:
p_PAR = Model1_Params('PAR');
[tOA_PAR,st_PAR] = A1_VarOAFeedin('PAR',p_PAR',rng_OAin,SigOS,[],[]);

% Plotting the pareto front for both circuits:
figure(2); clf
for i = 1:2 % for NAR and PAR
    if i == 1 %#ok<*MATCH2>
        tOA = tOA_NAR;
        st = st_NAR;
        colr = autumn(length(tOA));
        colr = colr(end:-1:1,:);
    elseif i == 2
        tOA = tOA_PAR;
        st = st_PAR;
        colr = winter(length(tOA));
        colr = colr(end:-1:1,:);
    end
    
    % plot line plot first
    hold on
    plot(tOA,st,'-','LineWidth',2,'Color',colr(end,:))
    hold off
    
    for j = 1:length(tOA)
        hold on
        plot(tOA(j),st(j),'o','MarkerSize',6,'Color',colr(j,:),'MarkerFaceColor',colr(j,:))
        hold off
    end
end
xlabel('Total oleic acid used (\mu M)')
ylabel('Switch time (h)')
axis([0,500,0,90])

% % Determining the optimal regime as that closest to point (0,0):
%     % ... for circuit with NAR:
%     normDis_N = sqrt((tOA_NAR./500).^2 + (st_NAR./90).^2);
%     Pos_MinDist_N = find(min(normDis_N) == normDis_N);
%     optReg_4Nar = [tOA_NAR(Pos_MinDist_N),st_NAR(Pos_MinDist_N)];
%       % highlighting this point on the plot:
%       hold on
%       plot(optReg_4Nar(1),optReg_4Nar(2),'ro','MarkerSize',8,'LineWidth',3)
%       hold off
%     
%     % ... for circuit with PAR:
%     normDis_P = sqrt((tOA_PAR/500).^2 + (st_PAR/90).^2);
%     Pos_MinDist_P = find(min(normDis_P) == normDis_P);
%     optReg_4Par = [tOA_PAR(Pos_MinDist_P),st_PAR(Pos_MinDist_P)];
%       % highlighting this point on the plot:
%       hold on
%       plot(optReg_4Par(1),optReg_4Par(2),'bo','MarkerSize',8,'LineWidth',3)
%       hold off

      
% -----------------------------------------------------------------------
% Generate Figure 3(c), i.e. time course dynamics for an optimal regime
% ------------------------------------------------------------------------

% % Simulating the optimal induction regime ...
%     % ... for the circuit with NAR:
%     p_N = [p_NAR; rng_OAin(Pos_MinDist_N); SigOS];
%     [t_n,x_n] = m_SimIndDyn('NAR',3,p_N,100);
%     
%     % ... for the circuit with PAR:
%     p_P = [p_PAR; rng_OAin(Pos_MinDist_P); SigOS];
%     [t_p,x_p] = m_SimIndDyn('PAR',4,p_P,100);

% Highlight a selected induction regime on the plot of Obj1 vs Obj2:
    posIR = 28; % position of inudction regime
    OAfeedin = rng_OAin(posIR);
    SigOS = 0.3;
    hold on
    plot(tOA_NAR(posIR),st_NAR(posIR),'ro','MarkerSize',8,'LineWidth',3)
    plot(tOA_PAR(posIR),st_PAR(posIR),'bo','MarkerSize',8,'LineWidth',3)
    hold off
    
% Simulating time courses for both circuits, for given induction regime:
    p_N = [p_NAR; OAfeedin; SigOS];
    p_P = [p_PAR; OAfeedin; SigOS];
    [t_n,x_n] = m_SimIndDyn('NAR',3,p_N,100);
    [t_p,x_p] = m_SimIndDyn('PAR',4,p_P,100);
    
% Plotting the time courses for ...
    % ... circuit with NAR:
    figure(3); clf
    subplot(1,3,1) % of FadR and oleic acid
        plot(t_n,x_n(:,1),'r-','LineWidth',2)
        hold on
        plot(t_n,x_n(:,6),'k-','LineWidth',2)
        plot([0,100],[1/p_N(18),1/p_N(18)],'k--')
        hold off
        xlabel('Time (h)'); ylabel('Concentration (\mu M)')
        legend('FadR','OA','Location','NorthEast')
        title('Circuit with NAR')
        set(gca,'YScale','log')
        ylim([5e-4,1e-1])
    subplot(1,3,2) % of total FadR
        plot(t_n,x_n(:,1)+x_n(:,4),'k-','LineWidth',2)
        hold on
        plot(t_n,x_n(:,4),'k--','LineWidth',2)
        hold off
        xlabel('Time (h)'); ylabel('Concentration (\mu M)')
        legend('Total FadR','Complex','Location','SouthEast')
    subplot(1,3,3) % of acyl-CoA
        plot(t_n,x_n(:,3),'k-','LineWidth',2)
        xlabel('Time (h)'); ylabel('Acyl-CoA (\mu M)')
    
    % ... circuit with PAR:
    figure(4); clf
    subplot(1,3,1) % of FadR and oleic acid
        plot(t_p,x_p(:,1),'b-','LineWidth',2)
        hold on
        plot(t_p,x_p(:,6),'k-','LineWidth',2)
        plot([0,100],[1/p_P(18),1/p_P(18)],'k--')
        hold off
        xlabel('Time (h)'); ylabel('Concentration (\mu M)')
        legend('FadR','OA','Location','NorthEast')
        title('Circuit with PAR')
        set(gca,'YScale','log')
        ylim([5e-4,1e-1])
    subplot(1,3,2) % of total FadR
        plot(t_p,x_p(:,1)+x_p(:,4),'k-','LineWidth',2)
        hold on
        plot(t_p,x_p(:,4),'k--','LineWidth',2)
        hold off
        xlabel('Time (h)'); ylabel('Concentration (\mu M)')
        legend('Total FadR','Complex','Location','NorthEast')
    subplot(1,3,3) % of acyl-CoA
        plot(t_p,x_p(:,3),'k-','LineWidth',2)
        xlabel('Time (h)'); ylabel('Acyl-CoA (\mu M)')

