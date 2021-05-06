% ------------------------------------------------------------------------
% STUDY HOW TUNING OA FEED-IN AFFECTS TOTAL INDUCER USED AND SWITCH TIME
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 12/11/2020
% ------------------------------------------------------------------------


function [TotOA,SwitchTime,rng_OAin] = A1_VarOAFeedin(arch,params,rng_OAin,SigOS,plotID,optID)

%%
% --- Set parameters -----------------------------------------------------

% Circuit topology:
% arch = 'NAR';

% Hours of culture time:
hoursCult = 100;

% Define parameters of control system:
% p = FAU_DynCont_Params(arch);
p = params;

% Define OA feed-in rate:
% rng_OAin = logspace(log10(0.25),log10(5),50);

% Define fraction of 1/Kp for which induction feed-in should stop:
% SigOS = 0.3;


% --- Run simulations ----------------------------------------------------

TotOA = zeros(length(rng_OAin),1);
SwitchTime = zeros(length(rng_OAin),1);

% w = waitbar(0,'Please wait ...');
for i = 1:length(rng_OAin)
    % Define full set of parameters:
    p_t = [p, rng_OAin(i), SigOS];
    
    % Run simulations:
    [t,x] = m_SimIndDyn(arch,1,p_t,hoursCult);
    
    % Save values of OA feed-in flux vs time:
    OAinVsT = zeros(size(t));
    for j = 1:length(t)
        [~,OAin] = Model1_Dyn(t(j),x(j,:),p_t,arch);
        OAinVsT(j) = OAin;
    end
    
    % Record total amount of oleic acid used:
    TotOA(i) = trapz(t,OAinVsT);
    
    % Record time it took to reach 1/Kp:
    pos = min(find(x(:,1) <= (1/p(18)))); %#ok<*MXFND>
    if isempty(pos)
        SwitchTime(i) = hoursCult;
    else
        SwitchTime(i) = t(pos);
    end
    
%     waitbar(i/length(rng_OAin))
end
% close(w)


%%
if plotID == 1
    
    % Plotting
    figure(2)

    if strmatch(arch,'PAR') == 1 %#ok<*MATCH2>
        colr = winter(length(rng_OAin));
        colr = colr(end:-1:1,:);
    elseif strmatch(arch,'NAR') == 1
        colr = autumn(length(rng_OAin));
        colr = colr(end:-1:1,:);
    end

%     subplot(1,3,1)
%     for i = 1:length(rng_OAin)
%         hold on
%         plot(rng_OAin(i),TotOA(i),'o','MarkerSize',4,'Color',colr(i,:),'MarkerFaceColor',colr(i,:))
%         hold off
%     end
%     xlabel('Oleic acid feed-in rate (\mu M/h)')
%     ylabel('Total oleic acid used (\mu M)')
% 
%     subplot(1,3,2)
%     for i = 1:length(rng_OAin)
%         hold on
%         plot(rng_OAin(i),SwitchTime(i),'o','MarkerSize',4,'Color',colr(i,:),'MarkerFaceColor',colr(i,:))
%         hold off
%     end
%     xlabel('Oleic acid feed-in rate (\mu M/h)')
%     ylabel('Time for FadR to reach 1/K_p (h)')
% 
%     subplot(1,3,3)
    for i = 1:length(rng_OAin)
        hold on
        plot(TotOA(i),SwitchTime(i),'o','MarkerSize',4,'Color',colr(i,:),'MarkerFaceColor',colr(i,:))
        hold off
    end
    xlabel('Total oleic acid used (\mu M)')
    ylabel('Time for FadR to reach 1/K_p (h)')

end


%%

if optID == 1
    % Determining optimal induction regime and generating time course responses
    % of FadR and OA:

    Ob1vsObj2 = [TotOA{1}',SwitchTime{1}'];
    Ob1vsObj2(:,1) = Ob1vsObj2(:,1)/max(Ob1vsObj2(:,1));
    Ob1vsObj2(:,2) = Ob1vsObj2(:,2)/max(Ob1vsObj2(:,2));
    Dist2Utopia = zeros(size(TotOA{1}));
    for i = 1:length(TotOA{1})
        Dist2Utopia(i) = norm(Ob1vsObj2(i,:));
    end
    [~,min_p] = min(Dist2Utopia);

    opt_OAin = rng_OAin(min_p);

    % Run simulation of optimal OA feed-in rate:
    p_t = [p, opt_OAin, SigOS]; % Define full set of parameters
    fig = 3;
    figure(fig); clf
    [t,x] = Analysis_0_OAstims(arch,fig,p_t,hoursCult);
end

