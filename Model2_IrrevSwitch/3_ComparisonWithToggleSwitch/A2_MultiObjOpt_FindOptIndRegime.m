% ------------------------------------------------------------------------
% MULTIOBJECTIVE OPTIMIZATION TO FIND OPTIMAL INDUCTION REGIMES THAT
% MINIMIZE (I) TOTAL INDUCER USED, AND (II) SWITCH TIME
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 22/03/2021
% ------------------------------------------------------------------------

function A2_MultiObjOpt_FindOptIndRegime


%%
% --- Perform multiobjective optimization --------------------------------

% Define control system:
load('p_irrev')
p_nom = p_irrev; % parameters of irreversible switch
p_nom(16) = 0; % do not consider fixed inducer conc in media

% Run multiobjective optimization ...
hoursCult = 100; % length of induction/production process
lb = [1, 0.1]; % for parameters [Ifeed,Texp] (feed-in flux and exposure time)
ub = [50, 5];
options = optimoptions(@gamultiobj,'Display','iter','ParetoFraction',0.6,'PlotFcn','gaplotpareto','PopulationSize',100,'UseParallel',true);

% ... and generating the Pareto front:
tic
[IndParams_opt,ObjVals_opt,exitflag] = gamultiobj(@(p)m_MultiObjs(p,p_nom,hoursCult),length(lb),[],[],[],[],lb,ub,[],options);
toc

disp(exitflag)

% Save solutions:
% save('IndParams','IndParams')
% save('ObjVals','ObjVals')


% --- Plotting -----------------------------------------------------------

IndParams = IndParams_opt;
ObjVals = ObjVals_opt;

% Delete solutions with switch time above an acceptable time limit (10hrs):
del1 = find(ObjVals(:,2) > 99);
% del2 = find(ObjVals(:,2) < 5);
ObjVals(del1,:) = [];
IndParams(del1,:) = [];

% Sort Pareto solutions in descending order of Obj1 (total OA used):
[~,idx] = sort(ObjVals(:,1),'ascend');
ObjVals = ObjVals(idx,:);
IndParams = IndParams(idx,:);

% Save solutions:
% save('IndParams','IndParams')
% save('ObjVals','ObjVals')

figure(7); clf

% Plotting Pareto front of solutions:
subplot(3,1,1)
plot(ObjVals(:,1),ObjVals(:,2),'ko-','MarkerSize',4,'MarkerFaceColor','k')
xlabel('Total IPTG used (\mu M)');
ylabel('Switch time (h)')
set(gca,'XScale','log')

% Bar plot of values of inducer feed-in rates along Pareto:
subplot(3,1,2)
% bar(1:length(IndParams(:,1)),IndParams(idx,1))
plot(ObjVals(:,1),IndParams(:,1),'ko-','MarkerSize',4,'MarkerFaceColor',[0.7,0.7,0.7])
ylabel('IPTG feed-in flux (\mu M/h)')
% xlabel('Pareto solution')
xlabel('Total IPTG used (\mu M)')
set(gca,'XScale','log')

% Bar plot of values of OAxt along Pareto:
subplot(3,1,3)
% bar(1:length(IndParams(:,2)),IndParams(idx,2))
plot(ObjVals(:,1),IndParams(:,2),'ko-','MarkerSize',4,'MarkerFaceColor',[0.7,0.7,0.7])
ylabel('Exposure time (h)')
% xlabel('Pareto solution')
xlabel('Total IPTG used (\mu M)')
set(gca,'XScale','log')

