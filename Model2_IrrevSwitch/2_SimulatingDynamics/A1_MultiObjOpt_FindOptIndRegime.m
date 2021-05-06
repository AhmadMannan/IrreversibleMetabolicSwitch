% ------------------------------------------------------------------------
% PERFORMING MULTIOBJECTIVE OPTIMIZATION TO FIND OPTIMAL INDUCTION REGIMES
% THAT MINIMIZE (I) TOTAL INDUCER USED, AND (II) SWITCH TIME
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 24/11/2020
% ------------------------------------------------------------------------

function A1_MultiObjOpt_FindOptIndRegime


%%
% --- Perform multiobjective optimization --------------------------------

% Define control system:
p0 = Model2_Dyn_Params; % parameters of irreversible switch
comp = 'COMP';

% Run multiobjective optimization ...
hoursCult = 100; % length of induction/production process
lb = [1, 0.1]; % for parameters [OAin,OAxt] (feed-in flux and exposure time)
ub = [50, 5];
options = optimoptions(@gamultiobj,'Display','iter','ParetoFraction',0.6,'PlotFcn','gaplotpareto','PopulationSize',100,'UseParallel',true);

% ... and generating the Pareto front:
tic
[IndParams,ObjVals,exitflag] = gamultiobj(@(p)m_MultiObjs(p,p0,comp,hoursCult),length(lb),[],[],[],[],lb,ub,[],options);
toc

disp(exitflag)


%% --- Plotting -----------------------------------------------------------

% Delete solutions with switch time above an acceptable time limit (10hrs):
del = find(ObjVals(:,2) > 24);
ObjVals(del,:) = [];
IndParams(del,:) = [];

% Save solutions:
% save('IndParams','IndParams')
% save('ObjVals','ObjVals')

% Sort Pareto solutions in descending order of Obj1 (total OA used):
[~,idx] = sort(ObjVals(:,1),'descend');


figure(1); clf

% Plotting Pareto front of solutions:
subplot(3,1,1)
plot(ObjVals(idx,1),ObjVals(idx,2),'ko-','MarkerSize',4,'MarkerFaceColor','k')
xlabel('Total OA used (\mu M)');
ylabel('Switch time (h)')
set(gca,'XScale','log')

% Bar plot of values of OAin along Pareto:
subplot(3,1,2)
% bar(1:length(IndParams(:,1)),IndParams(idx,1))
plot(ObjVals(idx,1),IndParams(idx,1),'ko-','MarkerSize',4,'MarkerFaceColor',[0.7,0.7,0.7])
ylabel('OA feed-in flux (\mu M/h)')
% xlabel('Pareto solution')
xlabel('Total OA used (\mu M)')
set(gca,'XScale','log')

% Bar plot of values of OAxt along Pareto:
subplot(3,1,3)
% bar(1:length(IndParams(:,2)),IndParams(idx,2))
plot(ObjVals(idx,1),IndParams(idx,2),'ko-','MarkerSize',4,'MarkerFaceColor',[0.7,0.7,0.7])
ylabel('Exposure time (h)')
% xlabel('Pareto solution')
xlabel('Total OA used (\mu M)')
set(gca,'XScale','log')

