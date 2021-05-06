% ------------------------------------------------------------------------
% TIME COURSE SIMULATION OF THE OPTIMAL INDUCTION REGIME
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 24/11/2020
% ------------------------------------------------------------------------

function A2_TimeCourseSim

%%
% Define optimal induction regime:
load('IndParams')
load('ObjVals')
% optPos = find(min(ObjVals(:,1)) == ObjVals);
pos = 5;
ip = IndParams(pos,:);

% Indicate this point in Figure 1 plots:
figure(1)
subplot(3,1,1)
hold on
plot(ObjVals(pos,1),ObjVals(pos,2),'ro','MarkerSize',6,'MarkerFaceColor','r')
hold off
subplot(3,1,2)
hold on
plot(ObjVals(pos,1),IndParams(pos,1),'ro','MarkerSize',6,'MarkerFaceColor','r')
hold off
subplot(3,1,3)
hold on
plot(ObjVals(pos,1),IndParams(pos,2),'ro','MarkerSize',6,'MarkerFaceColor','r')
hold off
%%
% Define optimal parameters:
pOpt = Model2_Dyn_Params;
pOpt([23,24]) = ip;

% Simulate time course dynamics of this solution:
[t,x] = m_TimeCourseDynSim(pOpt,'COMP',100);

%
% Plot time course dynamics:
figure(3); clf

subplot(2,1,1) % plot time course of OA
plot(t,x(:,7),'k-','LineWidth',2)
xlabel('Time (h)'); ylabel('Oleic acid (\mu M)')
set(gca,'XScale','log')
xlim([4e-2,100])

subplot(2,1,2) % plot time course of FadR, TetR and complex (C)
plot(t,x(:,1),'-','LineWidth',2,'Color','b') % FadR
hold on
plot(t,x(:,6),'-','LineWidth',2,'Color','r') % TetR
plot(t,x(:,4),'-','LineWidth',2,'Color',[0.7,0.7,0.7]) % complex
plot([4e-2,100],[0.05,0.05],'k--')
hold off
xlabel('Time (h)'); ylabel('Concentration (\mu M)')
legend('FadR','TetR','Complex','Location','NorthWest')
set(gca,'XScale','log','YScale','log')
xlim([4e-2,100])
ylim([5e-3,3])
