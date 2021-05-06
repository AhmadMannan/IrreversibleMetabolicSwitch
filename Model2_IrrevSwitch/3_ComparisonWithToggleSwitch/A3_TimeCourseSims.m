% ------------------------------------------------------------------------
% RUN DYNAMICS OF INDUCING
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 22/03/2021
% ------------------------------------------------------------------------

function A3_TimeCourseSims


% Define optimal induction regime:
load('p_irrev')
p_nom = p_irrev;
p_nom(16) = 0; % do not consider fixed inducer conc in media
% p_nom(12) = 0.011;

% Find optimal solution closest to ideal/utopia point:
load('IndParams')
load('ObjVals')
objv = ObjVals.^2;
objv = sqrt(sum(objv,2));
pos = find(objv == min(objv));
ip = IndParams(pos,:);

% Indicate this point in Figure 1 plots:
figure(7)
subplot(3,1,1)
hold on
plot(ObjVals(pos,1),ObjVals(pos,2),'ro','MarkerSize',6,'MarkerFaceColor','r')
hold off
set(gca,'YScale','log')
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
pOpt = p_nom;
pOpt(16) = 0;
pOpt([3,4]) = ip;

% Simulate time course dynamics of this solution:
[t,x] = m_TimeCourseDynSim(pOpt,100);

% Output growth rates:
I_ff = [];
GR = [];
for i = 1:length(t)
    [~,~,~,i_ff,gr] = Model3_TS_Dyn(t(i),x(i,:),pOpt);
    I_ff = [I_ff; i_ff];
    GR = [GR; gr];
end


% Plot time course dynamics:
figure(8); clf

subplot(2,2,1) % plot of inducer feed-in flux over time:
plot(t,I_ff,'k-','LineWidth',2)
xlabel('Time (h)'); ylabel('IPTG feed-in flux (\mu M/h)')
set(gca,'XScale','log')
xlim([4e-2,100])

subplot(2,2,2) % plot time course of inducer
plot(t,x(:,1),'k-','LineWidth',2)
xlabel('Time (h)'); ylabel('IPTG (\mu M)')
set(gca,'XScale','log')
xlim([4e-2,100])

subplot(2,2,3) % plot time course of LacI, TetR and complex (C)
plot(t,x(:,3),'-','LineWidth',2,'Color','b') % FadR
hold on
plot(t,x(:,5),'-','LineWidth',2,'Color','r') % TetR
plot(t,x(:,4),'-','LineWidth',2,'Color',[0.7,0.7,0.7]) % complex
plot([4e-2,100],[0.05,0.05],'k--')
hold off
xlabel('Time (h)'); ylabel('Concentration (\mu M)')
legend('LacI','TetR','Complex','Location','NorthWest')
set(gca,'XScale','log','YScale','log')
xlim([4e-2,100])
% ylim([5e-3,3])

subplot(2,2,4) % plot of growth rate over time
plot(t,GR,'k-','LineWidth',2)
xlabel('Time (h)')
ylabel('Growth rate (1/h)')
set(gca,'XScale','log')
xlim([4e-2,100])
