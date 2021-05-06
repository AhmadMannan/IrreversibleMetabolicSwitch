p0 = Model2_Dyn_Params; % parameters of irreversible switch
comp = 'COMP';

% Run multiobjective optimization:
hoursCult = 100; % length of induction/production process
OAxt = 7;
OAin = logspace(0.1,0.5,20);
col = gray(length(OAin) + 1);
col(end,:) = [];
col = col(end:-1:1,:);

figure(1) %;clf
w = waitbar(0,'Please wait ...');
for i = 1:length(OAin)
    [Obj,t,x] = m_MultiObjs([OAin(i),OAxt],p0,comp,hoursCult);
    
    % plot time course:
    subplot(1,2,1)
    ylim([1e-6,100])
    plot(t,x(:,1),'r-','LineWidth',2)
    hold on
    plot(t,x(:,7),'k-','LineWidth',2)
    plot([0,hoursCult],[0.05,0.05],'k--')
    hold off
    set(gca,'YScale','lin')
    
    % plot objectives
    subplot(1,2,2)
    hold on
    plot(Obj(1),Obj(2),'o','MarkerSize',6,'Color',col(i,:),'MarkerFaceColor',col(i,:))
    hold off
    
    pause(0.3)
    waitbar(i/length(OAin))
end
close(w)
xlabel('Obj1 = Total OA used (\mu M)')
ylabel('Obj2 = Switch time (h)')
