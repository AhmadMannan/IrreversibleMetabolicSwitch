% ------------------------------------------------------------------------
% PLOT SENSITIVITY INDICES AND HIGHLIGHT THOSE SIGNIF. DIFF. FROM DUMMY
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% ------------------------------------------------------------------------

% figure(1);clf

% Preparing results to plot average and 1 standard deviation from search
% repeats:
% ... first order sensitivities:
Si = squeeze(s_KM.rangeSi);
Si_avg = mean(Si,2);
Si_std = std(Si,0,2);
% ... total order sensitivities:
% Sti = squeeze(s_KM.rangeSti);
% Sti_avg = mean(Sti,2);
% Sti_std = std(Sti,0,2);

% Plot:
% bar(1:length(Sti_avg),Sti_avg,'FaceColor',[0.8,0.8,0.8])
bar(1:length(Si_avg),Si_avg,'FaceColor',[0.4,0.4,0.4])
hold on
% errorbar(1:length(Sti_avg),Sti_avg,Sti_std,'k.')
for i = 1:length(Si(1,:))
    plot(1:length(Si_avg),Si(:,i),'ko','MarkerSize',6,'LineWidth',1)
end
errorbar(1:length(Si_avg),Si_avg,Si_std,'k.')
hold off
set(gca,'XTickLabel',efast_var,'YScale','log')
xlabel('Kinetic model parameters')
ylabel('Measured sensitivity')
% legend('Total order S_{Ti}','First order S_{i}','Location','SouthEast')
legend('First order S_{i}','Location','SouthEast')
ylim([0 1])

% Indentifying parameters with sensitivities significantly different from
% dummy parameter:
p_sigDiff = find(s_KM.p_Si < 0.01);
hold on
text(p_sigDiff,Si_avg(p_sigDiff) + Si_std(p_sigDiff),'*','FontSize',20,'VerticalAlignment','bottom','HorizontalAlignment','center');
hold off

% title('GSA with eFAST')