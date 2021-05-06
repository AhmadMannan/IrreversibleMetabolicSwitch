% ------------------------------------------------------------------------
% Plotting mean and var of eFAST sensitivity indices based on re-sampling
% ------------------------------------------------------------------------
% Author: Ahmad A. Mannan
% Date: 12/12/2018
% ------------------------------------------------------------------------

function Results_BarPlot_SensIdx(s_LV,figNum)


% Defining results needed for calculations:
rangeSi = s_LV.rangeSi;
rangeSti = s_LV.rangeSti;
p_Si = s_LV.p_Si;
p_Sti = s_LV.p_Sti;
[NumModelPars,~,NumModelVars] = size(s_LV.Si);

% Calculating mean and 2*standard deviation of sensitivity indices, for
% each of the model variables:
Si_MnV = [];
Sti_MnV = [];
for i = 1:NumModelVars % # of vars in model
    Si_t = squeeze(rangeSi(:,:,:,i));
    Sti_t = squeeze(rangeSti(:,:,:,i));
    
    % Calculating the mean and var index value over re-samples:
    Si_MnV{i} = [mean(Si_t,2), 2*std(Si_t,0,2)]; %#ok<*AGROW>
    Sti_MnV{i} = [mean(Sti_t,2), 2*std(Sti_t,0,2)];
end

% Bar plot of mean and stdev of 1st- and total-order sensitivity indices:
figure(figNum); clf

for i = 1:NumModelVars
    subplot(1,2,i)
    
    % Plotting mean index and 2*stdev errorbars:
    bar(Sti_MnV{i}(:,1))
    hold on
    bar(Si_MnV{i}(:,1))
    errorbar(Si_MnV{i}(:,1), Si_MnV{i}(:,2),'.')
    errorbar(Sti_MnV{i}(:,1), Sti_MnV{i}(:,2),'.')
    hold off
    legend('total-order','first-order')
    
    % Labelling parameters that give significantly different sensitivity
    % index than index of dummy parameter:
    for j = 1:NumModelPars-1
        if p_Si(j) < 0.05 % Assuming significance level of 0.05
            hold on
            text(j,Si_MnV{i}(j,1)+0.5*Si_MnV{i}(j,2),'*','FontSize',20,'VerticalAlignment','bottom','HorizontalAlignment','center');
            hold off
        end
        if p_Sti(j) < 0.05 % Assuming significance level of 0.05
            hold on
            text(j,Sti_MnV{i}(j,1)+0.5*Sti_MnV{i}(j,2),'*','FontSize',20,'VerticalAlignment','bottom','HorizontalAlignment','center');
            hold off
        end
    end
    
    % Labels:
    set(gca,'XTickLabel',{'\alpha','\beta','\sigma','\delta','dummy'})
    title(['eFAST Sensitivities for var',num2str(i)])
    ylim([0,1.1])
end

