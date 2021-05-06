% ------------------------------------------------------------------------
% STUDY HOW CONTROL SYSTEM PARAMETERS AFFECT INDUCTION REGIME
% (THE PARETO FRONT)
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 12/11/2020
% ------------------------------------------------------------------------


function A2_ParamsEffectOnObj1VsObj2Curve(p_pos,fig)

%%
% #################
%   fig = 3;
%   p_pos = 2; % Which parameter are we tuning?
% #################

% Names of parameters:
p_names = {'b_R','a_R','K_R','b_D','a_D','K_D','kcatD','Km_D','kcatB','Km_B','B','kf','kr','a_g','K_g'};

% Define range of scaling to nominal parameter values:
% rng = logspace(-0.5,0.5,3);
rng = logspace(-2,0,5);

arch = {'NAR','PAR'};
figure(fig); clf
for a = 1:2
    % Defining set of parameters:
    params = Model1_Params(arch{a});
    params = params';
    
    % Defining coloring for parameters:
    if strmatch(arch{a},'NAR') == 1 %#ok<*MATCH2>
        colr = autumn(length(rng));
        colr = colr(end:-1:1,:);
    elseif strmatch(arch{a},'PAR') == 1
        colr = winter(length(rng));
        colr = colr(end:-1:1,:);
    end
    
    % Generating and plotting the Pareto front for each parameter set:
    w = waitbar(0,['Analysis for ',arch{a},', param ',p_names{p_pos},'...']);
    for i = 1:length(rng)
        % Define parameter change:
        p_t = params;
        p_t(p_pos) = params(p_pos) * rng(i);
        
        % Generate Pareto front:
        rng_OAin = logspace(log10(0.25),log10(5),35); % ... range of OA feed-in flux
        SigOS = 0.3; % Define fraction of 1/Kp for which induction feed-in should stop
        [TotOA_t,SwitchTime_t] = A1_VarOAFeedin(arch{a},p_t,rng_OAin,SigOS,0,0);
        
        % Plotting Pareto front:
        figure(fig)
        % subplot(1,2,a)
        % if i == 1
        %     cla
        % end
        hold on
        plot(TotOA_t,SwitchTime_t,'-','Color',colr(i,:),'LineWidth',3)
        hold off
        xlabel('Total oleic acid used (\mu M)')
        ylabel('Time for FadR to reach 1/K_p (h)')
        title([arch{a},' - Param ',p_names(p_pos)])
        
        waitbar(i/length(rng))
    end
    close(w)
end
