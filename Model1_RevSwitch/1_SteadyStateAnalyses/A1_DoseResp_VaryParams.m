% ------------------------------------------------------------------------
%           PLOTTING CHANGE IN STEADY STATE DOSE-RESPONSE CURVE
%                        FOR CHANGES IN PARAMETERS
%             Continuous culture of system controlling growth
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 04/11/2020
% ------------------------------------------------------------------------

function A1_DoseResp_VaryParams(arch,FigPlotNum)

% --- Define parameters, which to vary and how ---------------------------

% Nominal (optimally fitted) model parameters:
% params = FAU_MetEng_Params;

% Global parameters to generate dose-response curves:
% ###-###-###-###-###-###-###-###-###-###-###-###-###-###
  if isempty(arch) == 1
      arch = 'NAR';                  % <<<|ADJUST|###########
  end
  if isempty(FigPlotNum) == 1
      FigPlotNum = 1;                % <<<|ADJUST|###########
  end
  params = Model1_Params(arch);
  rangeOA = logspace(-4.5,-1.5,1000); % <<<|ADJUST|###########
  numCurves = 9;                 % <<<|ADJUST|###########
  scaleP_min = 0.1; % 10% of nominal value
  scaleP_max = 10;  % 20x nominal value
  
  % Which params to vary?
%   p_names = {'b_D','b_R','a_D','a_R','K_D','K_R'};
%   p_pos =   [   6 ,   2 ,   7 ,   3 ,   8 ,  4  ];
  p_names = {'b_D','b_R','a_D','a_R','K_D','K_R'};
  p_pos =   [   4 ,   1 ,   5 ,   2 ,   6 ,  3  ];
% ###-###-###-###-###-###-###-###-###-###-###-###-###-###


% --- Generate Dose-Response for Varying Parameters ----------------------

tic
figure(FigPlotNum); clf
w = waitbar(0,'Please wait ...');
for p = 1:length(p_names) % Varying each parameter in turn, from list of p_names.
    % FigPlotNum = p + 1;
    
    ParOI = params(p_pos(p)); % Define the value of the parameter of interest
%     paramRng = [logspace(log10(ParOI*scaleP_min),log10(ParOI*scaleP_max),numCurves),ParOI]; % Nominal param value at end
    paramRng = [logspace(log10(scaleP_min),log10(scaleP_max),numCurves)*ParOI,ParOI];
    
    % Curve colors for plotting:
    colrs = gray(numCurves+1);
    colrs(end,:) = [];
    colrs = colrs(numCurves:-1:1,:);
    
    for i = 1:length(paramRng)
        % Define adjusted set of parameters:
        params_t = params; % Reset to set of nominal (fitted) values
        params_t(p_pos(p)) = paramRng(i); % Define new value
        
        % Generating the dose-response curves:
        VarSSs = m_DoseResp(params_t,arch,rangeOA,[],[]);
        
        % Reorder results (in descending order of ssFadR value) for plots:
        [~,idx] = sort(VarSSs(:,2),'descend');
        VarSSs = VarSSs(idx,:);
        
        % Plotting steady state dose-response of FadR and FadD:
        % figure(FigPlotNum)
        % if i == 1; clf; end
        
        % Finding stable and unstable steady states:
        % stable_p   = find(VarSSs(:,end) == 1);
        unstable_p = find(VarSSs(:,end) == 0);
        
        % Plotting:
        subplot(ceil(length(p_names)/2),2,p)
        if i < length(paramRng)
            col = colrs(i,:);
        else
            col = 'b';
        end
        hold on
%         plot(VarSSs(:,1), VarSSs(:,2),'-o','Color',col,'MarkerFaceColor',col,'Markersize',4);
%         plot(VarSSs(unstable_p,1), VarSSs(unstable_p,2),'ro','MarkerFaceColor','r','Markersize',4);
        plot(VarSSs(:,1), VarSSs(:,2),'-','Color',col,'LineWidth',3);
        plot(VarSSs(unstable_p,1), VarSSs(unstable_p,2),'r-','LineWidth',3);
        hold off
        xlabel('Inducing OA (\mu M)'); ylabel('Free FadR (\mu M)')
        title ([arch,' - Varying ',p_names{p}])
        set(gca,'Xscale','log')
        set(gca,'YScale','lin')
        
        % subplot(1,2,2)
        % if i < length(paramRng)
        %     col = colrs(i,:);
        % else
        %     col = 'b';
        % end
        % hold on
        % plot(VarSSs(:,1), VarSSs(:,3),'-o','Color',col,'MarkerFaceColor',col,'Markersize',4);
        % plot(VarSSs(unstable_p,1), VarSSs(unstable_p,3),'ro','MarkerFaceColor','r','Markersize',4);
        % hold off
        % xlabel('Inducing OA (\mu M)'); ylabel('FadD (\mu M)')
        % set(gca,'Xscale','log')
    end
    xlim([min(rangeOA),max(rangeOA)])
    set(gca,'YScale','log')
    ylim([1e-7 1])
    
    waitbar(p/length(p_names))
end
toc
close(w)
