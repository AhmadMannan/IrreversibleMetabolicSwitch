% ------------------------------------------------------------------------
%    GENERATING PLOTS OF HOW TUNING PARAMETERS AFFECTS BISTABLE RANGE,
%   INDUCTION THRESHOLD AND REVERSION THRESHOLD OF DOSE-RESPONSE CURVE
%             Continuous culture of system controlling growth
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 10/11/2020
% ------------------------------------------------------------------------

function A2_Plots_ParamsAffectBistRng %(Yparpos,Xparpos)


%% -----------------------------------------------------------------------
% [I] Varying each parameter in turn for each circuit topology
% ------------------------------------------------------------------------

% Define values for scaling parameter nominal values:
rngSP = logspace(-1,log10(10),30);
rngSP_DR_pos = ceil((length(rngSP)/4))*[1,2,3];
rngSP_DR = rngSP(rngSP_DR_pos); % used for plotting three samples of the dose-response curves

% Define parameters:
p_names = {'a_R','b_D','a_D'};
p_pos = [2,4,5];

arch = {'NAR','PAR'};
for i = 1:2 % arch = NAR or PAR
    % Define parameters:
    p = Model1_Params(arch{i});
    
    % Define output matrix:
    initM = zeros(length(rngSP),length(p_pos));
    BR{i} = initM; % Bistable range
    IT{i} = initM; % Induction threshold
    RT{i} = initM; % Reversion threshold
    
    for j = 1:length(p_pos) % looking at varying particular param
        p_t = p;
        
        b = waitbar(0,['Varying param ',num2str(j),' for circuit with ',arch{i},'...']);
        for k = 1:length(rngSP) % spaning over values for param
            p_t(p_pos(j)) = p(p_pos(j)) * rngSP(k);
            
            % Generate the dose-response curve for FadR vs OA, and
            % calculate the BR, IT and RT:
            [BR_t,IT_t,RT_t,DR_t] = m_QuantDoseResp_BR_IT_RT(p_t,arch{i},[]);
            
            % Save results:
            BR{i}(k,j) = BR_t;
            if BR_t == 0
                IT{i}(k,j) = NaN;
                RT{i}(k,j) = NaN;
            else
                IT{i}(k,j) = IT_t;
                RT{i}(k,j) = RT_t;
            end
            
            % Saving dose-resp from 3 samples of variation in params:
            z = find(rngSP(k) == rngSP_DR);
            if isempty(z) == 0 %#ok<*EFIND>
                DR{i}{j}{z} = DR_t;
            end
            
            waitbar(k/length(rngSP))
        end
        close(b)
    end
end

%%
% --- Plotting -----------------------------------------------------------

% Plots of the effect of varying parameters on BR, IT and RT for circuit
% with NAR (red) and PAR (blue) ... 
figure(1); clf

% ... varying bD effect on bistable range:
subplot(3,3,1)
plot(rngSP,BR{1}(:,2),'r-','LineWidth',2) % circuit with NAR
hold on
plot(rngSP,BR{2}(:,2),'b-','LineWidth',2) % circuit with PAR
% marking the three sampled param values used to plot dose-response curves:
plot(rngSP(rngSP_DR_pos),BR{1}(rngSP_DR_pos,2),'rx','MarkerSize',8)
plot(rngSP(rngSP_DR_pos),BR{2}(rngSP_DR_pos,2),'bx','MarkerSize',8)
hold off
xlabel('Scale of FadD expression leakiness b_D')
ylabel('Bistable range (\mu M)')
set(gca,'YScale','log','XScale','log')
ylim([3e-5,7e-2])
% ... varying bD effect on induction threshold:
subplot(3,3,2)
plot(rngSP,IT{1}(:,2),'r-','LineWidth',2) % circuit with NAR
hold on
plot(rngSP,IT{2}(:,2),'b-','LineWidth',2) % circuit with PAR
% marking the three sampled param values used to plot dose-response curves:
plot(rngSP(rngSP_DR_pos),IT{1}(rngSP_DR_pos,2),'rx','MarkerSize',8)
plot(rngSP(rngSP_DR_pos),IT{2}(rngSP_DR_pos,2),'bx','MarkerSize',8)
hold off
xlabel('Scale of FadD expression leakiness b_D')
ylabel('Induction threshold (\mu M)')
set(gca,'YScale','log','XScale','log')
ylim([3e-5,7e-2])
% ... varying aD effect on reversion threshold:
subplot(3,3,3)
plot(rngSP,RT{1}(:,3),'r-','LineWidth',2) % circuit with NAR
hold on
plot(rngSP,RT{2}(:,3),'b-','LineWidth',2) % circuit with PAR
% marking the three sampled param values used to plot dose-response curves:
plot(rngSP(rngSP_DR_pos),RT{1}(rngSP_DR_pos,3),'rx','MarkerSize',8)
plot(rngSP(rngSP_DR_pos),RT{2}(rngSP_DR_pos,3),'bx','MarkerSize',8)
hold off
xlabel('Scale of FadD promoter strength a_D')
ylabel('Reversion threshold (\mu M)')
set(gca,'YScale','log','XScale','log')
ylim([3e-5,7e-2])


% Plots of dose-response curve from three samples of the parameters:
for p = 1:3
    if p == 1 || p == 2
        par = 2;
    elseif p == 3
        par = 3;
    end
    subplot(3,3,p+3) % ... for circuit with NAR
    for i = 1:length(rngSP_DR)
        drCurve = DR{1}{par}{i};
        unstPart = find(drCurve(:,3)==0);
        hold on
        plot(drCurve(:,1),drCurve(:,2),'k-','LineWidth',2)
        plot(drCurve(unstPart,1),drCurve(unstPart,2),'r-','LineWidth',2)
        hold off
    end
    set(gca,'XScale','log','YScale','log'); xlim([1e-4,1e-2])
    xlabel('Oleic acid (\mu M)'); ylabel('FadR (\mu M)'); title(['NAR - vary ',p_names{par}])
    subplot(3,3,p+6) % ... for circuit with PAR
    for i = 1:length(rngSP_DR)
        drCurve = DR{2}{par}{i};
        unstPart = find(drCurve(:,3)==0);
        hold on
        plot(drCurve(:,1),drCurve(:,2),'k-','LineWidth',2)
        plot(drCurve(unstPart,1),drCurve(unstPart,2),'r-','LineWidth',2)
        hold off
    end
    set(gca,'XScale','log','YScale','log'); xlim([1e-4,1e-2])
    xlabel('Oleic acid (\mu M)'); ylabel('FadR (\mu M)'); title(['PAR - vary ',p_names{par}])
end

%%

%% -----------------------------------------------------------------------
% [II] Varying sT for each circuit topology
% ------------------------------------------------------------------------

% Define range of values of sT:
sT_val = linspace(0,0.95,30);
sTval_DR_pos = ceil((length(sT_val)/4))*[1,2,3];
sTval_DR = sT_val(sTval_DR_pos); % used for plotting three samples of the dose-response curves

% Define parameters:
p_names = {'s_T'};
p_pos = 17;

arch = {'NAR','PAR'};
for i = 1:2 % arch = NAR or PAR
    % Define parameters:
    p = Model1_Params(arch{i});
    
    % Define output matrix:
    initM = zeros(length(sT_val),length(p_pos));
    BR{i} = initM; % Bistable range
    IT{i} = initM; % Induction threshold
    RT{i} = initM; % Reversion threshold
    
    for j = 1:length(p_pos) % looking at varying particular param
        p_t = p;
        
        b = waitbar(0,['Varying param ',p_names{j},' for circuit with ',arch{i},'...']);
        for k = 1:length(sT_val) % spaning over values for sT
            p_t(p_pos(j)) = sT_val(k);
            
            % Generate the dose-response curve for FadR vs OA, and
            % calculate the BR, IT and RT:
            [BR_t,IT_t,RT_t,DR_t] = m_QuantDoseResp_BR_IT_RT(p_t,arch{i},[]);
            
            % Save results:
            BR{i}(k,j) = BR_t;
            if BR_t == 0
                IT{i}(k,j) = NaN;
                RT{i}(k,j) = NaN;
            else
                IT{i}(k,j) = IT_t;
                RT{i}(k,j) = RT_t;
            end
            
            % Saving dose-resp from 3 samples of variation in params:
            z = find(sT_val(k) == sTval_DR);
            if isempty(z) == 0 %#ok<*EFIND>
                DR{i}{j}{z} = DR_t;
            end
            
            waitbar(k/length(sT_val))
        end
        close(b)
    end
end

%%
% Plots of the effect of varying sT on BR, IT and RT for circuit with NAR
% (red) and PAR (blue) ...

% Plot of effect of varying sT on bistable range, for both circuits:
figure(2); clf
subplot(2,2,1)
plot(sT_val,IT{1},'r-','LineWidth',2) % circuit with NAR
hold on
plot(sT_val,IT{2},'b-','LineWidth',2) % circuit with PAR
% marking the three sampled values of sT used to plot the dose-resp curves:
plot(sT_val(sTval_DR_pos),IT{1}(sTval_DR_pos),'rx','MarkerSize',8)
plot(sT_val(sTval_DR_pos),IT{2}(sTval_DR_pos),'bx','MarkerSize',8)
hold off
xlabel('Severity of growth attenuation s_T')
ylabel('Induction threshold (\mu M)')
set(gca,'YScale','log','XScale','lin')
% ylim([3e-5,7e-2])

subplot(2,2,2)
plot(sT_val,RT{1},'r-','LineWidth',2) % circuit with NAR
hold on
plot(sT_val,RT{2},'b-','LineWidth',2) % circuit with PAR
% marking the three sampled values of sT used to plot the dose-resp curves:
plot(sT_val(sTval_DR_pos),RT{1}(sTval_DR_pos),'rx','MarkerSize',8)
plot(sT_val(sTval_DR_pos),RT{2}(sTval_DR_pos),'bx','MarkerSize',8)
hold off
xlabel('Severity of growth attenuation s_T')
ylabel('Reversion threshold (\mu M)')
set(gca,'YScale','log','XScale','lin')
% ylim([3e-5,7e-2])

% Plots of dose-response curve from three samples of the value of sT:
subplot(2,2,3) % ... for circuit with NAR
for i = 1:length(sTval_DR)
    drCurve = DR{1}{1}{i};
    unstPart = find(drCurve(:,3)==0);
    hold on
    plot(drCurve(:,1),drCurve(:,2),'k-','LineWidth',2)
    plot(drCurve(unstPart,1),drCurve(unstPart,2),'r-','LineWidth',2)
    hold off
end
set(gca,'XScale','log','YScale','log'); xlim([1e-4,1e-1])
xlabel('Oleic acid (\mu M)'); ylabel('FadR (\mu M)'); title('NAR - vary s_T')

subplot(2,2,4) % ... for circuit with PAR
for i = 1:length(sTval_DR)
    drCurve = DR{2}{1}{i};
    unstPart = find(drCurve(:,3)==0);
    hold on
    plot(drCurve(:,1),drCurve(:,2),'k-','LineWidth',2)
    plot(drCurve(unstPart,1),drCurve(unstPart,2),'r-','LineWidth',2)
    hold off
end
set(gca,'XScale','log','YScale','log'); xlim([1e-4,1e-1])
xlabel('Oleic acid (\mu M)'); ylabel('FadR (\mu M)'); title('PAR - vary s_T')

