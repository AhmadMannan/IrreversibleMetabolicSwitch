% ------------------------------------------------------------------------
% STUDY OF HOW PARETO FRONT IS AFFECTED BY TUNING CONTROL SYSTEM PARAMS.
%  - Obj_1 = total OA used (uM)
%  - Obj_2 = switch time (h)
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 25/11/2020
% ------------------------------------------------------------------------

function A4_HowParetoAffectedByTuning_bD


%%
% Loading induction parameters of the Pareto front:
load('IndParams') %#ok<*LOAD>
ip = IndParams; %#ok<*NODEF>

% Load parameters of irreversible switch:
load('Params_Eng_IrrevSw')
p0 = Params_Eng_IrrevSw;

% Considering only tuning in the experimentally-accessible parameters:
% p_Names = {'b_R','a_R','K_R','b_D','a_D','K_D','b_T','a_T','K_{Ri}','K_T'};
% p_Pos = [1,2,3,4,5,6,18,19,20,21];
p_Names = {'b_D'};
p_Pos = 4;

% Define range of tuning each parameter:
% rngTuning = linspace(0.75,1.25,3);
rngTuning = linspace(1,10,10);

% Calculating Obj1 and Obj2 for each induction regime along Pareto for
% each parameter tuning:
Objs = [];
for i = 1:length(p_Pos)
    w = waitbar(0,['Tuning param ',p_Names{i}]);
    
    Objs{i} = zeros(length(ip(:,1)),2*length(rngTuning)); %#ok<*AGROW>
    for j = 1:length(rngTuning)
        % Define the system parameters:
        sc = ones(size(p0));
        sc(p_Pos(i)) = rngTuning(j);
        p_t = p0 .* sc;
        
        % Recalculate objectives for each induction regime:
        for k = 1:length(ip(:,1))
            % Define induction regime (OAin and OAxt)
            % p_t2 = p_t;
            % p_t2([23,24]) = IndParams(k,:);
            
            % Calculate and save Obj1 and Obj2:
            obj_t = m_MultiObjs(ip(k,:),p_t,'COMP',100);
            Objs{i}(k,[2*j-1,2*j]) = obj_t;
        end
        
        waitbar(j/length(rngTuning))
    end
    close(w)
end

%
% Plotting:
figure(3); clf
col = gray(length(rngTuning)+1);
col(end,:) = [];
col = col(end:-1:1,:);

numSP = sqrt(length(p_Pos)); % number of subplots

for i = 1:length(p_Pos)
    subplot(ceil(numSP),round(numSP),i)
    
    [~,idx] = sort(Objs{i}(:,3),'descend');
    
    for j = 1:length(rngTuning)
        hold on
        plot(Objs{i}(idx,2*j-1),Objs{i}(idx,2*j),'o-','MarkerSize',6,'Color',col(j,:),'MarkerFaceColor',col(j,:))
        hold off
    end
    xlabel('Obj_1'); ylabel('Obj_2'); title(p_Names{i})
%     if isempty(find(i == [1,5,6,8,9])) == 0
%         ylim([4.14,4.25])
%     elseif isempty(find(i == [3,7,10])) == 0
%         ylim([3.5,5])
%     else
%         ylim([3,6])
%     end
    set(gca,'XScale','log','YScale','lin')
end
ylim([0,5])
