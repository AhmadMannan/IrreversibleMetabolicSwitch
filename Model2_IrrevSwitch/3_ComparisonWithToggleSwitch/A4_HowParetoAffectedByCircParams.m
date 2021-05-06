% ------------------------------------------------------------------------
% STUDY OF HOW PARETO FRONT IS AFFECTED BY TUNING CONTROL SYSTEM PARAMS.
%  - Obj_1 = total OA used (uM)
%  - Obj_2 = switch time (h)
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 23/03/2021
% ------------------------------------------------------------------------

function A4_HowParetoAffectedByCircParams


%%
% Loading induction parameters of the Pareto front:
load('IndParams') %#ok<*LOAD>
ip = IndParams; %#ok<*NODEF>

% Load parameters of irreversible switch:
load('p_irrev')
p0 = p_irrev;
p0(16) = 0;

% Considering only tuning in the experimentally-accessible parameters:
p_Names = {'b_L','a_L','K_L','b_T','a_T','K_T'};
p_Pos = [8,9,10,11,12,13];

% Define range of tuning each parameter:
rngTuning = linspace(0.9,1.1,3);

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
            obj_t = m_MultiObjs(ip(k,:),p_t,100);
            Objs{i}(k,[2*j-1,2*j]) = obj_t;
        end
        
        waitbar(j/length(rngTuning))
    end
    close(w)
end

%%
% Plotting:
figure(9); clf

col = gray(length(rngTuning)+1);
col(end,:) = [];
col = col(end:-1:1,:);

numSP = sqrt(length(p_Pos)); % number of subplots

for i = 1:length(p_Pos)
    subplot(round(numSP),ceil(numSP),i)
    
    [~,idx] = sort(Objs{i}(:,3),'descend');
    
    for j = 1:length(rngTuning)
        if j == 2
            hold on
            plot(Objs{i}(idx,2*j-1),Objs{i}(idx,2*j),'o-','MarkerSize',6,'Color','b')
            hold off
        else
            hold on
            plot(Objs{i}(idx,2*j-1),Objs{i}(idx,2*j),'o-','MarkerSize',6,'Color',col(j,:))
            hold off
        end
    end
    xlabel('Obj_1'); ylabel('Obj_2'); title(p_Names{i})
    set(gca,'XScale','log','YScale','lin')
    xlim([11,17])
end
