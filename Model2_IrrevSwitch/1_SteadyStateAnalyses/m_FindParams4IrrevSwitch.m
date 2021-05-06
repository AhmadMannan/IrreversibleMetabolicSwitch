% ------------------------------------------------------------------------
% DETERMINE PARAMETER REGIME WHERE NUMBER OF STEADY STATES = 3 WHEN OA = 0
% where the switch is irreversible.
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 19/11/2020
% ------------------------------------------------------------------------

function [P3SS,rng] = m_FindParams4IrrevSwitch(comp,fig)


%%
% ------------------------------------------------------------------------
% [i] Explore over 4 parameters to find parameter regime where switch is
%     irreversible. Do for system with PAR with COMP and NONCOMP binding
%     and system w/o PAR (NOPAR).
% ------------------------------------------------------------------------
% comp = 'COMP';

% Define nominal parameter values:
params = Model2_Params;

% Define range of scaling to params:
numVals = 71;
numPs = 4;
pNames = {'a_R','a_T','K_{Ri}','K_T'};
rng{1} = logspace(-4,3,numVals); % scale a_R
rng{2} = logspace(-4,3,numVals); % scale a_T
rng{3} = logspace(-4,3,numVals); % scale K_Ri
rng{4} = logspace(-4,3,numVals); % scale K_T
pscale = zeros(numVals^numPs,numPs);
count = 0;
for i = 1:numVals
    for j = 1:numVals
        for k = 1:numVals
            for l = 1:numVals
                count = count + 1;
                pscale(count,:) = [rng{1}(i),rng{2}(j),rng{3}(k),rng{4}(l)];
            end
        end
    end
end

NumSS = zeros(numVals^numPs,numPs+1);
parfor i = 1:length(pscale(:,1))
    % Tracker:
    if mod(i,10000) == 0
        disp(i)
    end
    
    psc = ones(size(params));
    % Params [a_R  a_T  K_Ri  K_T]:
    psc(     [2,   19,  20,   21 ]) = pscale(i,:);
    p_t = params .* psc;
    
    % Calculate the number of steady states for param values:
    numSS = m_NumRssAt0OA(p_t,comp);
    
    % Save param scaling and num of FadR steady states:
    if numSS > 1
        NumSS(i,:) = [pscale(i,:), numSS]; %#ok<*AGROW>
    end
end

% Find parameter scalings that achieve irreversible bistable switch with 3
% steady states at OA = 0:
P3SS = NumSS(find(NumSS(:,end) > 1),1:end-1); %#ok<*FNDSB>

% Select a parameter set where switch was found to behave irreversibly:
SelParTuning = P3SS(round(length(P3SS(:,1))/3),:);
psc = ones(size(params)); % parameter scaling
psc([2,19,20,21]) = SelParTuning;
p_IS = params .* psc; % params for irreversible switch

% Creating an image of the design space and the portion where switch is
% irreversible:
figure(fig); clf
colormap(gray)
sppos = zeros(numPs); sppos(1:numPs^2) = 1:numPs^2; sppos = sppos'; % subplot positions
for i = 1:numPs-1
    for j = i+1:numPs
        % Defining parameters that give irreversibility:
        m_t = ones(numVals);
        for k = 1:length(P3SS(:,1))
            m_t(find(P3SS(k,i) == rng{i}), find(P3SS(k,j) == rng{j})) = 0;
        end
        
        % Define subplots:
        subplot(numPs,numPs,sppos(i,j))
        imagesc(m_t,[-2,1])
        set(gca,'YDir','normal')
        xlabel(pNames{j}); ylabel(pNames{i})
        
        % Add point in plot of system with nominal params:
        hold on
        plot(find(1 == rng{j}),find(1 == rng{i}),'ws','MarkerSize',6,'MarkerFaceColor',[0.2,0.2,0.2]);
        hold off
        
%         % Add point in plot of selected param for irreversible switch:
%         hold on
%         plot(find(SelParTuning(j) == rng{j}),find(SelParTuning(i) == rng{i}),'rs','MarkerSize',6,'MarkerFaceColor','r');
%         hold off
        
        if i == 1 && j == 2
            title(comp)
        end
    end
end

