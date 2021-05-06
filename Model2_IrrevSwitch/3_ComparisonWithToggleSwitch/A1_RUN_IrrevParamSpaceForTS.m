% ------------------------------------------------------------------------
% SHOW PARAMETER SPACE WHERE TOGGLE SWITCH IS BISTABLE
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 19/03/2021
% ------------------------------------------------------------------------

function p_irrev = A1_RUN_IrrevParamSpaceForTS


%%
% --- [A] Bistable design space ------------------------------------------

% Define nominal parameter values:
params = Model3_TS_Dyn_Params;

% Define range of scaling to params:
numVals = 15;
numPs = 4;
pNames = {'a_L','a_T','K_L','K_T'};
rng{1} = logspace(-3,3,numVals); % scale a_L
rng{2} = logspace(-3,3,numVals); % scale a_T
rng{3} = logspace(-3,3,numVals); % scale K_L
rng{4} = logspace(-3,3,numVals); % scale K_T
totparamsets = length(rng{1})*length(rng{2})*length(rng{3})*length(rng{4});
pscale = zeros(totparamsets,numPs);
count = 0;
for i = 1:length(rng{1})
    for j = 1:length(rng{2})
        for k = 1:length(rng{3})
            for l = 1:length(rng{4})
                count = count + 1;
                pscale(count,:) = [rng{1}(i),rng{2}(j),rng{3}(k),rng{4}(l)];
            end
        end
    end
end

BistID = zeros(numVals^numPs,numPs+1);
startup_STB % startup Sundials ODE solver for calculation of dose response curves
tic
parfor i = 1:length(pscale(:,1))
    % Tracker:
    if mod(i,10000) == 0
        disp(i)
    end
    
    psc = ones(size(params));
    % Params [a_L  a_T  K_L  K_T]:
    psc(     [9,   12,  10,  13 ]) = pscale(i,:);
    p_t = params .* psc;
    
    % Determine if the system is monostable (1), bistable (2) or
    % irreversible (3) for a given parameter set:
    hystID = m_Det_MonoBistIrrev(p_t);
    
    % Save param scaling and whether or not system is bistable/irrev.:
%     if hystID > 1
        BistID(i,:) = [pscale(i,:), hystID]; %#ok<*AGROW>
%     end
end
toc


%
% Image of projected shadows of parameter space where switch is monostable:
figure(3); clf
colormap(gray)
sppos = zeros(numPs); sppos(1:numPs^2) = 1:numPs^2; sppos = sppos'; % subplot positions
for i = 1:numPs-1
    for j = i+1:numPs
        % Defining parameters for reversible bistability:
        m_t = ones(numVals);
        for k = 1:length(BistID(:,1))
            if BistID(k,end) == 1
                m_t(find(BistID(k,i) == rng{i}), find(BistID(k,j) == rng{j})) = 0; %#ok<*FNDSB> % monostable
            end
        end
        
        % Define subplots:
        subplot(numPs,numPs,sppos(i,j))
        imagesc(m_t,[-2,1])
        set(gca,'YDir','normal')
        xlabel(pNames{j}); ylabel(pNames{i})
        
        if i == 1 && j == 2
            title('Param space for monostability')
        end
    end
end

% Image of projected shadows of parameter space where switch is bistable:
figure(4); clf
colormap(gray)
sppos = zeros(numPs); sppos(1:numPs^2) = 1:numPs^2; sppos = sppos'; % subplot positions
for i = 1:numPs-1
    for j = i+1:numPs
        % Defining parameters for reversible bistability:
        m_t = ones(numVals);
        for k = 1:length(BistID(:,1))
            if BistID(k,end) == 2
                m_t(find(BistID(k,i) == rng{i}), find(BistID(k,j) == rng{j})) = 0; % reversible bistable
            end
        end
        
        % Define subplots:
        subplot(numPs,numPs,sppos(i,j))
        imagesc(m_t,[-2,1])
        set(gca,'YDir','normal')
        xlabel(pNames{j}); ylabel(pNames{i})
        
        if i == 1 && j == 2
            title('Param space for bistability')
        end
    end
end

% Image of projected shadows of parameter space where switch is irrev.:
figure(5); clf
colormap(gray)
sppos = zeros(numPs); sppos(1:numPs^2) = 1:numPs^2; sppos = sppos'; % subplot positions
for i = 1:numPs-1
    for j = i+1:numPs
        % Defining parameters for reversible bistability:
        m_t = ones(numVals);
        for k = 1:length(BistID(:,1))
            if BistID(k,end) == 3
                m_t(find(BistID(k,i) == rng{i}), find(BistID(k,j) == rng{j})) = 0; % reversible bistable
            end
        end
        
        % Define subplots:
        subplot(numPs,numPs,sppos(i,j))
        imagesc(m_t,[-2,1])
        set(gca,'YDir','normal')
        xlabel(pNames{j}); ylabel(pNames{i})
        
        if i == 1 && j == 2
            title('Param space for irreversibility')
        end
    end
end


%

% ------------------------------------------------------------------------
% Selecting param regime for irrevesible TS, closest dist to nominal vals
% ------------------------------------------------------------------------

IrrDes = BistID(:,1:end-1);
dist2nom = (IrrDes - ones(size(IrrDes))).^2;
% dist2nom = IrrDes;
% dist2nom = dist2nom - ones(size(IrrDes));
% dist2nom = dist2nom.^2;
dist2nom = sqrt(sum(dist2nom,2));
dist2nom(find(BistID(:,end) < 3),:) = max(dist2nom);
minDistPos = find(dist2nom == min(dist2nom));
disp(minDistPos)

% Define scaling that achieves irreversible TS:
pscale = IrrDes(minDistPos(1),:);


% ------------------------------------------------------------------------
% Plotting param regime on param space plots
% ------------------------------------------------------------------------

% Plotting dose-response of irrev. TS with params:
params = Model3_TS_Dyn_Params; % set of nominal param values
params(16) = 1;
p_irrev = params;
p_irrev([9,12,10,13]) = params([9,12,10,13])' .* pscale;

% #############################
  save('p_irrev','p_irrev')
% #############################


% Plotting our selected parameter regime for irreversible switch on figs:
numPs = length(pscale(1,:));
sppos = zeros(numPs); sppos(1:numPs^2) = 1:numPs^2; sppos = sppos'; % subplot positions
figure(5)
for i = 1:numPs-1
    for j = i+1:numPs
        % Go to correct subplot:
        subplot(numPs,numPs,sppos(i,j))
        hold on
        plot(find(pscale(j) == rng{j}),find(pscale(i) == rng{i}),'ks','MarkerSize',6,'MarkerFaceColor','g');
        hold off
    end
end


%% ------------------------------------------------------------------------
% Plotting dose-response curve for selected parameter regime
% ------------------------------------------------------------------------

% Generate dose-response:
rngInd = logspace(-4,3,100);
[ss_inc,ss_dec] = m_DoseResp(p_irrev,rngInd);

figure(6); clf

subplot(2,1,1)
plot(ss_inc(:,1),ss_inc(:,1+3),'o-','LineWidth',2,'Color','k')
hold on
plot(ss_dec(:,1),ss_dec(:,1+3),'o-','LineWidth',2,'Color',[0.7,0.7,0.7])
hold off
xlabel('IPTG (\mu M)'); ylabel('LacI (\mu M)')
set(gca,'XScale','log','YScale','log')
legend('increasing doses','decreasing doses','Location','NorthEast')
axis([1e-4,1e2,1e-5,1e1])

subplot(2,1,2)
plot(ss_inc(:,1),ss_inc(:,1+5),'o-','LineWidth',2,'Color','k')
hold on
plot(ss_dec(:,1),ss_dec(:,1+5),'o-','LineWidth',2,'Color',[0.7,0.7,0.7])
hold off
xlabel('IPTG (\mu M)'); ylabel('TetR (\mu M)')
set(gca,'XScale','log','YScale','log')
legend('increasing doses','decreasing doses','Location','SouthEast')
