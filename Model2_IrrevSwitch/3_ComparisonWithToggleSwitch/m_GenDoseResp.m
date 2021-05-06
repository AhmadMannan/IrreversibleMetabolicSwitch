params = Model3_TS_Dyn_Params;
params(16) = 1;
rngInd = logspace(-3,3,30);
% rngInd = 100;
tic
[dr_inc,dr_dec] = m_DoseResp(params,rngInd);
toc

figure(1); clf

plot(dr_inc(:,1),dr_inc(:,1+3)/max(dr_inc(:,1+3)),'k-o')
hold on
plot(dr_dec(:,1),dr_dec(:,1+3)/max(dr_dec(:,1+3)),'r-o')
hold off
xlabel('IPTG (\mu M)')
ylabel('Relative level of LacI')
set(gca,'XScale','log','YScale','lin')


data = xlsread('Data_LacIvsIPTG_Xu2009.xlsx','B2:C46');
figure(1)
hold on
plot(data(:,1),data(:,2),'x')
hold off


%

ss0 = [];
ic0 = [0,0,1,0,0,1]; % [Ix,I,L,C,T,Eg]

% Define inducer conc:
ic0(1) = 0;

% Integrate till steady state (within tol):
iter = 1;
while iter < 100
    % integrate over fixed time span:
    [t_t,x_t] = ode15s(@(t,x)Model3_TS_Dyn(t,x,params),[0,100],ic0);

    % check max derivative value:
    dx = Model3_TS_Dyn(t_t(end),x_t(end,:),params);
    if max(abs(dx)) > 1e-4
        ic0 = x_t(end,:);
        iter = iter + 1;
    else
        break
    end
end
if iter == 100
    disp('Did not reach steady state!')
end

% Save steady state levels:
ss0 = x_t(end,:);


%

% Time course dynamics for different inducer levels of IPTG:
rngInd = logspace(-2,2,5);

tspan = [0,0.005];
ic = ss0;
results = [];
figure(2); clf
for i = 1:length(rngInd)
    
    ic0 = ic;
    ic0(1) = rngInd(i);
    
    [t_t,x_t] = ode15s(@(t,x)Model3_TS_Dyn(t,x,params),tspan,ic0);
    results{i} = [t_t,x_t];
    
    hold on
    plot(t_t,x_t(:,3),'-','LineWidth',2)
    hold off
end
xlabel('Time (h)')
ylabel('LacI level')
set(gca,'XScale','lin')
