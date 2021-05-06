% ------------------------------------------------------------------------
% STUDY HOW TUNING CONTROL PARAMETERS AFFECTS THE OBJ1 VS OBJ2 PARETO
% FRONT, FOR CIRCUITS WITH NAR AND PAR
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% 12/11/2020
% ------------------------------------------------------------------------

function A2_RUN

%%
% Names and positions of all parameters:
% par_pos = [1,2,3,4,5,6,8,10,12,13];
% p_names = {'b_R','a_R','K_R','b_D','a_D','K_D','kcatD','Km_D','kcatB','Km_B','B','kf','kr','a_g','K_g'};

par_pos = [1,2,3,4,5,6];
p_names = {'b_R','a_R','K_R','b_D','a_D','K_D','kcatD','Km_D','kcatB','Km_B','B','kf','kr','a_g','K_g'};

tic
for i = 1:length(par_pos)
    disp('------------')
    disp(p_names(i))
    disp('------------')
    toc
    
    A2_ParamsEffectOnObj1VsObj2Curve(par_pos(i),i+4)
    set(gca,'XScale','log')
    axis([20,500,0,90])
end

%%
% for i = 1:length(par_pos)
%     figure(i+4)
%     subplot(1,2,1); ylim([0,125]); xlim([0,500])
%     subplot(1,2,2); ylim([0,125]); xlim([0,500])
% end
