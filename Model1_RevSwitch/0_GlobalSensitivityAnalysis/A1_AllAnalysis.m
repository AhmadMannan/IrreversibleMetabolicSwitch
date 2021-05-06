% ------------------------------------------------------------------------
% PERFORMING eFAST ON ...
% ... BOTH CIRCUITS, ONE WITH NAR AND OTHER WITH PAR
% ... FOR EACH OF THE THREE DOSE-RESPONSE CHARACTERISTICS: BISTABLE RANGE
%     (BR), INDUCTION THRESHOLD (IT) AND REVERSION THRESHOLD (RT).
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% ------------------------------------------------------------------------

function A1_AllAnalysis


%%
% Plotting bar chart of sensitivity indices (based on eFAST) for ...
% ... both circuit topologies:
arch = {'NAR','PAR'};
% ... and for each of the three dose-response characteristics: bistable
% region (BR), induction threshold (IT), and reversion threshold (RT):
metric = {'IT','RT','BR'};

% Calculating sensitivities and plotting bar plot of first-order
% sensitivity indices:
count = 6;
s_KM = [];
for i = 1:length(arch)
    tic
    w = waitbar(0,['Looking at circuit with ',arch{i},'...']);
    for j = 1:length(metric)
        count = count + 1;
        s_KM{count} = A1_efast(arch{i},metric{j},count); %#ok<*AGROW>
        waitbar(j/length(metric))
    end
    close(w)
    toc
end
