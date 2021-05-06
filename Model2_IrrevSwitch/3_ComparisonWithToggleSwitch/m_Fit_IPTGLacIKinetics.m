function kOpts = m_Fit_IPTGLacIKinetics


options = optimoptions('fmincon','Display','iter');
lb = [1e2,1e2];
ub = [1e4,1e4];
[kOpts,fval,exitflag] = fmincon(@(p)m_Fit_objs(p),[800;1000],[],[],[],[],lb,ub,[],options);
disp(kOpts)
disp(fval)
disp(exitflag)

function obj = m_Fit_objs(p)
    
    % Define parameters:
%     p_nom = Model3_TS_Dyn_Params;
%     p_nom(17) = 1;
%     p_nom([6,7]) = p;
    
    % Simulations and data:
    [dataDR,dataK,DR,K] = m_TSmodelParameterisation(p);
    
    % Define objectives, i.e. least squares fitting:
    obj1 = sum((dataDR(:,2) - DR(:,2)).^2);
    obj2 = sum((dataK(:,2) - K(:,2)).^2);
    obj = obj1 + obj2;
    
end

end