function [muProp,Pprop] = gsfPropagateStep(muPrior,Pprior,Ad,Bd,u,q,Q)
muProp = cell(size(muPrior,1),length(q));
Pprop = cell(size(muPrior,1),length(q));
for ii = 1:size(muProp,1)
    
    mu = muPrior{ii};
    P = Pprior{ii};
    parfor jj = 1:size(muProp,2)
        
        % state propagation
        muProp{ii,jj} = Ad*mu + Bd*[0;u] + q(:,jj);
        
        % covariance propagation
        Pprop{ii,jj} = Ad*P*Ad' + Q;
        
    end
    
end