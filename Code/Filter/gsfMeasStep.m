function [muUpd,Pupd] = updateGSF(muProp,Pprop,yCurr,C,r,R,gamma)
muUpd = cell(size(muProp,1),length(q));
Pupd = cell(size(muProp,1),length(q));

for ii = 1:size(muUpd,1)
    mu = muUpd{ii};
    P = Pupd{ii};
end