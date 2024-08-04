function kg_local = kgLocal2D(p)

kg_local = zeros(14, 14);

if p==0
    kg_local = load('kgLocal2DP0.mat');
elseif p==1
    kg_local = load('kgLocal2DP1.mat');
elseif p==2
    kg_local = load('kgLocal2DP2.mat');    
elseif p==3
    kg_local = load('kgLocal2DP3.mat');    
elseif p==4
    kg_local = load('kgLocal2DP4.mat');
elseif p==5
    kg_local = load('kgLocal2DP5.mat');
elseif p==10
    kg_local = load('kgLocal2DP10.mat');
else
end
end