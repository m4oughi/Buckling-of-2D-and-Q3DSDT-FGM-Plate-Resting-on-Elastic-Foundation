function kg_local = kgLocal3D(p)

kg_local = zeros(18, 18);

if p==0
    kg_local = load('kgLocal3DP0.mat');
elseif p==1
    kg_local = load('kgLocal3DP1.mat');
elseif p==2
    kg_local = load('kgLocal3DP2.mat');
elseif p==3
    kg_local = load('kgLocal3DP3.mat');
elseif p==4
    kg_local = load('kgLocal3DP4.mat');
elseif p==5
    kg_local = load('kgLocal3DP5.mat');
elseif p==10
    kg_local = load('kgLocal3DP10.mat');
else
end
end