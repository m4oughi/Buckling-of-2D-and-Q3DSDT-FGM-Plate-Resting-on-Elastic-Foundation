function kw_local = kwLocal2D(p)

kw_local = zeros(14, 14);

if p==0
    kw_local = load('kwLocal2DP0.mat');
elseif p==1
    kw_local = load('kwLocal2DP1.mat');
elseif p==2
    kw_local = load('kwLocal2DP2.mat');    
elseif p==3
    kw_local = load('kwLocal2DP3.mat');    
elseif p==4
    kw_local = load('kwLocal2DP4.mat');
elseif p==5
    kw_local = load('kwLocal2DP5.mat');
elseif p==10
    kw_local = load('kwLocal2DP10.mat');
else
end
end