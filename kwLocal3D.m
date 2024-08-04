function kw_local = kwLocal3D(p)

kw_local = zeros(18, 18);

if p==0
    kw_local = load('kwLocal3DP0.mat');
elseif p==1
    kw_local = load('kwLocal3DP1.mat');
elseif p==2
    kw_local = load('kwLocal3DP2.mat');    
elseif p==3
    kw_local = load('kwLocal3DP3.mat');    
elseif p==4
    kw_local = load('kwLocal3DP4.mat');
elseif p==5
    kw_local = load('kwLocal3DP5.mat');
elseif p==10
    kw_local = load('kwLocal3DP10.mat');
else
end
end