function ks_local = ksLocal3D(p)

ks_local = zeros(18, 18);

if p==0
    ks_local = load('ksLocal3DP0.mat');
elseif p==1
    ks_local = load('ksLocal3DP1.mat');
elseif p==2
    ks_local = load('ksLocal3DP2.mat');    
elseif p==3
    ks_local = load('ksLocal3DP3.mat');    
elseif p==4
    ks_local = load('ksLocal3DP4.mat');
elseif p==5
    ks_local = load('ksLocal3DP5.mat');
elseif p==10
    ks_local = load('ksLocal3DP10.mat');
else
end
end