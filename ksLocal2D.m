function ks_local = ksLocal2D(p)

ks_local = zeros(14, 14);

if p==0
    ks_local = load('ksLocal2DP0.mat');
elseif p==1
    ks_local = load('ksLocal2DP1.mat');
elseif p==2
    ks_local = load('ksLocal2DP2.mat');    
elseif p==3
    ks_local = load('ksLocal2DP3.mat');    
elseif p==4
    ks_local = load('ksLocal2DP4.mat');
elseif p==5
    ks_local = load('ksLocal2DP5.mat');
elseif p==10
    ks_local = load('ksLocal2DP10.mat');
else
end
end