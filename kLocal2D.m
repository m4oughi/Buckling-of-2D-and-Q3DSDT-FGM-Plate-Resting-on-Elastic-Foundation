function k_local = kLocal2D(p)

k_local = zeros(14, 14);

if p==0
    k_local = load('kLocal2DP0.mat');
elseif p==1
    k_local = load('kLocal2DP1.mat');
elseif p==2
    k_local = load('kLocal2DP2.mat');    
elseif p==3
    k_local = load('kLocal2DP3.mat');    
elseif p==4
    k_local = load('kLocal2DP4.mat');
elseif p==5
    k_local = load('kLocal2DP5.mat');
elseif p==10
    k_local = load('kLocal2DP10.mat');
else
end
end

