function k_local = kLocal3D(p)

k_local = zeros(18, 18);

if p==0
    k_local = load('kLocal3DP0.mat');
elseif p==1
    k_local = load('kLocal3DP1.mat');
elseif p==2
    k_local = load('kLocal3DP2.mat');    
elseif p==3
    k_local = load('kLocal3DP3.mat');    
elseif p==4
    k_local = load('kLocal3DP4.mat');
elseif p==5
    k_local = load('kLocal3DP5.mat');
elseif p==10
    k_local = load('kLocal3DP10.mat');
else
end
end

