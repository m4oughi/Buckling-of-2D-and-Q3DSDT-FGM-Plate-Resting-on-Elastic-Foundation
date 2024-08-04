function kg_local = kgLocal3D(p, Eci, Emi, nooi, bi, total_ai, hi, ni, axi, ayi, axyi, kw_bari, ks_bari)

if p==0
    data = load('kgLocal3DP0.mat');
elseif p==1
    data = load('kgLocal3DP1.mat');
elseif p==2
    data = load('kgLocal3DP2.mat');
elseif p==3
    data = load('kgLocal3DP3.mat');
elseif p==4
    data = load('kgLocal3DP4.mat');
elseif p==5
    data = load('kgLocal3DP5.mat');
elseif p==10
    data = load('kgLocal3DP10.mat');
else
end

Ec = sym('Ec');
Em = sym('Em');
noo = sym('noo');
b = sym('b');
total_a = sym('total_a');
h = sym('h');
n = sym('n');
ax = sym('ax');
ay = sym('ay');
axy = sym('axy');
kw_bar = sym('kw_bar');
ks_bar = sym('ks_bar');

kg_local = data.kg_local;
values = {Eci, Emi, nooi, bi, total_ai, hi, ni, axi, ayi, axyi, kw_bari, ks_bari};

kg_local = subs(kg_local, {Ec, Em, noo, b, total_a, h, n, ax, ay, axy, kw_bar, ks_bar}, values);
kg_local = double(kg_local);
end