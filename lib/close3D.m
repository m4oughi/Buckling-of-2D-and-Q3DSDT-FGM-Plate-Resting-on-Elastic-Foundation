clear all
clc

syms x y z
syms Ec Em Roc Rom b total_a h n ax ay axy kw_bar ks_bar

noo = 0.3;
m = 2;
p = 2;

%% B) Requirements
% --------   Width of Each Strip  --------%
a = total_a/n;

disp("The p in this loop is equal to " + num2str(p))

E = (Ec-Em)*(((z/h)+0.5)^p)+Em;
Ro = (Roc-Rom)*(((z/h)+0.5)^p)+Rom;

D = [
    [ (E*(noo - 1))/((2*noo - 1)*(noo + 1)),      -(E*noo)/((2*noo - 1)*(noo + 1)),      -(E*noo)/((2*noo - 1)*(noo + 1)),             0,             0,             0]
    [      -(E*noo)/((2*noo - 1)*(noo + 1)), (E*(noo - 1))/((2*noo - 1)*(noo + 1)),      -(E*noo)/((2*noo - 1)*(noo + 1)),             0,             0,             0]
    [      -(E*noo)/((2*noo - 1)*(noo + 1)),      -(E*noo)/((2*noo - 1)*(noo + 1)), (E*(noo - 1))/((2*noo - 1)*(noo + 1)),             0,             0,             0]
    [                                     0,                                     0,                                     0, E/(2*noo + 2),             0,             0]
    [                                     0,                                     0,                                     0,             0, E/(2*noo + 2),             0]
    [                                     0,                                     0,                                     0,             0,             0, E/(2*noo + 2)]
    ];

Dm = (Em*h^3)/(12*(1-noo^2));

Kw = (kw_bar*Dm)/(b^4);
Ks = (ks_bar*Dm)/(b^2);

% ----  shape functions  ----%
N1 = (2*x-a)*(x-a)/(a^2);
N2 = (4*x*(a-x)/(a^2));
N3 = (x*(2*x-a))/(a^2);

F1 = 1-(3*(x/a)^2)+(2*(x/a)^3);
F3 = (3*(x/a)^2)-(2*(x/a)^3);
H1 = x*(1-2*(x/a)+(x/a)^2);
H3 = x*(-(x/a)+(x/a)^2);

f = (((h/pi)*sinh(pi/h*z))-z)/(floor(cosh((pi/2)-1)));

%% C) Calculating Strain Matrix
B = [];
N = [];
U = [];
V = [];
W = [];

for j=1:1:m
    s = sin(j.*pi.*y./b);
    c = cos(j.*pi.*y./b);
    
    Bm = [
        diff(N1, x)*s   0   -z*diff(F1, x, 2)*s   -z*diff(H1, x, 2)*s   -f*diff(F1, x, 2)*s   -f*diff(H1, x, 2)*s   0   0            diff(N2, x)*s   0            diff(N3, x)*s   0   -z*diff(F3, x, 2)*s   -z*diff(H3, x, 2)*s   -f*diff(F3, x, 2)*s   -f*diff(H3, x, 2)*s   0   0;
        0   N1*diff(c, y)   -z*F1*diff(s, y, 2)   -z*H1*diff(s, y, 2)   -f*F1*diff(s, y, 2)   -f*H1*diff(s, y, 2)   0   0            0   N2*diff(c, y)            0   N3*diff(c, y)   -z*F3*diff(s, y, 2)   -z*H3*diff(s, y, 2)   -f*F3*diff(s, y, 2)   -f*H3*diff(s, y, 2)   0   0;
        0   0   0   0   0   0   -diff(f, z, 2)*F1*s   -diff(f, z, 2)*H1*s          0   0          0   0   0   0   0   0   -diff(f, z, 2)*F3*s   -diff(f, z, 2)*H3*s;
        N1*diff(s, y)   diff(N1, x)*c   -2*z*diff(F1, x)*diff(s, y)   -2*z*diff(H1, x)*diff(s, y)   -2*f*diff(F1, x)*diff(s, y)   -2*f*diff(H1, x)*diff(s, y)   0   0            N2*diff(s, y)   diff(N2, x)*c            N3*diff(s, y)   diff(N3, x)*c   -2*z*diff(F3, x)*diff(s, y)   -2*z*diff(H3, x)*diff(s, y)   -2*f*diff(F3, x)*diff(s, y)   -2*f*diff(H3, x)*diff(s, y)   0   0;
        0   0   0   0   (1-diff(f, z))*F1*diff(s, y)   (1-diff(f, z))*H1*diff(s, y)   (1-diff(f, z))*F1*diff(s, y)   (1-diff(f, z))*H1*diff(s, y)            0   0            0   0   0   0   (1-diff(f, z))*F3*diff(s, y)   (1-diff(f, z))*H3*diff(s, y)   (1-diff(f, z))*F3*diff(s, y)   (1-diff(f, z))*H3*diff(s, y);
        0   0   0   0   (1-diff(f, z))*diff(F1, x)*s   (1-diff(f, z))*diff(H1, x)*s   (1-diff(f, z))*diff(F1, x)*s   (1-diff(f, z))*diff(H1, x)*s            0   0            0   0   0   0   (1-diff(f, z))*diff(F3, x)*s   (1-diff(f, z))*diff(H3, x)*s   (1-diff(f, z))*diff(F3, x)*s   (1-diff(f, z))*diff(H3, x)*s;
        ];
    
    Nm = [
        N1*s   0   -z*diff(F1, x)*s   -z*diff(H1, x)*s   -f*diff(F1, x)*s   -f*diff(H1, x)*s   0   0            N2*s   0          N3*s   0   -z*diff(F3, x)*s   -z*diff(H3, x)*s   -f*diff(F3, x)*s   -f*diff(H3, x)*s   0   0;
        0   N1*c   -z*F1*diff(s, y)   -z*H1*diff(s, y)   -f*F1*diff(s, y)   -f*H1*diff(s, y)   0   0            0   N2*c          0   N3*c   -z*F3*diff(s, y)   -z*H3*diff(s, y)   -f*F3*diff(s, y)   -f*H3*diff(s, y)   0   0;
        0   0   F1*s   H1*s   F1*s   H1*s   (1-diff(f, z))*F1*s   (1-diff(f, z))*H1*s                           0   0             0   0   F3*s   H3*s   F3*s   H3*s   (1-diff(f, z))*F3*s   (1-diff(f, z))*H3*s;
        ];
    
    Um = [N1*s   0   -z*diff(F1, x)*s   -z*diff(H1, x)*s   -f*diff(F1, x)*s   -f*diff(H1, x)*s   0   0            N2*s   0          N3*s   0   -z*diff(F3, x)*s   -z*diff(H3, x)*s   -f*diff(F3, x)*s   -f*diff(H3, x)*s   0   0;];
    Vm = [0   N1*c   -z*F1*diff(s, y)   -z*H1*diff(s, y)   -f*F1*diff(s, y)   -f*H1*diff(s, y)   0   0            0   N2*c          0   N3*c   -z*F3*diff(s, y)   -z*H3*diff(s, y)   -f*F3*diff(s, y)   -f*H3*diff(s, y)   0   0;];
    Wm = [0   0   F1*s   H1*s   F1*s   H1*s   (1-diff(f, z))*F1*s   (1-diff(f, z))*H1*s                           0   0             0   0   F3*s   H3*s   F3*s   H3*s   (1-diff(f, z))*F3*s   (1-diff(f, z))*H3*s;];
    
    B = [B, Bm];
    N = [N, Nm];
    U = [U, Um];
    V = [V, Vm];
    W = [W, Wm];
end
U = U.';
V = V.';
W = W.';
disp("B, N, U, V, W is generated!");

%% H) Calculating Stiffness Matrix
%{
k_local0 = B' * D * B;
disp("B' * D * B is calculated")
k_local1 = int(k_local0, x, 0, a);
disp("int(k_local, x) is calculated")
clear k_local0
k_local2 = int(k_local1, y, 0, b);
disp("int(k_local, y) is calculated")
clear k_local1
k_local = int(k_local2, z, -0.5*h, 0.5*h);
disp("int(k_local, z) is calculated")
clear k_local2

txt1 = "kLocal3DP" + num2str(p)+ "m" + num2str(m) + ".mat";
save(txt1, 'k_local')

clear k_local
%}

%% H) Calculating Geometric Stiffness Matrix (18*18)
%{
kg_local_ax  = ax* ( diff(U, x).*diff(U, x)' + diff(V, x).*diff(V, x)' + diff(W, x).*diff(W, x)' );
disp("kg_local_ax")
kg_local_ay  = ay* ( diff(U, y).* diff(U, y)' + diff(V, y).* diff(V, y)' + diff(W, y).* diff(W, y)' );
disp("kg_local_ay")
kg_local_axy = axy*( diff(U, x).*diff(U, y)' + diff(U, x).*diff(U, y)' + diff(V, x).*diff(V, y)' + diff(V, x).*diff(V, y)' + diff(W, x).* diff(W, y)' + diff(W, x).*diff(W, y)');
disp("kg_local_axy")
kg_local00 = kg_local_ax + kg_local_ay + kg_local_axy;
disp("kg_local00")

kg_local0 = int(kg_local00, x, 0, a);
disp("kg_local0")
clear kg_local00

kg_local1 = int(kg_local0, y, 0, b);
disp("kg_local1")
clear kg_local0

kg_local = int(kg_local1, z, -0.5*h, 0.5*h);
disp("kg_local")
clear kg_local1


txt2 = "kgLocal3DP" + num2str(p) + "m" + num2str(m) + ".mat";
save(txt2, 'kg_local')
disp("kg_local is saved!")
clear kg_local
%}

%% Elastic Foundation
%{
%Winkler
kw_local0 = subs(N'*Kw*N, z, 0);
kw_local1 = int(kw_local0, x, 0, a);
clear kw_local0
kw_local = int(kw_local1, y, 0, b);
clear kw_local1
disp("Winkler stiffness matrix is calculated!")

txt3 = "kwLocal3DP" + num2str(p) + "m" + num2str(m) + ".mat";
save(txt3, 'kw_local')
disp("kw_local is saved!")
clear kw_local
%}
%
% Pasternak
Nx = diff(N, x);
disp("Nx is calculated");

ksx_local0 = subs(Nx'*Ks*Nx, z, 0);
ksx_local1 = int(ksx_local0, x, 0, a);
clear ksx_local0
ksx_local = int(ksx_local1, y, 0, b);
clear ksx_local1
disp("ksx_local");

Ny = diff(N, y);
disp("Ny is calculated");

ksy_local0 = subs(Ny'*Ks*Ny, z, 0);
ksy_local1 = int(ksy_local0, x, 0, a);
clear ksy_local0
ksy_local = int(ksy_local1, y, 0, b);
clear ksy_local1
disp("ksy_local");

ks_local = ksx_local+ksy_local;
disp("Sum of ksx_local and ksy_local is calculated!")

txt4 = "ksLocal3DP"+num2str(p)+ "m" + num2str(m) + ".mat";
save(txt4, 'ks_local')
disp("ks_local is saved!")
clear ksx_local
clear ksy_local
clear ks_local
%}
%{
m_local0 = int( N.'*Ro*N , x, 0, a);
m_local1 = int(m_local0, y, 0, b);
clear m_local0
m_local = int(m_local1, z, -0.5*h, 0.5*h);
clear m_local1
disp("m_local")

txt5 = "mLocal3DP"+num2str(p)+ "m" + num2str(m) + ".mat";
save(txt5, 'm_local')
disp("m_local is saved")
%}
beep on
beep
