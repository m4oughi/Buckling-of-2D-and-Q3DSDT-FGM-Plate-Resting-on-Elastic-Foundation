clear all
clc

syms x y z
syms Ec Em b total_a h n ax ay axy kw_bar ks_bar

noo = 0.3;
m = 1;


%% B) Requirements
% --------   Width of Each Strip  --------%
a = total_a/n;

% --------   Flextural Rigidity  --------%
syms z x y


E = (Ec-Em)*(((z/h)+0.5)^10)+Em;
D = [
    [ (E*(noo - 1))/((2*noo - 1)*(noo + 1)),      -(E*noo)/((2*noo - 1)*(noo + 1)),      -(E*noo)/((2*noo - 1)*(noo + 1)),             0,             0,             0]
    [      -(E*noo)/((2*noo - 1)*(noo + 1)), (E*(noo - 1))/((2*noo - 1)*(noo + 1)),      -(E*noo)/((2*noo - 1)*(noo + 1)),             0,             0,             0]
    [      -(E*noo)/((2*noo - 1)*(noo + 1)),      -(E*noo)/((2*noo - 1)*(noo + 1)), (E*(noo - 1))/((2*noo - 1)*(noo + 1)),             0,             0,             0]
    [                                     0,                                     0,                                     0, E/(2*noo + 2),             0,             0]
    [                                     0,                                     0,                                     0,             0, E/(2*noo + 2),             0]
    [                                     0,                                     0,                                     0,             0,             0, E/(2*noo + 2)]
    ];


% ----  shape functions  ----%
syms z x y
N1 = (2*x-a)*(x-a)/(a^2);
N2 = (4*x*(a-x)/(a^2));
N3 = (x*(2*x-a))/(a^2);

F1 = 1-(3*(x/a)^2)+(2*(x/a)^3);
F3 = (3*(x/a)^2)-(2*(x/a)^3);
H1 = x*(1-2*(x/a)+(x/a)^2);
H3 = x*(-(x/a)+(x/a)^2);

s = sin(m*pi*y/b);
c = cos(m*pi*y/b);

% f = z*(1-((4/3)*((z^2)/(h^2))));
f = (((h/pi)*sinh(pi/h*z))-z)/(floor(cosh((pi/2)-1)));


%% C) Calculating Strain Matrix (6*18)
Bm = [
    diff(N1, x)*s   0   -z*diff(F1, x, 2)*s   -z*diff(H1, x, 2)*s   -f*diff(F1, x, 2)*s   -f*diff(H1, x, 2)*s   0   0            diff(N2, x)*s   0            diff(N3, x)*s   0   -z*diff(F3, x, 2)*s   -z*diff(H3, x, 2)*s   -f*diff(F3, x, 2)*s   -f*diff(H3, x, 2)*s   0   0;
    0   N1*diff(c, y)   -z*F1*diff(s, y, 2)   -z*H1*diff(s, y, 2)   -f*F1*diff(s, y, 2)   -f*H1*diff(s, y, 2)   0   0            0   N2*diff(c, y)            0   N3*diff(c, y)   -z*F3*diff(s, y, 2)   -z*H3*diff(s, y, 2)   -f*F3*diff(s, y, 2)   -f*H3*diff(s, y, 2)   0   0;
    0   0   0   0   0   0   -diff(f, z, 2)*F1*s   -diff(f, z, 2)*H1*s          0   0          0   0   0   0   0   0   -diff(f, z, 2)*F3*s   -diff(f, z, 2)*H3*s;
    N1*diff(s, y)   diff(N1, x)*c   -2*z*diff(F1, x)*diff(s, y)   -2*z*diff(H1, x)*diff(s, y)   -2*f*diff(F1, x)*diff(s, y)   -2*f*diff(H1, x)*diff(s, y)   0   0            N2*diff(s, y)   diff(N2, x)*c            N3*diff(s, y)   diff(N3, x)*c   -2*z*diff(F3, x)*diff(s, y)   -2*z*diff(H3, x)*diff(s, y)   -2*f*diff(F3, x)*diff(s, y)   -2*f*diff(H3, x)*diff(s, y)   0   0;
    0   0   0   0   (1-diff(f, z))*F1*diff(s, y)   (1-diff(f, z))*H1*diff(s, y)   (1-diff(f, z))*F1*diff(s, y)   (1-diff(f, z))*H1*diff(s, y)            0   0            0   0   0   0   (1-diff(f, z))*F3*diff(s, y)   (1-diff(f, z))*H3*diff(s, y)   (1-diff(f, z))*F3*diff(s, y)   (1-diff(f, z))*H3*diff(s, y);
    0   0   0   0   (1-diff(f, z))*diff(F1, x)*s   (1-diff(f, z))*diff(H1, x)*s   (1-diff(f, z))*diff(F1, x)*s   (1-diff(f, z))*diff(H1, x)*s            0   0            0   0   0   0   (1-diff(f, z))*diff(F3, x)*s   (1-diff(f, z))*diff(H3, x)*s   (1-diff(f, z))*diff(F3, x)*s   (1-diff(f, z))*diff(H3, x)*s;
    ];
disp("Bm")


%%  D) Calculating Stiffness Matrix for Each Strip (18*18)
k_local = int(int(int( Bm' * D * Bm, x, 0, a), y, 0, b), z, -0.5*h, 0.5*h);
disp("k_local")


%% Elastic Foundation
Dm = (Em*h^3)/(12*(1-noo^2));

Kw = (kw_bar*Dm)/(b^4);
Ks = (ks_bar*Dm)/(b^2);

N = [
    N1*s   0   -z*diff(F1, x)*s   -z*diff(H1, x)*s   -f*diff(F1, x)*s   -f*diff(H1, x)*s   0   0            N2*s   0          N3*s   0   -z*diff(F3, x)*s   -z*diff(H3, x)*s   -f*diff(F3, x)*s   -f*diff(H3, x)*s   0   0;
    0   N1*c   -z*F1*diff(s, y)   -z*H1*diff(s, y)   -f*F1*diff(s, y)   -f*H1*diff(s, y)   0   0            0   N2*c          0   N3*c   -z*F3*diff(s, y)   -z*H3*diff(s, y)   -f*F3*diff(s, y)   -f*H3*diff(s, y)   0   0;
    0   0   F1*s   H1*s   F1*s   H1*s   (1-diff(f, z))*F1*s   (1-diff(f, z))*H1*s                           0   0             0   0   F3*s   H3*s   F3*s   H3*s   (1-diff(f, z))*F3*s   (1-diff(f, z))*H3*s;
    ];

kw_local = int(int( subs(N'*Kw*N, z, 0), x, 0, a), y, 0, b);
disp("kw_global")


Nx = diff(N, x);
Ny = diff(N, y);
ksx_local = int(int( subs(Nx'*Ks*Nx, z, 0), x, 0, a), y, 0, b);
ksy_local = int(int( subs(Ny'*Ks*Ny, z, 0), x, 0, a), y, 0, b);
ks_local = ksx_local + ksy_local;

disp("ks_global")


%% G) the displacements field (18*1)
U = [N1*s;   0;   -z*diff(F1, x)*s;   -z*diff(H1, x)*s;   -f*diff(F1, x)*s;   -f*diff(H1, x)*s;   0;   0;            N2*s;   0;          N3*s;   0;   -z*diff(F3, x)*s;   -z*diff(H3, x)*s;   -f*diff(F3, x)*s;   -f*diff(H3, x)*s;   0;   0;];
V = [0;   N1*c;   -z*F1*diff(s, y);   -z*H1*diff(s, y);   -f*F1*diff(s, y);   -f*H1*diff(s, y);   0;   0;            0;   N2*c;          0;   N3*c;   -z*F3*diff(s, y);   -z*H3*diff(s, y);   -f*F3*diff(s, y);   -f*H3*diff(s, y);   0;   0;];
W = [0;   0;   F1*s;   H1*s;   F1*s;   H1*s;   (1-diff(f, z))*F1*s;   (1-diff(f, z))*H1*s;                           0;   0;             0;   0;   F3*s;   H3*s;   F3*s;   H3*s;   (1-diff(f, z))*F3*s;   (1-diff(f, z))*H3*s;];
disp("U V W")


%% H) Calculating Geometric Stiffness Matrix (18*18)
kg_local = int( int( int( ax*( diff(U, x).*diff(U, x)' + diff(V, x).*diff(V, x)' + diff(W, x).*diff(W, x)' ) + ay*( diff(U, y).* diff(U, y)' + diff(V, y).* diff(V, y)' + diff(W, y).* diff(W, y)' ) + axy*( diff(U, x).*diff(U, y)' + diff(U, x).*diff(U, y)' + diff(V, x).*diff(V, y)' + diff(V, x).*diff(V, y)' + diff(W, x).* diff(W, y)' + diff(W, x).*diff(W, y)' ), x, 0, a), y, 0, b), z, -0.5*h, 0.5*h);
disp("kg_local")

save('kLocal3DP10.mat', 'k_local')
save('kgLocal3DP10.mat', 'kg_local')
save('kwLocal3DP10.mat', 'kw_local')
save('ksLocal3DP10.mat', 'ks_local')

beep on
beep