% 0.10 --> 11.0779 11.1796 12.0950
% 0.15 --> 10.6069 10.7076 11.6141
% 0.20 --> 10.0137 10.1133 11.0095
clc
clear
disp("THREE")
%% A) Inputs
filename = 'inputs.xlsx';
inputs   = readtable(filename);


% --------  Modulus of Elasticity --------%
Ec = table2array(inputs(1, 3));
Em = table2array(inputs(2, 3));

% --------     Poisson's Ratio    --------%
noo = table2array(inputs(3, 3));

% --------     Shear's Modulus    --------%
G = table2array(inputs(4, 3));

% --------        Dimensions      --------%
b       = table2array(inputs(5, 3));  % Length
total_a = table2array(inputs(6, 3));  % Width
h       = table2array(inputs(7, 3));  % height

% --------       Mode Number      --------%
m = table2array(inputs(8, 3));

% --------    Number of Nodal lines    --------%
n = table2array(inputs(9, 3));
 

% --------       Geometric        --------%
ax  = table2array(inputs(10, 3));
ay  = table2array(inputs(11, 3));
axy = table2array(inputs(12, 3));

% -------- Elastic Foundation ---------%
kw_bar = table2array(inputs(13, 3));
ks_bar = table2array(inputs(14, 3));


disp(n)
%% B) Requirements
% --------   Width of Each Strip  --------%
a = total_a/n;

% --------   Flextural Rigidity  --------%
syms z x y
E = (Ec-Em)*(((z/h)+0.5)^0)+Em;
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
k_local = double(int(int(int( Bm' * D * Bm, x, 0, a), y, 0, b), z, -0.5*h, 0.5*h));
disp("k_local")


%% E) Calculating Global Stiffness Matrix (10*n+8)*(10*n+8)
k__global = sym(zeros(10*n+8, 10*n+8));
for j = 1:1:n
    k___global = sym(zeros(10*n+8, 10*n+8));
    k___global(10*j-9:10*j+8, 10*j-9:10*j+8) = k_local;
    k__global = k__global + k___global;
end
k_global = double(k__global);
disp("k_global")


%% F) Applying S-S Boundary Condition for Global Stiffness Matrix (10*n+4)*(10*n+4)
%
% S-S
%END
k_global(10*n+7, :) = [];
k_global(:, 10*n+7) = [];
k_global(10*n+5, :) = [];
k_global(:, 10*n+5) = [];
k_global(10*n+3, :) = [];
k_global(:, 10*n+3) = [];
k_global(10*n+2, :) = [];
k_global(:, 10*n+2) = [];
k_global(10*n+1, :) = [];
k_global(:, 10*n+1) = [];

% FIRST
k_global(7, :) = [];
k_global(:, 7) = [];
k_global(5, :) = [];
k_global(:, 5) = [];
k_global(3, :) = [];
k_global(:, 3) = [];
k_global(2, :) = [];
k_global(:, 2) = [];
k_global(1, :) = [];
k_global(:, 1) = [];
disp("hazf")
%}
%{
% C-C
%END
k_global(10*n+8, :) = [];
k_global(:, 10*n+8) = [];
k_global(10*n+7, :) = [];
k_global(:, 10*n+7) = [];
k_global(10*n+6, :) = [];
k_global(:, 10*n+6) = [];
k_global(10*n+5, :) = [];
k_global(:, 10*n+5) = [];
k_global(10*n+4, :) = [];
k_global(:, 10*n+4) = [];
k_global(10*n+3, :) = [];
k_global(:, 10*n+3) = [];
k_global(10*n+2, :) = [];
k_global(:, 10*n+2) = [];
k_global(10*n+1, :) = [];
k_global(:, 10*n+1) = [];

% FIRST
k_global(8, :) = [];
k_global(:, 8) = [];
k_global(7, :) = [];
k_global(:, 7) = [];
k_global(6, :) = [];
k_global(:, 6) = [];
k_global(5, :) = [];
k_global(:, 5) = [];
k_global(4, :) = [];
k_global(:, 4) = [];
k_global(3, :) = [];
k_global(:, 3) = [];
k_global(2, :) = [];
k_global(:, 2) = [];
k_global(1, :) = [];
k_global(:, 1) = [];
disp("hazf")
%}
%{
%S-C
%END
k_global(10*n+7, :) = [];
k_global(:, 10*n+7) = [];
k_global(10*n+5, :) = [];
k_global(:, 10*n+5) = [];
k_global(10*n+3, :) = [];
k_global(:, 10*n+3) = [];
k_global(10*n+2, :) = [];
k_global(:, 10*n+2) = [];
k_global(10*n+1, :) = [];
k_global(:, 10*n+1) = [];

% FIRST
k_global(8, :) = [];
k_global(:, 8) = [];
k_global(7, :) = [];
k_global(:, 7) = [];
k_global(6, :) = [];
k_global(:, 6) = [];
k_global(5, :) = [];
k_global(:, 5) = [];
k_global(4, :) = [];
k_global(:, 4) = [];
k_global(3, :) = [];
k_global(:, 3) = [];
k_global(2, :) = [];
k_global(:, 2) = [];
k_global(1, :) = [];
k_global(:, 1) = [];
disp("hazf")
%}
%{
% S-F
% FIRST
k_global(7, :) = [];
k_global(:, 7) = [];
k_global(5, :) = [];
k_global(:, 5) = [];
k_global(3, :) = [];
k_global(:, 3) = [];
k_global(2, :) = [];
k_global(:, 2) = [];
k_global(1, :) = [];
k_global(:, 1) = [];
disp("hazf")
%}
%{
%C-F
% FIRST
k_global(8, :) = [];
k_global(:, 8) = [];
k_global(7, :) = [];
k_global(:, 7) = [];
k_global(6, :) = [];
k_global(:, 6) = [];
k_global(5, :) = [];
k_global(:, 5) = [];
k_global(4, :) = [];
k_global(:, 4) = [];
k_global(3, :) = [];
k_global(:, 3) = [];
k_global(2, :) = [];
k_global(:, 2) = [];
k_global(1, :) = [];
k_global(:, 1) = [];
disp("hazf")
%}


%
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
kw_local = double(kw_local);
kw__global = sym(zeros(10*n+8, 10*n+8));
for i = 1:1:n
        kw___global = sym(zeros(10*n+8, 10*n+8));
        kw___global(10*i-9:10*i+8, 10*i-9:10*i+8) = kw_local;
        kw__global = kw__global + kw___global;
end
kw_global = double(kw__global);
disp("kw_global")

% S-S
%END

kw_global(10*n+7, :) = [];
kw_global(:, 10*n+7) = [];
kw_global(10*n+5, :) = [];
kw_global(:, 10*n+5) = [];
kw_global(10*n+3, :) = [];
kw_global(:, 10*n+3) = [];
kw_global(10*n+2, :) = [];
kw_global(:, 10*n+2) = [];
kw_global(10*n+1, :) = [];
kw_global(:, 10*n+1) = [];

% FIRST
kw_global(7, :) = [];
kw_global(:, 7) = [];
kw_global(5, :) = [];
kw_global(:, 5) = [];
kw_global(3, :) = [];
kw_global(:, 3) = [];
kw_global(2, :) = [];
kw_global(:, 2) = [];
kw_global(1, :) = [];
kw_global(:, 1) = [];
disp("hazf")

Nx = diff(N, x);
Ny = diff(N, y);
ksx_local = int(int( subs(Nx'*Ks*Nx, z, 0), x, 0, a), y, 0, b);
ksy_local = int(int( subs(Ny'*Ks*Ny, z, 0), x, 0, a), y, 0, b);
ks_local = double(ksx_local + ksy_local);
ks__global = sym(zeros(10*n+8, 10*n+8));
for i = 1:1:n
        ks___global = sym(zeros(10*n+8, 10*n+8));
        ks___global(10*i-9:10*i+8, 10*i-9:10*i+8) = ks_local;
        ks__global = ks__global + ks___global;
end
ks_global = double(ks__global);
disp("ks_global")

%END
ks_global(10*n+7, :) = [];
ks_global(:, 10*n+7) = [];
ks_global(10*n+5, :) = [];
ks_global(:, 10*n+5) = [];
ks_global(10*n+3, :) = [];
ks_global(:, 10*n+3) = [];
ks_global(10*n+2, :) = [];
ks_global(:, 10*n+2) = [];
ks_global(10*n+1, :) = [];
ks_global(:, 10*n+1) = [];

% FIRST
ks_global(7, :) = [];
ks_global(:, 7) = [];
ks_global(5, :) = [];
ks_global(:, 5) = [];
ks_global(3, :) = [];
ks_global(:, 3) = [];
ks_global(2, :) = [];
ks_global(:, 2) = [];
ks_global(1, :) = [];
ks_global(:, 1) = [];
disp("hazf")

%}

%% G) the displacements field (18*1)
U = [N1*s;   0;   -z*diff(F1, x)*s;   -z*diff(H1, x)*s;   -f*diff(F1, x)*s;   -f*diff(H1, x)*s;   0;   0;            N2*s;   0;          N3*s;   0;   -z*diff(F3, x)*s;   -z*diff(H3, x)*s;   -f*diff(F3, x)*s;   -f*diff(H3, x)*s;   0;   0;];
V = [0;   N1*c;   -z*F1*diff(s, y);   -z*H1*diff(s, y);   -f*F1*diff(s, y);   -f*H1*diff(s, y);   0;   0;            0;   N2*c;          0;   N3*c;   -z*F3*diff(s, y);   -z*H3*diff(s, y);   -f*F3*diff(s, y);   -f*H3*diff(s, y);   0;   0;];
W = [0;   0;   F1*s;   H1*s;   F1*s;   H1*s;   (1-diff(f, z))*F1*s;   (1-diff(f, z))*H1*s;                           0;   0;             0;   0;   F3*s;   H3*s;   F3*s;   H3*s;   (1-diff(f, z))*F3*s;   (1-diff(f, z))*H3*s;];
disp("U V W")


%% H) Calculating Geometric Stiffness Matrix (18*18)
kg_local = double(int( int( int( ax*( diff(U, x).*diff(U, x)' + diff(V, x).*diff(V, x)' + diff(W, x).*diff(W, x)' ) + ay*( diff(U, y).* diff(U, y)' + diff(V, y).* diff(V, y)' + diff(W, y).* diff(W, y)' ) + axy*( diff(U, x).*diff(U, y)' + diff(U, x).*diff(U, y)' + diff(V, x).*diff(V, y)' + diff(V, x).*diff(V, y)' + diff(W, x).* diff(W, y)' + diff(W, x).*diff(W, y)' ), x, 0, a), y, 0, b), z, -0.5*h, 0.5*h));
disp("kg_local")


%% I) Calculating Global Stiffness Matrix (10*n+8)*(10*n+8)
kg__global = sym(zeros(10*n+8, 10*n+8));
for j = 1:1:n
    kg___global = sym(zeros(10*n+8, 10*n+8));
    kg___global(10*j-9:10*j+8, 10*j-9:10*j+8) = kg_local;
    kg__global = kg__global + kg___global;
end
kg_global = double(kg__global);
disp("kg_global")


%% J) S-S Boundary Condition for Global Geometric Stiffness Matrix (10*n+4)*(10*n+4)
%
% S-S
%END
kg_global(10*n+7, :) = [];
kg_global(:, 10*n+7) = [];
kg_global(10*n+5, :) = [];
kg_global(:, 10*n+5) = [];
kg_global(10*n+3, :) = [];
kg_global(:, 10*n+3) = [];
kg_global(10*n+2, :) = [];
kg_global(:, 10*n+2) = [];
kg_global(10*n+1, :) = [];
kg_global(:, 10*n+1) = [];

% FIRST
kg_global(7, :) = [];
kg_global(:, 7) = [];
kg_global(5, :) = [];
kg_global(:, 5) = [];
kg_global(3, :) = [];
kg_global(:, 3) = [];
kg_global(2, :) = [];
kg_global(:, 2) = [];
kg_global(1, :) = [];
kg_global(:, 1) = [];
disp("hazf")
%}
%{
% C-C
%END
kg_global(10*n+8, :) = [];
kg_global(:, 10*n+8) = [];
kg_global(10*n+7, :) = [];
kg_global(:, 10*n+7) = [];
kg_global(10*n+6, :) = [];
kg_global(:, 10*n+6) = [];
kg_global(10*n+5, :) = [];
kg_global(:, 10*n+5) = [];
kg_global(10*n+4, :) = [];
kg_global(:, 10*n+4) = [];
kg_global(10*n+3, :) = [];
kg_global(:, 10*n+3) = [];
kg_global(10*n+2, :) = [];
kg_global(:, 10*n+2) = [];
kg_global(10*n+1, :) = [];
kg_global(:, 10*n+1) = [];

% FIRST
kg_global(8, :) = [];
kg_global(:, 8) = [];
kg_global(7, :) = [];
kg_global(:, 7) = [];
kg_global(6, :) = [];
kg_global(:, 6) = [];
kg_global(5, :) = [];
kg_global(:, 5) = [];
kg_global(4, :) = [];
kg_global(:, 4) = [];
kg_global(3, :) = [];
kg_global(:, 3) = [];
kg_global(2, :) = [];
kg_global(:, 2) = [];
kg_global(1, :) = [];
kg_global(:, 1) = [];
disp("hazf")
%}
%{
% C-S
%END
kg_global(10*n+7, :) = [];
kg_global(:, 10*n+7) = [];
kg_global(10*n+5, :) = [];
kg_global(:, 10*n+5) = [];
kg_global(10*n+3, :) = [];
kg_global(:, 10*n+3) = [];
kg_global(10*n+2, :) = [];
kg_global(:, 10*n+2) = [];
kg_global(10*n+1, :) = [];
kg_global(:, 10*n+1) = [];

% FIRST
kg_global(8, :) = [];
kg_global(:, 8) = [];
kg_global(7, :) = [];
kg_global(:, 7) = [];
kg_global(6, :) = [];
kg_global(:, 6) = [];
kg_global(5, :) = [];
kg_global(:, 5) = [];
kg_global(4, :) = [];
kg_global(:, 4) = [];
kg_global(3, :) = [];
kg_global(:, 3) = [];
kg_global(2, :) = [];
kg_global(:, 2) = [];
kg_global(1, :) = [];
kg_global(:, 1) = [];
disp("hazf")
%}
%{
% S-F
% FIRST
kg_global(7, :) = [];
kg_global(:, 7) = [];
kg_global(5, :) = [];
kg_global(:, 5) = [];
kg_global(3, :) = [];
kg_global(:, 3) = [];
kg_global(2, :) = [];
kg_global(:, 2) = [];
kg_global(1, :) = [];
kg_global(:, 1) = [];
disp("hazf")
%}
%{
%C-F
% FIRST
kg_global(8, :) = [];
kg_global(:, 8) = [];
kg_global(7, :) = [];
kg_global(:, 7) = [];
kg_global(6, :) = [];
kg_global(:, 6) = [];
kg_global(5, :) = [];
kg_global(:, 5) = [];
kg_global(4, :) = [];
kg_global(:, 4) = [];
kg_global(3, :) = [];
kg_global(:, 3) = [];
kg_global(2, :) = [];
kg_global(:, 2) = [];
kg_global(1, :) = [];
kg_global(:, 1) = [];
disp("hazf")
%}


%% K) Calculating Eigenvalue and Eigenvector
[Vec, Landa] = eig(k_global+kw_global+ks_global, kg_global);
% [Vec, Landa] = eig(k_global, kg_global);
disp("Landa") 


%% L) Critical Buckling Load Parameter(minimum)
Max_D = max(Landa);
dddd = min(Max_D);


%% M) Non-dimensionalization
% dddd*((total_a)^2)*12*(1-noo^2) / ((pi^2)*(h^2)*Em)
t1 = (total_a^2)/(h^2);
t2 = 12*(1-noo^2)/Em;
dddd*t1*t2
beep on
beep