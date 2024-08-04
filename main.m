clear all
clc
disp("TWO")
%% A) Inputs
fileName = 'inputs.xlsx';
inputs = readtable(fileName);

% --------  Modulus of Elasticity --------%
global p Ec Em noo G b total_a h m n ax ay axy kw_bar ks_bar;

p = 0;
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


%%  D) Calculating Stiffness Matrix for Each Strip (18*18)
k_local = double(kLocal2D(p, Ec, b, total_a, h, n));
disp("k_local")


%% E) Calculating Global Stiffness Matrix (10*n+8)*(10*n+8)
k__global = sym(zeros(8*n+6, 8*n+6));
for j = 1:1:n
    k___global = sym(zeros(8*n+6, 8*n+6));
    k___global(8*j-7:8*j+6, 8*j-7:8*j+6) = k_local;
    k__global = k__global + k___global;
end
k_global = double(k__global);
disp("k_global")


%% F) Applying S-S Boundary Condition for Global Stiffness Matrix (10*n+4)*(10*n+4)
%
% S-S
%END
k_global(8*n+5, :) = [];
k_global(:, 8*n+5) = [];
k_global(8*n+3, :) = [];
k_global(:, 8*n+3) = [];
k_global(8*n+2, :) = [];
k_global(:, 8*n+2) = [];
k_global(8*n+1, :) = [];
k_global(:, 8*n+1) = [];

% FIRST
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
k_global(8*n+6, :) = [];
k_global(:, 8*n+6) = [];
k_global(8*n+5, :) = [];
k_global(:, 8*n+5) = [];
k_global(8*n+4, :) = [];
k_global(:, 8*n+4) = [];
k_global(8*n+3, :) = [];
k_global(:, 8*n+3) = [];
k_global(8*n+2, :) = [];
k_global(:, 8*n+2) = [];
k_global(8*n+1, :) = [];
k_global(:, 8*n+1) = [];

% FIRST
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
k_global(8*n+5, :) = [];
k_global(:, 8*n+5) = [];
k_global(8*n+3, :) = [];
k_global(:, 8*n+3) = [];
k_global(8*n+2, :) = [];
k_global(:, 8*n+2) = [];
k_global(8*n+1, :) = [];
k_global(:, 8*n+1) = [];

% FIRST
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
kw_local = double(kwLocal2D(p, Em, b, total_a, h, n, kw_bar));

kw__global = sym(zeros(8*n+6, 8*n+6));
for j = 1:1:n
    kw___global = sym(zeros(8*n+6, 8*n+6));
    kw___global(8*j-7:8*j+6, 8*j-7:8*j+6) = kw_local;
    kw__global = kw__global + kw___global;
end
kw_global = double(kw__global);
disp("kw_global")

%
% S-S
%END
kw_global(8*n+5, :) = [];
kw_global(:, 8*n+5) = [];
kw_global(8*n+3, :) = [];
kw_global(:, 8*n+3) = [];
kw_global(8*n+2, :) = [];
kw_global(:, 8*n+2) = [];
kw_global(8*n+1, :) = [];
kw_global(:, 8*n+1) = [];

% FIRST
kw_global(5, :) = [];
kw_global(:, 5) = [];
kw_global(3, :) = [];
kw_global(:, 3) = [];
kw_global(2, :) = [];
kw_global(:, 2) = [];
kw_global(1, :) = [];
kw_global(:, 1) = [];
disp("hazf")
%}
%}
% S-C
%



ks_local = ksLocal2D(p, Em, b, total_a, h, n, ks_bar);

ks__global = sym(zeros(8*n+6, 8*n+6));
for j = 1:1:n
    ks___global = sym(zeros(8*n+6, 8*n+6));
    ks___global(8*j-7:8*j+6, 8*j-7:8*j+6) = ks_local;
    ks__global = ks__global + ks___global;
end
ks_global = double(ks__global);
disp("ks_global")

%
% S-S
%END
ks_global(8*n+5, :) = [];
ks_global(:, 8*n+5) = [];
ks_global(8*n+3, :) = [];
ks_global(:, 8*n+3) = [];
ks_global(8*n+2, :) = [];
ks_global(:, 8*n+2) = [];
ks_global(8*n+1, :) = [];
ks_global(:, 8*n+1) = [];

% FIRST
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


%% H) Calculating Geometric Stiffness Matrix (18*18)
kg_local = kgLocal2D(p, b, total_a, h, n, ax, ay);


%% I) Calculating Global Stiffness Matrix (10*n+8)*(10*n+8)
kg__global = sym(zeros(8*n+6, 8*n+6));
for j = 1:1:n
    kg___global = sym(zeros(8*n+6, 8*n+6));
    kg___global(8*j-7:8*j+6, 8*j-7:8*j+6) = kg_local;
    kg__global = kg__global + kg___global;
end
kg_global = double(kg__global);
disp("kg_global")


%% J) S-S Boundary Condition for Global Geometric Stiffness Matrix (10*n+4)*(10*n+4)
%
% S-S
%END
kg_global(8*n+5, :) = [];
kg_global(:, 8*n+5) = [];
kg_global(8*n+3, :) = [];
kg_global(:, 8*n+3) = [];
kg_global(8*n+2, :) = [];
kg_global(:, 8*n+2) = [];
kg_global(8*n+1, :) = [];
kg_global(:, 8*n+1) = [];

% FIRST
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
kg_global(8*n+6, :) = [];
kg_global(:, 8*n+6) = [];
kg_global(8*n+5, :) = [];
kg_global(:, 8*n+5) = [];
kg_global(8*n+4, :) = [];
kg_global(:, 8*n+4) = [];
kg_global(8*n+3, :) = [];
kg_global(:, 8*n+3) = [];
kg_global(8*n+2, :) = [];
kg_global(:, 8*n+2) = [];
kg_global(8*n+1, :) = [];
kg_global(:, 8*n+1) = [];

% FIRST
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
kg_global(8*n+5, :) = [];
kg_global(:, 8*n+5) = [];
kg_global(8*n+3, :) = [];
kg_global(:, 8*n+3) = [];
kg_global(8*n+2, :) = [];
kg_global(:, 8*n+2) = [];
kg_global(8*n+1, :) = [];
kg_global(:, 8*n+1) = [];

% FIRST
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