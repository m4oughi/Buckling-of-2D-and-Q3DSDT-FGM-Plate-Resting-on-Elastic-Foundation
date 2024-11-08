%{
m=1 n=1 1  2  3  4  5  6   7  8   9  10 11 12 13 14

m=1 n=2 1  2  3  4  5  6   7  8   9  10 11 12 13 14
        9 10 11 12 13 14   15 16  17 18 19 20 21 22

m=1 n=3 1  2  3  4  5  6   7  8   9  10 11 12 13 14
        9  10 11 12 13 14  15 16  17 18 19 20 21 22
        17 18 19 20 21 22  23 24  25 26 27 28 29 30
_____________________________________________________________________________

m=2 n=1 1  2  3  4  5  6   7  8    9  10 11 12 13 14 | 15 16 17 18 19 20 21 22   23 24 25 26 27 28

m=2 n=2 1  2  3  4  5  6   7  8    9  10 11 12 13 14 | 23 24 25 26 27 28   29 30   31 32 33 34 35 36
        9 10 11 12 13 14   15 16   17 18 19 20 21 22 | 31 32 33 34 35 36   37 38   39 40 41 42 43 44

m=2 n=3 1  2  3  4  5  6    7  8    9  10 11 12 13 14 | 31 32 33 34 35 36  37 38   39 40 41 42 43 44
        9  10 11 12 13 14   15 16   17 18 19 20 21 22 | 39 40 41 42 43 44  45 46   47 48 49 50 51 52
        17 18 19 20 21 22   23 24   25 26 27 28 29 30 | 47 48 49 50 51 52  53 54   55 56 57 58 59 60
%}
clear all
clc


%% A) Inputs
fileName = 'inputs.xlsx';
inputs = readtable(fileName);

if table2array(inputs(1, 3)) == 1 && table2array(inputs(1, 4)) == 0
    plateTheory = "2D";
    disp("You are employing the 2D Plate Theory")
elseif table2array(inputs(1, 3)) == 0 && table2array(inputs(1, 4)) == 1
    plateTheory = "Quasi-3D";
    disp("You are employing the Quasi-3D Plate Theory")
elseif table2array(inputs(1, 3)) == 0 && table2array(inputs(1, 4)) == 0
    error("You didn't choose any Plate Theory. Please choose it from 'inputs.xlsx'.")
elseif table2array(inputs(1, 3)) == 1 && table2array(inputs(1, 4)) == 1
    error("You cannot choose 2D and Quasi-3D as a Plate Theory in the same time! please modify the 'input.xlsx'")
else
    error("Please only choose 0 as negative and 1 as positive for Plate Theory in 'inputs.xlsx'!")
end


% --------  Modulus of Elasticity --------%
p = table2array(inputs(4, 3));
Rom = table2array(inputs(5, 3));
Roc = table2array(inputs(6, 3));
Ec = table2array(inputs(7, 3));
Em = table2array(inputs(8, 3));

% --------     Poisson's Ratio    --------%
noo = table2array(inputs(9, 3));

% --------        Dimensions      --------%
b       = table2array(inputs(11, 3));  % Length
total_a = table2array(inputs(12, 3));  % Width
h       = table2array(inputs(13, 3));  % height

% --------       Mode Number      --------%
m = table2array(inputs(15, 3));

% --------    Number of Nodal lines    --------%
n = table2array(inputs(16, 3));

% --------       Geometric        --------%
ax  = table2array(inputs(17, 3));
ay  = table2array(inputs(18, 3));
axy = table2array(inputs(19, 3));

% -------- Elastic Foundation ---------%
kw_bar = table2array(inputs(20, 3));
ks_bar = table2array(inputs(21, 3));

% -------- Boundary Conditions ---------%
SSSS = table2array(inputs(1, 6));
SSSC = table2array(inputs(2, 6));
SSSF = table2array(inputs(3, 6));
SCSC = table2array(inputs(4, 6));
SCSF = table2array(inputs(5, 6));
SFSF = table2array(inputs(6, 6));

if SSSS==1 && SSSC==0 && SSSF==0 && SCSC==0 && SCSF==0 && SFSF==0
    boundaryCondition = "SSSS";
    disp("You pick the SSSS boundary condition")
elseif SSSS==0 && SSSC==1 && SSSF==0 && SCSC==0 && SCSF==0 && SFSF==0
    boundaryCondition = "SSSC";
    disp("You pick the SSSC boundary condition")
elseif SSSS==0 && SSSC==0 && SSSF==1 && SCSC==0 && SCSF==0 && SFSF==0
    boundaryCondition = "SSSF";
    disp("You pick the SSSF boundary condition")
elseif SSSS==0 && SSSC==0 && SSSF==0 && SCSC==1 && SCSF==0 && SFSF==0
    boundaryCondition = "SCSC";
    disp("You pick the SCSC boundary condition")
elseif SSSS==0 && SSSC==0 && SSSF==0 && SCSC==0 && SCSF==1 && SFSF==0
    boundaryCondition = "SCSF";
    disp("You pick the SCSF boundary condition")
elseif SSSS==0 && SSSC==0 && SSSF==0 && SCSC==0 && SCSF==0 && SFSF==1
    boundaryCondition = "SFSF";
    disp("You pick the SFSF boundary condition")
else
    error("Please choose only one boundary condition from 'inputs.xlsx' in each run!")
end


disp(n)
%% B) Requirements
% --------   Width of Each Strip  --------%
a = total_a/n;

% --------   Flextural Rigidity  --------%
syms x y z
E = (Ec-Em)*(((z/h)+0.5)^p)+Em;
%
if plateTheory == "2D"
    D = [
        (E/(1-noo^2))        noo*(E/(1-noo^2))      0               0                0;
        noo*(E/(1-noo^2))    (E/(1-noo^2))          0               0                0;
        0                     0                     E/(2*(1+noo))   0                0;
        0                     0                     0               E/(2*(1+noo))    0;
        0                     0                     0               0                E/(2*(1+noo));
        ];
    
elseif plateTheory == "Quasi-3D"
    D = [
        [ (E*(noo - 1))/((2*noo - 1)*(noo + 1)),      -(E*noo)/((2*noo - 1)*(noo + 1)),      -(E*noo)/((2*noo - 1)*(noo + 1)),             0,             0,             0]
        [      -(E*noo)/((2*noo - 1)*(noo + 1)), (E*(noo - 1))/((2*noo - 1)*(noo + 1)),      -(E*noo)/((2*noo - 1)*(noo + 1)),             0,             0,             0]
        [      -(E*noo)/((2*noo - 1)*(noo + 1)),      -(E*noo)/((2*noo - 1)*(noo + 1)), (E*(noo - 1))/((2*noo - 1)*(noo + 1)),             0,             0,             0]
        [                                     0,                                     0,                                     0, E/(2*noo + 2),             0,             0]
        [                                     0,                                     0,                                     0,             0, E/(2*noo + 2),             0]
        [                                     0,                                     0,                                     0,             0,             0, E/(2*noo + 2)]
        ];
    
end

Dm = (Em*h^3)/(12*(1-noo^2));

Kw = (kw_bar*Dm)/(b^4);
Ks = (ks_bar*Dm)/(b^2);
%}

% ----  shape functions  ----%
syms z x y
N1 = (2*x-a)*(x-a)/(a^2);
N2 = (4*x*(a-x)/(a^2));
N3 = (x*(2*x-a))/(a^2);

F1 = 1-(3*(x/a)^2)+(2*(x/a)^3);
F3 = (3*(x/a)^2)-(2*(x/a)^3);
H1 = x*(1-2*(x/a)+(x/a)^2);
H3 = x*(-(x/a)+(x/a)^2);

z0 = int(E*z, z, -0.5*h, 0.5*h)/int(E, z, -0.5*h, 0.5*h);
f = (((h/pi)*sinh(pi/h*z))-z)/(floor(cosh((pi/2)-1)));

z_ = z-z0;
f_ = subs(f, z, (z-z0));


%% C) Calculating Strain and displacement fields Matrices
%
B = []; N = []; U = []; V = []; W = [];
if plateTheory == "2D"
    for i=1:1:m
        s = sin(i.*pi.*y./b);
        c = cos(i.*pi.*y./b);
        
        Bm = [
            diff(N1, x)*s   0   -z_*diff(F1, x, 2)*s   -z_*diff(H1, x, 2)*s   -f_*diff(F1, x, 2)*s   -f_*diff(H1, x, 2)*s            diff(N2, x)*s   0            diff(N3, x)*s   0   -z_*diff(F3, x, 2)*s   -z_*diff(H3, x, 2)*s   -f_*diff(F3, x, 2)*s   -f_*diff(H3, x, 2)*s;
            0   N1*diff(c, y)   -z_*F1*diff(s, y, 2)   -z_*H1*diff(s, y, 2)   -f_*F1*diff(s, y, 2)   -f_*H1*diff(s, y, 2)            0   N2*diff(c, y)            0   N3*diff(c, y)   -z_*F3*diff(s, y, 2)   -z_*H3*diff(s, y, 2)   -f_*F3*diff(s, y, 2)   -f_*H3*diff(s, y, 2);
            N1*diff(s, y)   diff(N1, x)*c   -2*z_*diff(F1, x)*diff(s, y)   -2*z_*diff(H1, x)*diff(s, y)   -2*f_*diff(F1, x)*diff(s, y)   -2*f_*diff(H1, x)*diff(s, y)            N2*diff(s, y)   diff(N2, x)*c            N3*diff(s, y)   diff(N3, x)*c   -2*z_*diff(F3, x)*diff(s, y)   -2*z_*diff(H3, x)*diff(s, y)   -2*f_*diff(F3, x)*diff(s, y)   -2*f_*diff(H3, x)*diff(s, y);
            0   0   0   0   (1-diff(f_, z))*F1*diff(s, y)   (1-diff(f_, z))*H1*diff(s, y)            0   0            0   0   0   0   (1-diff(f_, z))*F3*diff(s, y)   (1-diff(f_, z))*H3*diff(s, y);
            0   0   0   0   (1-diff(f_, z))*diff(F1, x)*s   (1-diff(f_, z))*diff(H1, x)*s            0   0            0   0   0   0   (1-diff(f_, z))*diff(F3, x)*s   (1-diff(f_, z))*diff(H3, x)*s;
            ];
        
        Nm = [
            N1*s   0   -z_*diff(F1, x)*s   -z_*diff(H1, x)*s   -f_*diff(F1, x)*s   -f_*diff(H1, x)*s            N2*s   0          N3*s   0   -z_*diff(F3, x)*s   -z_*diff(H3, x)*s   -f_*diff(F3, x)*s   -f_*diff(H3, x)*s;
            0   N1*c   -z_*F1*diff(s, y)   -z_*H1*diff(s, y)   -f_*F1*diff(s, y)   -f_*H1*diff(s, y)            0   N2*c          0   N3*c   -z_*F3*diff(s, y)   -z_*H3*diff(s, y)   -f_*F3*diff(s, y)   -f*H3*diff(s, y);
            0   0   F1*s   H1*s   F1*s   H1*s          0   0                                                0   0             F3*s   H3*s   F3*s   H3*s;
            ];
        
        Um = [N1*s   0   -z_*diff(F1, x)*s   -z_*diff(H1, x)*s   -f_*diff(F1, x)*s   -f_*diff(H1, x)*s            N2*s   0          N3*s   0   -z_*diff(F3, x)*s   -z_*diff(H3, x)*s   -f_*diff(F3, x)*s   -f_*diff(H3, x)*s;];
        Vm = [0   N1*c   -z_*F1*diff(s, y)   -z_*H1*diff(s, y)   -f_*F1*diff(s, y)   -f_*H1*diff(s, y)            0      N2*c       0   N3*c   -z_*F3*diff(s, y)   -z_*H3*diff(s, y)   -f_*F3*diff(s, y)   -f_*H3*diff(s, y);];
        Wm = [0   0   F1*s   H1*s   F1*s   H1*s                                                               0   0             0   0   F3*s   H3*s   F3*s   H3*s;];
        
        B = [B, Bm];
        N = [N, Nm];
        U = [U, Um];
        V = [V, Vm];
        W = [W, Wm];
    end
    
elseif plateTheory == "Quasi-3D"
    for j=1:1:m
        s = sin(j.*pi.*y./b);
        c = cos(j.*pi.*y./b);
        
        Bm = [
            diff(N1, x)*s   0   -z_*diff(F1, x, 2)*s   -z_*diff(H1, x, 2)*s   -f_*diff(F1, x, 2)*s   -f_*diff(H1, x, 2)*s   0   0            diff(N2, x)*s   0            diff(N3, x)*s   0   -z_*diff(F3, x, 2)*s   -z_*diff(H3, x, 2)*s   -f_*diff(F3, x, 2)*s   -f_*diff(H3, x, 2)*s   0   0;
            0   N1*diff(c, y)   -z_*F1*diff(s, y, 2)   -z_*H1*diff(s, y, 2)   -f_*F1*diff(s, y, 2)   -f_*H1*diff(s, y, 2)   0   0            0   N2*diff(c, y)            0   N3*diff(c, y)   -z_*F3*diff(s, y, 2)   -z_*H3*diff(s, y, 2)   -f_*F3*diff(s, y, 2)   -f*H3*diff(s, y, 2)   0   0;
            0   0   0   0   0   0   -diff(f_, z, 2)*F1*s   -diff(f_, z, 2)*H1*s          0   0          0   0   0   0   0   0   -diff(f_, z, 2)*F3*s   -diff(f_, z, 2)*H3*s;
            N1*diff(s, y)   diff(N1, x)*c   -2*z_*diff(F1, x)*diff(s, y)   -2*z_*diff(H1, x)*diff(s, y)   -2*f_*diff(F1, x)*diff(s, y)   -2*f_*diff(H1, x)*diff(s, y)   0   0            N2*diff(s, y)   diff(N2, x)*c            N3*diff(s, y)   diff(N3, x)*c   -2*z_*diff(F3, x)*diff(s, y)   -2*z_*diff(H3, x)*diff(s, y)   -2*f_*diff(F3, x)*diff(s, y)   -2*f_*diff(H3, x)*diff(s, y)   0   0;
            0   0   0   0   (1-diff(f_, z))*F1*diff(s, y)   (1-diff(f_, z))*H1*diff(s, y)   (1-diff(f_, z))*F1*diff(s, y)   (1-diff(f_, z))*H1*diff(s, y)            0   0            0   0   0   0   (1-diff(f_, z))*F3*diff(s, y)   (1-diff(f_, z))*H3*diff(s, y)   (1-diff(f_, z))*F3*diff(s, y)   (1-diff(f_, z))*H3*diff(s, y);
            0   0   0   0   (1-diff(f_, z))*diff(F1, x)*s   (1-diff(f_, z))*diff(H1, x)*s   (1-diff(f_, z))*diff(F1, x)*s   (1-diff(f_, z))*diff(H1, x)*s            0   0            0   0   0   0   (1-diff(f_, z))*diff(F3, x)*s   (1-diff(f_, z))*diff(H3, x)*s   (1-diff(f_, z))*diff(F3, x)*s   (1-diff(f_, z))*diff(H3, x)*s;
            ];
        
        Nm = [
            N1*s   0   -z_*diff(F1, x)*s   -z_*diff(H1, x)*s   -f_*diff(F1, x)*s   -f_*diff(H1, x)*s   0   0            N2*s   0          N3*s   0   -z_*diff(F3, x)*s   -z_*diff(H3, x)*s   -f_*diff(F3, x)*s   -f_*diff(H3, x)*s   0   0;
            0   N1*c   -z_*F1*diff(s, y)   -z_*H1*diff(s, y)   -f_*F1*diff(s, y)   -f_*H1*diff(s, y)   0   0            0   N2*c          0   N3*c   -z_*F3*diff(s, y)   -z_*H3*diff(s, y)   -f_*F3*diff(s, y)   -f_*H3*diff(s, y)   0   0;
            0   0   F1*s   H1*s   F1*s   H1*s   (1-diff(f_, z))*F1*s   (1-diff(f_, z))*H1*s                           0   0             0   0   F3*s   H3*s   F3*s   H3*s   (1-diff(f_, z))*F3*s   (1-diff(f_, z))*H3*s;
            ];
        
        Um = [N1*s   0   -z_*diff(F1, x)*s   -z_*diff(H1, x)*s   -f_*diff(F1, x)*s   -f_*diff(H1, x)*s   0   0            N2*s   0          N3*s   0   -z_*diff(F3, x)*s   -z_*diff(H3, x)*s   -f_*diff(F3, x)*s   -f_*diff(H3, x)*s   0   0;];
        Vm = [0   N1*c   -z_*F1*diff(s, y)   -z_*H1*diff(s, y)   -f_*F1*diff(s, y)   -f_*H1*diff(s, y)   0   0            0   N2*c          0   N3*c   -z_*F3*diff(s, y)   -z_*H3*diff(s, y)   -f_*F3*diff(s, y)   -f_*H3*diff(s, y)   0   0;];
        Wm = [0   0   F1*s   H1*s   F1*s   H1*s   (1-diff(f_, z))*F1*s   (1-diff(f_, z))*H1*s                           0   0             0   0   F3*s   H3*s   F3*s   H3*s   (1-diff(f_, z))*F3*s   (1-diff(f_, z))*H3*s;];
        
        B = [B, Bm];
        N = [N, Nm];
        U = [U, Um];
        V = [V, Vm];
        W = [W, Wm];
    end
end

U = U.';
V = V.';
W = W.';
disp("B, N, U, V, W is generated!");
%}

%%  D) Calculating Stiffness Matrix for Each Strip (18*18)
%
k_local = double(int(int(int( B.'*D* B, x, 0, a), y, 0, b), z, -0.5*h, 0.5*h));
disp("k_local")
%}

%% Elastic Foundation
%
% Winkler
kw_local = int(int( subs(N.'*Kw*N, z, 0), x, 0, a), y, 0, b);
kw_local = double(kw_local);

% Pasternak
Nx = diff(N, x);
Ny = diff(N, y);
ksx_local = int(int( subs(Nx.'*Ks*Nx, z, 0), x, 0, a), y, 0, b);
ksy_local = int(int( subs(Ny.'*Ks*Ny, z, 0), x, 0, a), y, 0, b);
ks_local = double(ksx_local + ksy_local);
%}

%% Generate General stiffness matrix
%
kLocal = k_local + kw_local + ks_local;
disp("Sum of local, winkler, and pasternak stiffness matrices is calculated!")
%}

%% H) Calculating Geometric Stiffness Matrix (18*18)
%
kg_local = double(int( int( int( ax*( diff(U, x).*diff(U, x).' + diff(V, x).*diff(V, x).' + diff(W, x).*diff(W, x).' ) + ay*( diff(U, y).* diff(U, y).' + diff(V, y).* diff(V, y).' + diff(W, y).* diff(W, y).' ) + axy*( diff(U, x).*diff(U, y).' + diff(U, x).*diff(U, y).' + diff(V, x).*diff(V, y).' + diff(V, x).*diff(V, y).' + diff(W, x).* diff(W, y).' + diff(W, x).*diff(W, y).' ), x, 0, a), y, 0, b), z, -0.5*h, 0.5*h));
disp("kg_local")
%}

%% I) Calculating The Global Stiffness Matrix
%
if plateTheory == "2D"
    kltokg1 = zeros(14, 14, m, m);
    for i = 1:1:m
        for j = 1:1:m
            kltokg1(:, :, i, j) = kLocal(14*i-13:14*i, 14*j-13:14*j);
        end
    end
    
    kltokg3 = zeros(8*n+6, 8*n+6, m, m);
    for i = 1:1:m
        for j = 1:1:m
            for t = 1:1:n
                kltokg2 = zeros(8*n+6, 8*n+6, m, m);
                kltokg2(8*t-7:8*t+6, 8*t-7:8*t+6, i, j) = kltokg1(:, :, i, j);
                kltokg3(:, :, i, j) = kltokg3(:, :, i, j) + kltokg2(:, :, i, j);
            end
        end
    end
    
    k_global = zeros(m*(8*n+6), m*(8*n+6));
    for i = 1:1:m
        for j = 1:1:m
            k_global((8*n+6)*i-(8*n+5):(8*n+6)*i, (8*n+6)*j-(8*n+5):(8*n+6)*j) = kltokg3(:, :, i, j);
        end
    end
    
elseif plateTheory == "Quasi-3D"
    kltokg1 = zeros(18, 18, m, m);
    for i = 1:1:m
        for j = 1:1:m
            kltokg1(:, :, i, j) = kLocal(18*i-17:18*i, 18*j-17:18*j);
        end
    end
    
    kltokg2 = zeros(10*n+8, 10*n+8, m, m);
    for i = 1:1:m
        for j = 1:1:m
            for k = 1:1:n
                kltokg3 = zeros(10*n+8, 10*n+8, m, m);
                kltokg3(10*k-9:10*k+8, 10*k-9:10*k+8, i, j) = kltokg1(:, :, i, j);
                kltokg2(:, :, i, j) = kltokg2(:, :, i, j) + kltokg3(:, :, i, j);
            end
        end
    end
    
    k_global = zeros(m*(10*n+8), m*(10*n+8));
    for i = 1:1:m
        for j = 1:1:m
            k_global((10*n+8)*i-(10*n+7):(10*n+8)*i, (10*n+8)*j-(10*n+7):(10*n+8)*j) = kltokg2(:, :, i, j);
        end
    end
end


disp("Global stiffness matrix is calculated!")
%}


%% Calculating The Global Geometric Stiffness Matrix
%
if plateTheory == "2D"
    kgtokgg1 = zeros(14, 14, m, m);
    for i = 1:1:m
        for j = 1:1:m
            kgtokgg1(:, :, i, j) = kg_local(14*i-13:14*i, 14*j-13:14*j);
        end
    end
    
    kgtokgg3 = zeros(8*n+6, 8*n+6, m, m);
    for i = 1:1:m
        for j = 1:1:m
            for t = 1:1:n
                kgtokgg2 = zeros(8*n+6, 8*n+6, m, m);
                kgtokgg2(8*t-7:8*t+6, 8*t-7:8*t+6, i, j) = kgtokgg1(:, :, i, j);
                kgtokgg3(:, :, i, j) = kgtokgg3(:, :, i, j) + kgtokgg2(:, :, i, j);
            end
        end
    end
    
    kg_global = zeros(m*(8*n+6), m*(8*n+6));
    for i = 1:1:m
        for j = 1:1:m
            kg_global((8*n+6)*i-(8*n+5):(8*n+6)*i, (8*n+6)*j-(8*n+5):(8*n+6)*j) = kgtokgg3(:, :, i, j);
        end
    end
elseif plateTheory == "Quasi-3D"
    kgtokgg1 = zeros(18, 18, m, m);
    for i = 1:1:m
        for j = 1:1:m
            kgtokgg1(:, :, i, j) = kg_local(18*i-17:18*i, 18*j-17:18*j);
        end
    end
    
    kgtokgg2 = zeros(10*n+8, 10*n+8, m, m);
    for i = 1:1:m
        for j = 1:1:m
            for k = 1:1:n
                kgtokgg3 = zeros(10*n+8, 10*n+8, m, m);
                kgtokgg3(10*k-9:10*k+8, 10*k-9:10*k+8, i, j) = kgtokgg1(:, :, i, j);
                kgtokgg2(:, :, i, j) = kgtokgg2(:, :, i, j) + kgtokgg3(:, :, i, j);
            end
        end
    end
    
    kg_global = zeros(m*(10*n+8), m*(10*n+8));
    for i = 1:1:m
        for j = 1:1:m
            kg_global((10*n+8)*i-(10*n+7):(10*n+8)*i, (10*n+8)*j-(10*n+7):(10*n+8)*j) = kgtokgg2(:, :, i, j);
        end
    end
end

disp("Global geometry stiffness matrix is calculated!")
%}


%% F) Applying Boundary Condition for Global Stiffness Matrix
rm = [];
if plateTheory == "2D"
    if boundaryCondition == "SSSS"
        for i = m:-1:1
            rm = [rm, (8*n+6)*i-1]; %
            rm = [rm, (8*n+6)*i-3];
            rm = [rm, (8*n+6)*i-4];
            rm = [rm, (8*n+6)*i-5];
            
            rm = [rm, (8*n+6)*(i-1)+5];
            rm = [rm, (8*n+6)*(i-1)+3];
            rm = [rm, (8*n+6)*(i-1)+2];
            rm = [rm, (8*n+6)*(i-1)+1];
        end
    elseif boundaryCondition == "SSSC"
        for i = m:-1:1
            rm = [rm, (8*n+6)*i-1];
            rm = [rm, (8*n+6)*i-3];
            rm = [rm, (8*n+6)*i-4];
            rm = [rm, (8*n+6)*i-5];
            
            rm = [rm, (8*n+6)*(i-1)+6];
            rm = [rm, (8*n+6)*(i-1)+5];
            rm = [rm, (8*n+6)*(i-1)+4];
            rm = [rm, (8*n+6)*(i-1)+3];
            rm = [rm, (8*n+6)*(i-1)+2];
            rm = [rm, (8*n+6)*(i-1)+1];
        end
    elseif boundaryCondition == "SCSC"
        for i = m:-1:1
            rm = [rm, (8*n+6)*i];
            rm = [rm, (8*n+6)*i-1];
            rm = [rm, (8*n+6)*i-2];
            rm = [rm, (8*n+6)*i-3];
            rm = [rm, (8*n+6)*i-4];
            rm = [rm, (8*n+6)*i-5];
            
            rm = [rm, (8*n+6)*(i-1)+6];
            rm = [rm, (8*n+6)*(i-1)+5];
            rm = [rm, (8*n+6)*(i-1)+4];
            rm = [rm, (8*n+6)*(i-1)+3];
            rm = [rm, (8*n+6)*(i-1)+2];
            rm = [rm, (8*n+6)*(i-1)+1];
        end
    elseif boundaryCondition == "SSSF"
        for i = m:-1:1
            rm = [rm, (8*n+6)*i-1];
            rm = [rm, (8*n+6)*i-3];
            rm = [rm, (8*n+6)*i-4];
            rm = [rm, (8*n+6)*i-5];
        end
    elseif boundaryCondition == "SCSF"
        for i = m:-1:1
            rm = [rm, (8*n+6)*i];
            rm = [rm, (8*n+6)*i-1];
            rm = [rm, (8*n+6)*i-2];
            rm = [rm, (8*n+6)*i-3];
            rm = [rm, (8*n+6)*i-4];
            rm = [rm, (8*n+6)*i-5];
        end
    elseif boundaryCondition == "SFSF"
        % pass
    end
    
    
    
elseif plateTheory == "Quasi-3D"
    if boundaryCondition == "SSSS"
        for i = m:-1:1
            rm = [rm, (10*n+8)*i-1];
            rm = [rm, (10*n+8)*i-3];
            rm = [rm, (10*n+8)*i-5];
            rm = [rm, (10*n+8)*i-6];
            rm = [rm, (10*n+8)*i-7];
            
            rm = [rm, (10*n+8)*(i-1)+7];
            rm = [rm, (10*n+8)*(i-1)+5];
            rm = [rm, (10*n+8)*(i-1)+3];
            rm = [rm, (10*n+8)*(i-1)+2];
            rm = [rm, (10*n+8)*(i-1)+1];
        end
    elseif boundaryCondition == "SSSC"
        for i = m:-1:1
            rm = [rm, (10*n+8)*i-1];
            rm = [rm, (10*n+8)*i-3];
            rm = [rm, (10*n+8)*i-5];
            rm = [rm, (10*n+8)*i-6];
            rm = [rm, (10*n+8)*i-7];
            
            rm = [rm, (10*n+8)*(i-1)+8];
            rm = [rm, (10*n+8)*(i-1)+7];
            rm = [rm, (10*n+8)*(i-1)+6];
            rm = [rm, (10*n+8)*(i-1)+5];
            rm = [rm, (10*n+8)*(i-1)+4];
            rm = [rm, (10*n+8)*(i-1)+3];
            rm = [rm, (10*n+8)*(i-1)+2];
            rm = [rm, (10*n+8)*(i-1)+1];
        end
    elseif boundaryCondition == "SCSC"
        for i = m:-1:1
            rm = [rm, (10*n+8)*i];
            rm = [rm, (10*n+8)*i-1];
            rm = [rm, (10*n+8)*i-2];
            rm = [rm, (10*n+8)*i-3];
            rm = [rm, (10*n+8)*i-4];
            rm = [rm, (10*n+8)*i-5];
            rm = [rm, (10*n+8)*i-6];
            rm = [rm, (10*n+8)*i-7];
            
            rm = [rm, (10*n+8)*(i-1)+8];
            rm = [rm, (10*n+8)*(i-1)+7];
            rm = [rm, (10*n+8)*(i-1)+6];
            rm = [rm, (10*n+8)*(i-1)+5];
            rm = [rm, (10*n+8)*(i-1)+4];
            rm = [rm, (10*n+8)*(i-1)+3];
            rm = [rm, (10*n+8)*(i-1)+2];
            rm = [rm, (10*n+8)*(i-1)+1];
        end
    elseif boundaryCondition == "SSSF"
        for i = m:-1:1
            rm = [rm, (10*n+8)*i-1];
            rm = [rm, (10*n+8)*i-3];
            rm = [rm, (10*n+8)*i-5];
            rm = [rm, (10*n+8)*i-6];
            rm = [rm, (10*n+8)*i-7];
        end
        
    elseif boundaryCondition == "SCSF"
        for i = m:-1:1
            rm = [rm, (10*n+8)*i];
            rm = [rm, (10*n+8)*i-1];
            rm = [rm, (10*n+8)*i-2];
            rm = [rm, (10*n+8)*i-3];
            rm = [rm, (10*n+8)*i-4];
            rm = [rm, (10*n+8)*i-5];
            rm = [rm, (10*n+8)*i-6];
            rm = [rm, (10*n+8)*i-7];
        end
    elseif boundaryCondition == "SFSF"
        % pass
    end
end

disp(rm)
%
s_rm = size(rm);
for i = 1:1:s_rm(2)
    k_global(rm(i), :) = [];
    kg_global(rm(i), :) = [];
    
    k_global(:, rm(i)) = [];
    kg_global(:, rm(i)) = [];
end
%}

%
%% K) Calculating Eigenvalue and Eigenvector
[Vec, Landa] = eig(k_global, kg_global);
disp("Eigenvalue")


%% L) Critical Buckling Load Parameter(minimum)
Max_D = max(Landa);
dddd = min(Max_D);


%% M) Non-dimensionalization
bbb = dddd*((total_a)^2)*12*(1-noo^2) / ((pi^2)*(h^2)*Em);
fprintf('%.5f\n', bbb);
t1 = (total_a^2)/(h^2);
t2 = 12*(1-noo^2)/Em;
dddd*t1*t2
beep on
beep
%}