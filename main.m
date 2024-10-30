clear all
clc

% Add the subfolder to the path
addpath(fullfile(pwd, 'src'));


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
global p m Ec Em Roc Rom noo b total_a h n ax ay axy kw_bar ks_bar;

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


%%  D) Calculating Stiffness Matrix for Each Strip
if plateTheory == "2D"
    k_local = kLocal2D(m, p, Ec, Em, noo, b, total_a, h, n, ax, ay, axy, kw_bar, ks_bar);
    disp("k_local")
    
elseif plateTheory == "Quasi-3D"
    k_local = kLocal3D(m, p, Ec, Em, noo, b, total_a, h, n, ax, ay, axy, kw_bar, ks_bar);
    disp("k_local")
end

%% H) Calculating Geometric Stiffness Matrix (14*14 for 2D && 18*18 for Quasi-3D) and Global Geometric Stiffness Matrix (8*n+6)(8*n+6) for 2D and (10*n+8)*(10*n+8) for Quasi-3D
if plateTheory == "2D"
    kg_local = kgLocal2D(m, p, Ec, Em, noo, b, total_a, h, n, ax, ay, axy, kw_bar, ks_bar);
    disp("kg_local")
    
elseif plateTheory == "Quasi-3D"
    kg_local = kgLocal3D(m, p, Ec, Em, noo, b, total_a, h, n, ax, ay, axy, kw_bar, ks_bar);
    disp("kg_local")
end

%% Calculating Elastic Foundation Stiffness Matrices (Winkler and Pasternak)
if plateTheory == "2D"
    % Winkler
    kw_local = kwLocal2D(m, p, Ec, Em, noo, b, total_a, h, n, ax, ay, axy, kw_bar, ks_bar);
    disp("kw_local")
    
    % Pasternak
    ks_local = ksLocal2D(m, p, Ec, Em, noo, b, total_a, h, n, ax, ay, axy, kw_bar, ks_bar);
    disp("ks_local")
    
elseif plateTheory == "Quasi-3D"
    % Winkler
    kw_local = kwLocal3D(m, p, Ec, Em, noo, b, total_a, h, n, ax, ay, axy, kw_bar, ks_bar);
    disp("kw_local")
    
    % Pasternak
    ks_local = ksLocal3D(m, p, Ec, Em, noo, b, total_a, h, n, ax, ay, axy, kw_bar, ks_bar);
    disp("ks_local")
end

%% Generate General stiffness matrix
%
kLocal = k_local + kw_local + ks_local;
disp("Sum of local, winkler, and pasternak stiffness matrices is calculated!")
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