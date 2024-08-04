clear all
clc

%% A) Inputs
fileName = 'inputs.xlsx';
inputs = readtable(fileName);

if table2array(inputs(1, 2)) == "1" && table2array(inputs(1, 3)) == "0"
    plateTheory = "2D";
    disp("You are employing the 2D Plate Theory")
elseif table2array(inputs(1, 2)) == "0" && table2array(inputs(1, 3)) == "1"
    plateTheory = "Quasi-3D";
    disp("You are employing the Quasi-3D Plate Theory")
elseif table2array(inputs(1, 2)) == "0" && table2array(inputs(1, 3)) == "0"
    error("You didn't choose any Plate Theory. Please choose it from 'inputs.xlsx'.")
elseif table2array(inputs(1, 2)) == "1" && table2array(inputs(1, 3)) == "1"
    error("You cannot choose 2D and Quasi-3D as a Plate Theory in the same time! please modify the 'input.xlsx'")
else
    error("Please only choose 0 as negative and 1 as positive for Plate Theory in 'inputs.xlsx'!")
end


% --------  Modulus of Elasticity --------%
global p Ec Em noo G b total_a h m n ax ay axy kw_bar ks_bar;
 
p = table2array(inputs(4, 3));
Ec = table2array(inputs(5, 3));
Em = table2array(inputs(6, 3));

% --------     Poisson's Ratio    --------%
noo = table2array(inputs(7, 3));

% --------     Shear's Modulus    --------%
G = table2array(inputs(8, 3));

% --------        Dimensions      --------%
b       = table2array(inputs(9, 3));  % Length
total_a = table2array(inputs(10, 3));  % Width
h       = table2array(inputs(11, 3));  % height

% --------       Mode Number      --------%
m = 1; % constant

% --------    Number of Nodal lines    --------%
n = table2array(inputs(13, 3));

% --------       Geometric        --------%
ax  = table2array(inputs(14, 3));
ay  = table2array(inputs(15, 3));
axy = table2array(inputs(16, 3));

% -------- Elastic Foundation ---------%
kw_bar = table2array(inputs(17, 3));
ks_bar = table2array(inputs(18, 3));


% -------- Boundary Conditions ---------%
SSSS = table2array(inputs(1, 6));
SSSC = table2array(inputs(2, 6));
SSSF = table2array(inputs(3, 6));
SCSC = table2array(inputs(4, 6));
SCSF = table2array(inputs(5, 6));
SFSF = table2array(inputs(7, 6));

if SSSS==1 && SSSC==0 && SSSF==0 && SCSC==0 && SCSF==0 && SFSF==0
    boundaryCondition = "SSSS";
elseif SSSS==0 && SSSC==1 && SSSF==0 && SCSC==0 && SCSF==0 && SFSF==0
    boundaryCondition = "SSSC";
elseif SSSS==0 && SSSC==0 && SSSF==1 && SCSC==0 && SCSF==0 && SFSF==0
    boundaryCondition = "SSSF";
elseif SSSS==0 && SSSC==0 && SSSF==0 && SCSC==1 && SCSF==0 && SFSF==0
    boundaryCondition = "SCSC";
elseif SSSS==0 && SSSC==0 && SSSF==0 && SCSC==0 && SCSF==1 && SFSF==0
    boundaryCondition = "SCSF";
elseif SSSS==0 && SSSC==0 && SSSF==0 && SCSC==0 && SCSF==0 && SFSF==1
    boundaryCondition = "SFSF";
else
    error("Please choose only one boundary condition from 'inputs.xlsx' in each run!")
end


%%  D) Calculating Stiffness Matrix for Each Strip (14*14 for 2D && 18*18 for Quasi-3D)
if plateTheory == "2D"
    k_local = double(kLocal2D(p));
    disp("k_local")
    
    % Using Connectivity between strips to generate k global
    k__global = sym(zeros(8*n+6, 8*n+6));
    for j = 1:1:n
        k___global = sym(zeros(8*n+6, 8*n+6));
        k___global(8*j-7:8*j+6, 8*j-7:8*j+6) = k_local;
        k__global = k__global + k___global;
    end
    k_global = double(k__global);
    disp("k_global")
    
    
elseif plateTheory == "Quasi-3D"
    k_local = double(kLocal3D(p));
    disp("k_local")
    
    % Using Connectivity between strips to generate k global
    k__global = sym(zeros(10*n+8, 10*n+8));
    for j = 1:1:n
        k___global = sym(zeros(10*n+8, 10*n+8));
        k___global(10*j-9:10*j+8, 10*j-9:10*j+8) = k_local;
        k__global = k__global + k___global;
    end
    k_global = double(k__global);
    disp("k_global")
end




%% H) Calculating Geometric Stiffness Matrix (14*14 for 2D && 18*18 for Quasi-3D) and Global Geometric Stiffness Matrix (8*n+6)(8*n+6) for 2D and (10*n+8)*(10*n+8) for Quasi-3D

if plateTheory == "2D"
    kg_local = kgLocal2D(p);
    disp("kg_local")
    
    % Using Connectivity between strips to generate kg global
    kg__global = sym(zeros(8*n+6, 8*n+6));
    for j = 1:1:n
        kg___global = sym(zeros(8*n+6, 8*n+6));
        kg___global(8*j-7:8*j+6, 8*j-7:8*j+6) = kg_local;
        kg__global = kg__global + kg___global;
    end
    kg_global = double(kg__global);
    disp("kg_global")
    
    
elseif plateTheory == "Quasi-3D"
    kg_local = kgLocal3D(p);
    disp("kg_local")
    
    % Using Connectivity between strips to generate kg global
    kg__global = sym(zeros(10*n+8, 10*n+8));
    for j = 1:1:n
        kg___global = sym(zeros(10*n+8, 10*n+8));
        kg___global(10*j-9:10*j+8, 10*j-9:10*j+8) = kg_local;
        kg__global = kg__global + kg___global;
    end
    kg_global = double(kg__global);
    disp("kg_global")
end



%% Calculating Elastic Foundation Stiffness Matrices (Winkler and Pasternak)

if plateTheory == "2D"
    % Winkler
    kw_local = double(kwLocal2D(p));
    disp("kw_local")
    
    % Using Connectivity between strips to generate kw global
    kw__global = sym(zeros(8*n+6, 8*n+6));
    for j = 1:1:n
        kw___global = sym(zeros(8*n+6, 8*n+6));
        kw___global(8*j-7:8*j+6, 8*j-7:8*j+6) = kw_local;
        kw__global = kw__global + kw___global;
    end
    kw_global = double(kw__global);
    disp("kw_global")
    
    
    
    % Pasternak
    ks_local = double(ksLocal2D(p));
    disp("ks_local")
    
    % Using Connectivity between strips to generate ks global
    ks__global = sym(zeros(8*n+6, 8*n+6));
    for j = 1:1:n
        ks___global = sym(zeros(8*n+6, 8*n+6));
        ks___global(8*j-7:8*j+6, 8*j-7:8*j+6) = ks_local;
        ks__global = ks__global + ks___global;
    end
    ks_global = double(ks__global);
    disp("ks_global")
    
    
elseif plateTheory == "Quasi-3D"
    kw_local = double(kwLocal3D(p));
    disp("kw_local")
    
    % Using Connectivity between strips to generate kw global
    kw__global = sym(zeros(10*n+8, 10*n+8));
    for i = 1:1:n
        kw___global = sym(zeros(10*n+8, 10*n+8));
        kw___global(10*i-9:10*i+8, 10*i-9:10*i+8) = kw_local;
        kw__global = kw__global + kw___global;
    end
    kw_global = double(kw__global);
    disp("kw_global")
    
    % Pasternak
    ks_local = double(ksLocal3D(p));
    disp("ks_local")
    
    % Using Connectivity between strips to generate ks global
    ks__global = sym(zeros(10*n+8, 10*n+8));
    for i = 1:1:n
        ks___global = sym(zeros(10*n+8, 10*n+8));
        ks___global(10*i-9:10*i+8, 10*i-9:10*i+8) = ks_local;
        ks__global = ks__global + ks___global;
    end
    ks_global = double(ks__global);
    disp("ks_global")
end






%% K) Calculating Eigenvalue and Eigenvector
[Vec, Landa] = eig(k_global+kw_global+ks_global, kg_global);
disp("Eigenvalue")


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
%}