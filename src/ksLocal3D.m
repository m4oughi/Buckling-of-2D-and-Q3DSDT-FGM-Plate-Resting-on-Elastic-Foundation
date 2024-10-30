function ks_local = ksLocal3D(m, p, Eci, Emi, nooi, bi, total_ai, hi, ni, axi, ayi, axyi, kw_bari, ks_bari)

currentFolder = pwd;

if m==1
    if p==0
        filePath = fullfile(currentFolder, 'src', '3D', 'm1', 'p0', 'ksLocal.mat');
        
    elseif p==1
        filePath = fullfile(currentFolder, 'src', '3D', 'm1', 'p1', 'ksLocal.mat');
        
    elseif p==2
        filePath = fullfile(currentFolder, 'src', '3D', 'm1', 'p2', 'ksLocal.mat');
        
    elseif p==3
        filePath = fullfile(currentFolder, 'src', '3D', 'm1', 'p3', 'ksLocal.mat');
        
    elseif p==4
        filePath = fullfile(currentFolder, 'src', '3D', 'm1', 'p4', 'ksLocal.mat');
        
    elseif p==5
        filePath = fullfile(currentFolder, 'src', '3D', 'm1', 'p5', 'ksLocal.mat');
        
    elseif p==10
        filePath = fullfile(currentFolder, 'src', '3D', 'm1', 'p10', 'ksLocal.mat');
        
    else
        %pass
    end
    
    
elseif m==2
    if p==0
        filePath = fullfile(currentFolder, 'src', '3D', 'm2', 'p0', 'ksLocal.mat');
        
    elseif p==1
        filePath = fullfile(currentFolder, 'src', '3D', 'm2', 'p1', 'ksLocal.mat');
        
    elseif p==2
        filePath = fullfile(currentFolder, 'src', '3D', 'm2', 'p2', 'ksLocal.mat');
        
    elseif p==3
        filePath = fullfile(currentFolder, 'src', '3D', 'm2', 'p3', 'ksLocal.mat');
        
    elseif p==4
        filePath = fullfile(currentFolder, 'src', '3D', 'm2', 'p4', 'ksLocal.mat');
        
    elseif p==5
        filePath = fullfile(currentFolder, 'src', '3D', 'm2', 'p5', 'ksLocal.mat');
        
    elseif p==10
        filePath = fullfile(currentFolder, 'src', '3D', 'm2', 'p10', 'ksLocal.mat');
        
    else
        % Pass
    end
    
elseif m==3
    if p==0
        filePath = fullfile(currentFolder, 'src', '3D', 'm3', 'p0', 'ksLocal.mat');
        
    elseif p==1
        filePath = fullfile(currentFolder, 'src', '3D', 'm3', 'p1', 'ksLocal.mat');
        
    elseif p==2
        filePath = fullfile(currentFolder, 'src', '3D', 'm3', 'p2', 'ksLocal.mat');
        
    elseif p==3
        filePath = fullfile(currentFolder, 'src', '3D', 'm3', 'p3', 'ksLocal.mat');
        
    elseif p==4
        filePath = fullfile(currentFolder, 'src', '3D', 'm3', 'p4', 'ksLocal.mat');
        
    elseif p==5
        filePath = fullfile(currentFolder, 'src', '3D', 'm3', 'p5', 'ksLocal.mat');
        
    elseif p==10
        filePath = fullfile(currentFolder, 'src', '3D', 'm3', 'p10', 'ksLocal.mat');
    else
        % pass
    end
else
    % pass
end

if exist(filePath, 'file') == 2
    data = load(filePath);
else
    error('File does not exist: %s', filePath);
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

ks_local = data.ks_local;
values = {Eci, Emi, nooi, bi, total_ai, hi, ni, axi, ayi, axyi, kw_bari, ks_bari};

ks_local = subs(ks_local, {Ec, Em, noo, b, total_a, h, n, ax, ay, axy, kw_bar, ks_bar}, values);
ks_local = double(ks_local);
end