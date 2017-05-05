clear all
%% Extracellular Field

% muscle fibre properties
offset = 0;
fLen = offset+60; % total fibre length (mm)
zi = offset+35; % location of NMJ (mm)
A =96; % from Farina, D. (2001)
B =-90; % from Farina, D. (2001)
lambda=1; % to match implementation in Farina, D. (2001)
v = 30; % MFCV (mm/s)
L1 = 35;
L2 = 25;
a = 0.1;
sigma_i = 0.1;

% numerical settings
zinc = 0.1; % spatial increment

% Initialize arrays
z(int8(fLen/zinc)) = zeros;
Vm(int8(fLen/zinc)) = zeros;
dVm_dz(int8(fLen/zinc)) = zeros;
Im(int8(fLen/zinc)) = zeros;
k = 1;
z(1) = offset; 

% Compute spatial membrane potential and current (no propogation)
while z(k) < fLen
    % Voltage
    Vm(k) = A*z(k)^3*exp(-lambda*z(k))-B;
    dVm_dz(k) = A*z(k)^2*exp(-lambda*z(k))*(3-lambda*z(k)); 
    % Current
    Im(k) = lambda*z(k)*(6-6*lambda+lambda^2*z(k)^2)*exp(-lambda*z(k));
    
    z(k+1) = z(k) + zinc;
    k = k+1;
end

k = 1;
j = 1;
z_R(1) = zi ;

% Compute propogation of current
% spatial domain right of end plate (zi)
while z_R(k) <= zi + L2;
    % temporal domain 
    for t = 0:0.005:.25
        if z_R(k) >= -L1/2 || z_R(k) <= L2/2 
            phi_R(k,j) = A*(z_R(k)-zi-v*t)^2*exp(-lambda*(z_R(k)-zi-v*t))*(3-lambda*(z_R(k)-zi-v*t));
        else
            phi(k,j) = 0;
        end
        j = j+1;
    end     
    z_R(k+1) = z_R(k)+zinc;
    k = k+1; 
    j = 1;
end

% "hack" to avoid high values which seem to computed at any location
% "prior" to the start of the AP
for j = 1:size(phi_R,2)
    for k = 2:size(phi_R,1) 
        if phi_R(k,j) > phi_R(k-1,j)  
            phi_R(1:k-1,j) = 0;
            break
        end
    end
end

k = 1;
j = 1;
z_L(1) = 0;

% spatial domain left of the end plate (zi)
while z_L(k) <= zi;
    % temporal domain 
    for t = 0:0.01:.5
        phi_L(k,j) = -A*(-z_L(k)+zi-v*t)^2*exp(-lambda*(-z_L(k)+zi-v*t))*(3-lambda*(-z_L(k)+zi-v*t)); 
        j = j+1;
    end     
    z_L(k+1) = z_L(k)+zinc;
    k = k+1;
    j = 1;
end

% "hack" to avoid high values which seem to computed at any location
% "prior" to the start of the AP
for j = 1:size(phi_L,2)
    for k = size(phi_L,1)-1:-1:1 
    if phi_L(k-1,j) < phi_L(k,j)  
       phi_L(k:size(phi_L,1),j) = 0;
       break
    end
    end
end

%% Plotting

tinc = 30;

opt_grid = 'on';
opt_hold = 'on';
splotx = 2;
sploty = 2;

fnum = 1;
x = z_L(1:size(phi_L,1));
y = phi_L(:,tinc);
ftitle = 'Phi';
xtitle = 'Distance (mm)';
ytitle = 'Voltage (mV)';

plotxy(x,y, fnum, ftitle, xtitle, ytitle, opt_grid, opt_hold, ...
    splotx, sploty)

fnum = 1;
x = z_R(1:size(phi_R,1));
y = phi_R(:,tinc);
plotxy(x,y, fnum, ftitle, xtitle, ytitle, opt_grid, opt_hold, ...
    splotx, sploty)

fnum = 2;
zloc = 0.6*(fLen);
tol = 0.05;
% find corresponding index
for i = 1:length(z)
   if abs(z(i)-zloc) < tol  
        zloc_index = i;
        break
   end
end

if zloc < zi
    y = phi_L(zloc_index,:);
    
elseif zloc > zi
    y = phi_R(zloc_index-length(z_L),:);
end
x = 0:0.01:.5;
ftitle = 'Phi';
xtitle = 'Time (s)';
ytitle = 'Voltage (mV)';
plotxy(x,y, fnum, ftitle, xtitle, ytitle, opt_grid, opt_hold, ...
    splotx, sploty)

x = z(1:size(Vm,2));
y = Vm;
ftitle = 'Voltage';
xtitle = 'Distance (mm)';
ytitle = 'Voltage (mV)';
fnum = 3;
plotxy(x,y, fnum, ftitle, xtitle, ytitle, opt_grid, opt_hold, ...
    splotx, sploty)

y = Im;
ftitle = 'Current';
xtitle = 'Distance (mm)';
ytitle = 'Current (mV)';
fnum = 4;
plotxy(x,y, fnum, ftitle, xtitle, ytitle, opt_grid, opt_hold, ...
    splotx, sploty)

