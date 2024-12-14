%% Constants and initializing
Tinf = 25; % Ambient temperature (°C)
Tb = 37; % Body temperature (°C)
h = 20; % Heat transfer coefficient (W/m^2·K)
R = 0.09; % Radius (cm)
theta = pi; % Maximum angle (rad)
nr = 101; % Number of nodes in the radial direction
Ntheta = 101; % Number of nodes in the angular direction
dr = R / (nr - 1);
dtheta = theta / (Ntheta - 1);
omega = 1.2;

% Blood perfusion and metabolic heat
wb = 0.00018; % Blood perfusion rate (1/s)
qm = 450; % Metabolic heat generation (W/m^3)
Rob = 1000; % Density of blood (kg/m^3)
cb = 4200; % Specific heat of blood (J/kg·K)

% Muscle properties
Rom = 1085; % Density of muscle (kg/m^3)
cm = 3800; % Specific heat of muscle (J/kg·K)
km = 0.42; % Thermal conductivity of muscle (W/m·K)

% tumour properties
Rob_T = 920;
Tb_T = 37.15;
qm_T = 0;
km_T = 0.56;
wb_T= 0;
cb_T = 3000;

% Initialize temperature array
T = Tb * ones(nr, Ntheta);
T_old = T;
error_in = 100;

% Set boundary conditions
T(1, :) = Tb; % Center (r = 0)
T(:, 1) = Tb; % theta = 0
T(:, end) = Tb; % theta = p
K = Rob*wb*cb;
K_T = Rob_T*wb_T*cb_T;

% Iterative solution parameters
tol = 1e-6;
max_iter = 10000;
iter = 0;

%% define the tumor geometry and location
% info required: rc - distance of tumor center from domain center
               % rt - radius of tumor
               % theta_c - angle of tumor center from horizontal axis
rc = 0.045;
rt = 0.03;
rc_node = round(rc/dr);
theta_c = pi/2;
theta_c_node = round(theta_c/dtheta);

% define searching limits for finding the nodes within the tumor

min_theta = round((theta_c - asin(rt/rc))/dtheta); % ceil rounds off to the greater integer 
max_theta = ceil((theta_c + asin(rt/rc))/dtheta);
min_r = round((rc-rt)/dr);
max_r = ceil((rc+rt)/dr);

%% Main time-stepping loop
while error_in > tol && iter < max_iter 
    for i = 2:nr-1
        for j = 2:Ntheta-1
            r = (i-1) * dr;
            theta = (j-1) * dtheta;

            coeff = K + 2*km/dr^2 + 2*km/(r^2*dtheta^2);
            coeff_T = K_T + 2*km_T/dr^2 + 2*km_T/(r^2*dtheta^2);

            dist = ((r*cos(theta) - rc*cos(theta_c))^2 + (r*sin(theta) - rc*sin(theta_c))^2)^(1/2);

            % Finite difference approximation of the bioheat equation
            if dist < rt
                T(i,j) = (1/coeff_T)*((km_T/dr^2 + km_T/(2*r*dr)) * T(i+1,j) + (km_T/dr^2 - km_T/(2*r*dr)) * T(i-1,j) + km_T/(r^2*dtheta^2) * (T(i,j+1) + T(i,j-1)) + K_T*Tb_T + qm_T);
            else
                T(i,j) = (1/coeff)*((km/dr^2 + km/(2*r*dr)) * T(i+1,j) + (km/dr^2 - km/(2*r*dr)) * T(i-1,j) + km/(r^2*dtheta^2) * (T(i,j+1) + T(i,j-1)) + K*Tb + qm);
            end

            % Apply SOR
            T(i,j) = (1-omega) * T_old(i,j) + omega * T(i,j);
        end
    end

    % boundary conditions
    T(end,:) = (1/(1/dr + h/km)) * (T(end-1,:)/dr + Tinf*h/km);

    error_in = max(max(abs(T_old - T)));
    iter = iter + 1;
    T_old = T;
end