clear;
clc;
close;

% Physical Parameters
rho = 1;  % Density kg/m^3
cp = 10;  % Specific heat capacity J/(K*kg)
k = 0.12;  % Thermal conductivity W/(m*K)
nu = 1e-2;  % Kinematic viscosity m^2/s
Gamma = k/cp;  % Heat diffusion coefficient kg/(m*s)
L = 10;  % Channel length m
H = 1;  % Channel height m
Pe = 16.5;  % Peclet number
u_mean = Pe*Gamma/(2*H*rho);  % Average speed m/s
Re = 2*H*u_mean/nu;  % Reynolds number

% Mesh Parameters
nx = 50;  % Coarse mesh x refinement
ny = 5;  % Coarse mesh y refinement
dx = L/(nx-1);
dy = H/(ny-1);
N = nx*ny; % DOF number

% Boundary Temperature Parameters
T_in = 50;  % Inlet temperature °C
T_wall = 100;  % Wall temperature °C

% Velocity Field Function
u_x = @(y) 6*u_mean*(y/H)*(1-(y/H));

% Initialize coefficient matrix A and right-hand side vector b
A = zeros(N, N);
b = zeros(N, 1);

% Fill in the Equations for the internal CV
% Diffusion coefficient
D_e = Gamma*dy/dx;  % Eastern
D_w = Gamma*dy/dx;  % Western
D_n = Gamma*dx/dy;  % Northern
D_s = Gamma*dx/dy;  % Southern

for i = 2:(ny-1)
% Advection coefficient
y_e = H - (i-1)*dy;  % y_eastern
ux_e = u_x(y_e);  % Mean velocity
F_e = rho * ux_e * dy;  % Eastbound convection flux
F_w = F_e;  % Incompressible flow

% Coefficients of Algebraic Equations
a_E = D_e + max(-F_e,0);
a_W = D_w + max(F_w,0);
a_N = D_n;
a_S = D_s;
a_P = a_E + a_W + a_N + a_S;

    for j = 2:(nx-1)
        n = (i-1)*nx + j;

        % Fill matrix A (The internal CV has no source term, so b(n) is zero
        A(n,n) = a_P;
        A(n,n+1) = -a_E;
        A(n,n-1) = -a_W;
        A(n,n+nx) = -a_N;
        A(n,n-nx) = -a_S;
    end
end

% Fill boundary conditions at the inlet (Dirichlet T=T_in)
for i = 2:ny-1
    n = (i-1)*nx + 1;
    A(n,:) = 0;
    A(n,n) = 1;
    b(n) = T_in;
end

% Fill the boundary conditions at the outlet (Neumann)
for i = 2:ny-1
    n = (i-1)*nx + nx;

    % Advection coefficient
    y_e = H - (i-1)*dy;  % y_eastern
    ux_e = u_x(y_e);  % Mean velocity
    F_e = rho * ux_e * dy;  % Eastbound convection flux
    F_w = F_e;  % Incompressible flow
    
    % Coefficients of Algebraic Equations
    a_W = D_w + max(F_w,0);
    a_N = D_n;
    a_S = D_s;
    a_P = a_W + a_N + a_S;

    A(n,n) = a_P;
    A(n,n+1) = 0;
    A(n,n-1) = -a_W;
    A(n, n+nx) = -a_N;
    A(n, n-nx) = -a_S;
    b(n) = 0;
end

% Fill boundary conditions at the walls including 4 corner CVs (Dirichlet T=T_wall)
for j = 1:nx
    % Lower wall surface
    n = (ny-1)*nx + j;
    A(n,:) = 0;
    A(n,n) = 1;
    b(n) = T_wall;
    
    % Upper wall surface
    n = j;
    A(n,:) = 0;
    A(n,n) = 1;
    b(n) = T_wall;
end

% Solve the linear equation system
T = A\b;

T_2D = reshape(T, nx, ny)';
x_coords = linspace(0, L, nx);   
y_coords = linspace(0, H, ny);    
[X, Y] = meshgrid(x_coords, y_coords);

% Visualization and Validation of Temperature Fields
figure('Name', 'Coarse-mesh temperature field');
pcolor(X, Y, T_2D);   
shading interp;
colorbar;               
xlabel('x (m)');       
ylabel('y (m)');        
%title(['Channel temperature field（nx=' num2str(nx) ', ny=' num2str(ny) '）']);
axis equal tight;

%{
%%% 12.
% Verify inlet temperature (plotting the temperature curve at x=0)
figure('Name', 'Inlet Temperature Validation');
plot(y_coords, T_2D(:, 1), 'b-o'); 
ylim([T_in-1, T_in+1]);             
xlabel('y (m)');
ylabel('Temperature (°C)');
title('Temperature Distribution at the Inlet Boundary');
grid on;

% Verify wall surface temperature (plotting the temperature curve at y=0)
figure('Name', 'Lower Wall Temperature Verification');
plot(x_coords, T_2D(1, :), 'r-o'); 
ylim([T_wall-1, T_wall+1]);        
xlabel('x (m)');
ylabel('Temperature (°C)');
title('Temperature distribution on the lower wall surface');
grid on;


%%% 13. 
%visualizartion of the outlet temperature
figure('Name', 'Visualization of the outlet temperature');
plot(y_coords, T_2D(:, nx), 'r-o'); 
%ylim([T_wall-1, T_wall+1]);        
xlabel('y (m)');
ylabel('Temperature (°C)');
title('Temperature distribution on the outlet');
grid on;
%visualization of the centerline temperature
figure('Name', 'Visualization of the centerline temperature');
plot(x_coords, T_2D(ceil(ny/2),:), 'r-o'); 
%ylim([T_wall-1, T_wall+1]);        
xlabel('x (m)');
ylabel('Temperature (°C)');
title('Temperature distribution on the centerline');
grid on;
%}