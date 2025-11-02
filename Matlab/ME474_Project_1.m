clear;
clc;
close all;

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
u_x = @(y) 6*u_mean.*(y/H).*(1-(y/H));

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
xlabel('x (m)', 'FontSize', 15);       
ylabel('y (m)', 'FontSize', 15);        
%title(['Channel temperature field（nx=' num2str(nx) ', ny=' num2str(ny) '）']);
axis equal tight;

%%% 12.
figure('Name', 'Boundary conditions validation');

% Verify inlet temperature (plotting the temperature curve at x=0)
subplot(2, 2, 1);
plot(y_coords, T_2D(:, 1), 'b-o', 'MarkerFaceColor', 'blue');             
xlabel('y (m)');
ylabel('Temperature (°C)');
title('Inlet');
grid on;

% Verify visualizartion of the outlet temperature
subplot(2, 2, 2);
plot(x_coords, T_2D(2,:), 'b');  
hold on;
plot(x_coords, T_2D(ny-2,:), 'b');  
xlabel('x (m)');
ylabel('Temperature (°C)');
title('Streamwise Temperature');
grid on;

% Verify wall surface temperature (plotting the temperature curve at y=0)
subplot(2, 2, 3);
plot(x_coords, T_2D(ny, :), 'r-o');        
xlabel('x (m)');
ylabel('Temperature (°C)');
title('Wall (y = 0)');
grid on;

% Verify wall surface temperature (plotting the temperature curve at y=H)
subplot(2, 2, 4);
plot(x_coords, T_2D(1,:), 'r-o');       
xlabel('x (m)');
ylabel('Temperature (°C)');
title('Wall (y = H)');
grid on;

%%% 13.
figure();
plot(y_coords, T_2D(:, nx), 'b-o', 'MarkerFaceColor', 'blue', 'LineWidth', 2);             
xlabel('y (m)',  'FontSize', 20);
ylabel('Temperature (°C)', 'FontSize', 20);
%title('Outlet');
grid on;

figure();
plot(x_coords, T_2D(ceil(ny/2),:), 'r-o', 'LineWidth', 2);
hold on;
T_mean_weighted =@(T) u_x(y_coords) * T / sum(u_x(y_coords));
T_mean_weighted_values = zeros(size(x_coords, 2), 1);

for i=1:size(x_coords, 2)
    T_mean_weighted_values(i) = T_mean_weighted(T_2D(:,i));
end

plot(x_coords, T_mean_weighted_values, 'b-o', 'LineWidth', 2);
xlabel('x (m)', 'FontSize', 20);
ylabel('Temperature (°C)', 'FontSize', 20);
%title('Wall (y = H)');
grid on;

x_entry_index = 0;
for i=1:size(x_coords, 2)
    if T_2D(ceil(ny/2),i) >= 0.9 * T_wall
        x_entry_index = i;
        break
    end
end

xline(x_coords(x_entry_index), '-', 'LineWidth', 3);
legend('Centerline temperature', 'Weighted mean temperature', 'Entry length', 'FontSize', 20)
disp(x_coords(x_entry_index));

%%% 14.
function [error, residual, it, sol] = sor(A, b, w, guess, tol)
    n = size(b,1);
    sol = guess;
    solnew = zeros(n,1);

    error = [];
    residual = [];

    it = 0;
    var1 = tol + 1;
    var2 = tol + 1;
    while var1 > tol || var2 > tol
        for i=1:n
            solnew(i) = (b(i) - A(i, i+1:end)*sol(i+1:end) - A(i, 1:i-1)*solnew(1:i-1)) / A(i, i);
        end

        solold = sol;
        sol = (1 - w) * sol + w * solnew;
        var1 = sqrt(sum((sol - solold).^2)) / sqrt(sum((solold).^2));
        var2 = sqrt(sum((A*sol - b).^2)) / sqrt(sum((diag(A)' * sol).^2));

        error(end + 1) = var1;
        residual(end + 1) = var2;
        it = it + 1;
    end
end

guess = 50*ones(N, 1);

for j = 1:nx
    % Lower wall surface
    n = (ny-1)*nx + j;
    guess(n) = T_wall;
    
    % Upper wall surface
    n = j;
    guess(n) = T_wall;
end

[error1, residual1, it1, sol] = sor(A, b, 1, guess, 1e-5);
[error2, residual2, it2, sol] = sor(A, b, 1.5, guess, 1e-5);

figure();
iter1 = 1:1:it1;
iter2 = 1:1:it2;
plot(iter1, error1, 'LineWidth', 2);
hold on;
plot(iter2, error2, 'LineWidth', 2);
grid on;
xlabel('Iterations', 'FontSize', 20);
ylabel('Relative iteration error', 'FontSize', 20);
legend('Relative error (w=1)', 'Relative error (w=1.5)', 'FontSize', 20);

figure();
plot(iter1, residual1, 'LineWidth', 2);
hold on;
plot(iter2, residual2, 'LineWidth', 2);
grid on;
xlabel('Iterations', 'FontSize', 20);
ylabel('Normalized residuals', 'FontSize', 20);
legend('Normalized residual (w=1)', 'Normalized residual (w=1.5)', 'FontSize', 20);

fprintf('Number of iterations with Gauss-Seidel (w=1): %d\n', it1);
fprintf('Number of iterations with over_relaxation (w=1.5): %d\n', it2);