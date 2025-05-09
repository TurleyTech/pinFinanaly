% fin_temperature_simulation.m
% Simulates the temperature distribution along a brass fin using the shooting method
% with Runge-Kutta 5th order and Euler methods for solving the ODE

% Parameters and constants
L = 0.02;         % Fin length (m)
d = 0.001;        % Fin diameter (m)
k = 110;          % Thermal conductivity of brass (W/mK)
h = 20;           % Convective heat transfer coeff (W/m^2K)
T_ambient = 303.15;   % Ambient temperature (K) [30°C]
T_base = 348.15;      % Base temperature (K) [75°C]
P = pi * d;       % Perimeter (m)
A_c = pi * d^2 / 4;   % Cross-sectional area (m^2)
m2 = h * P / (k * A_c);  % m^2 parameter (1/m^2)

% Define ODE
function dydx = fin_ode(x, y)
    dydx = [y(2); m2 * y(1)];
end

% Euler method
function [y1_end, y2_end] = integrate_euler(slope_guess, h_step)
    global L T_base T_ambient
    N = ceil(L / h_step);
    y1 = T_base - T_ambient;
    y2 = slope_guess;
    x = 0;
    for i = 1:N
        if x + h_step > L
            h_step = L - x;
        end
        dy = fin_ode(x, [y1; y2]);
        y1 = y1 + h_step * dy(1);
        y2 = y2 + h_step * dy(2);
        x  = x + h_step;
    end
    y1_end = y1;
    y2_end = y2;
end

% Runge-Kutta 5th-order method
function [y1_end, y2_end] = integrate_rk5(slope_guess, h_step)
    global L T_base T_ambient
    N = ceil(L / h_step);
    y1 = T_base - T_ambient;
    y2 = slope_guess;
    x = 0;
    for i = 1:N
        if x + h_step > L
            h = L - x;
        else
            h = h_step;
        end
        y = [y1; y2];
        k1 = fin_ode(x, y);
        k2 = fin_ode(x + h/5, y + h/5 * k1);
        k3 = fin_ode(x + 3*h/10, y + h * (3/40*k1 + 9/40*k2));
        k4 = fin_ode(x + 3*h/5, y + h * (3/10*k1 - 9/10*k2 + 6/5*k3));
        k5 = fin_ode(x + h, y + h * (-11/54*k1 + 5/2*k2 - 70/27*k3 + 35/27*k4));
        k6 = fin_ode(x + 7*h/8, y + h * (1631/55296*k1 + 175/512*k2 + ...
            575/13824*k3 + 44275/110592*k4 + 253/4096*k5));
        y = y + h * (37/378*k1 + 250/621*k3 + 125/594*k4 + 512/1771*k6);
        x = x + h;
        y1 = y(1);
        y2 = y(2);
    end
    y1_end = y1;
    y2_end = y2;
end

% Shooting method using secant iteration
function [s_solution, iter] = shooting_secant(h_step, method)
    if strcmpi(method, 'Euler')
        integrator = @integrate_euler;
    else
        integrator = @integrate_rk5;
    end
    s1 = -100;
    s2 = -1000;
    [~, y2_L1] = integrator(s1, h_step);
    [~, y2_L2] = integrator(s2, h_step);
    iter = 0;
    while abs(y2_L2) > 1e-6 && iter < 50
        s_new = s2 - y2_L2 * (s2 - s1) / (y2_L2 - y2_L1);
        s1 = s2;
        s2 = s_new;
        y2_L1 = y2_L2;
        [~, y2_L2] = integrator(s2, h_step);
        iter = iter + 1;
    end
    s_solution = s2;
end

% Main execution
global L T_ambient T_base
h = 1e-4;
[slope_sol, iterations] = shooting_secant(h, 'RK5');
fprintf('Converged initial slope = %.3f K/m, in %d iterations.\n', slope_sol, iterations);

% Calculate temperature distribution
N = 200;
x_vals = linspace(0, L, N);
theta_vals = zeros(size(x_vals));
theta_vals(1) = T_base - T_ambient;
y2 = slope_sol;
x = 0;
for j = 2:N
    h_seg = x_vals(j) - x;
    y = [theta_vals(j-1); y2];
    k1 = fin_ode(x, y);
    k2 = fin_ode(x + h_seg/5, y + h_seg/5 * k1);
    k3 = fin_ode(x + 3*h_seg/10, y + h_seg * (3/40*k1 + 9/40*k2));
    k4 = fin_ode(x + 3*h_seg/5, y + h_seg * (3/10*k1 - 9/10*k2 + 6/5*k3));
    k5 = fin_ode(x + h_seg, y + h_seg * (-11/54*k1 + 5/2*k2 - 70/27*k3 + 35/27*k4));
    k6 = fin_ode(x + 7*h_seg/8, y + h_seg * (1631/55296*k1 + 175/512*k2 + ...
        575/13824*k3 + 44275/110592*k4 + 253/4096*k5));
    y_next = y + h_seg * (37/378*k1 + 250/621*k3 + 125/594*k4 + 512/1771*k6);
    theta_vals(j) = y_next(1);
    y2 = y_next(2);
    x = x_vals(j);
end
T_vals = theta_vals + T_ambient;

% Optional: Plot
figure;
plot(x_vals, T_vals, 'b-', 'LineWidth', 2);
xlabel('Fin Length (m)');
ylabel('Temperature (K)');
title('Temperature Distribution Along the Fin');
grid on;
