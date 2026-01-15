%% REDUCED_MASS - Reduced Mass and Mass Isotropy Analysis
%
% This script computes the equivalent reduced mass of the Balanced
% Iso-Spring Oscillator (BISO) system by analyzing the kinetic energy
% of all moving masses in the four pantograph mechanisms.
%
% The reduced mass is derived from the total kinetic energy using:
%   E_kinetic = (1/2) * m_reduced * v^2
%
% Mass isotropy is evaluated by computing the center of mass error
% normalized to the orbit radius.
%
% Author: Daniel Abraham Elmaleh
% Project: 2D Balanced Iso-Spring Oscillator (BISO)
%
%% Variable Definitions
%
%   m               - Point mass at each pantograph joint [kg]
%   pos_mi          - Position of the i-th mass (i = 1,...,8) over time
%   vitesses        - Velocity matrix of all masses at 360 time steps
%   masse_reduite   - Equivalent reduced mass of the system [kg]
%   r               - Effective orbit radius [m]
%   alpha           - Internal angle at the center pivot of a rhombus [rad]
%   f               - Natural oscillation frequency [Hz]
%   rho_min         - Minimum orbit radius [m]
%   rho_max         - Maximum orbit radius [m]
%   w               - Angular velocity of oscillating mass [rad/s]
%   a, l, L         - Geometric distance parameters [m]

%% Clear Environment

clear;
clc;
close all;

%% Constants and Vector Definitions

% Orbit radius range
rho_min = 0.3e-3;   % [m]
rho_max = 0.4e-3;   % [m]

% Parameter ranges
NUM_POINTS = 15;
f = linspace(1, 15, NUM_POINTS);        % Frequency range [Hz]
r = linspace(rho_min, rho_max, NUM_POINTS);  % Radius range [m]

% Create meshgrid for surface plots
[R, F] = meshgrid(r, f);

%% Kinetic Energy and Reduced Mass Computation

% Initialize matrices
kinetic_energy_table = zeros(NUM_POINTS, NUM_POINTS);
reduced_mass         = zeros(NUM_POINTS, NUM_POINTS);
isotropy_error       = zeros(NUM_POINTS, NUM_POINTS);

% Compute kinetic energy and isotropy for each (frequency, radius) pair
for i = 1:NUM_POINTS
    for j = 1:NUM_POINTS
        kinetic_energy_table(i, j) = compute_kinetic_energy(f(i), r(j));
        isotropy_error(i, j) = compute_mass_isotropy_error(f(i), r(j));

        % Normalize isotropy error to radius (percentage)
        isotropy_error(i, j) = (isotropy_error(i, j) / r(j)) * 100;
    end
end

% Derive reduced mass from kinetic energy: E_c = (1/2) * m * v^2 = (1/2) * m * (w*r)^2
% Therefore: m = 2*E_c / (w*r)^2 = 2*E_c / (4*pi^2*f^2*r^2)
for i = 1:NUM_POINTS
    for j = 1:NUM_POINTS
        reduced_mass(i, j) = (2 * kinetic_energy_table(i, j)) / ...
                            (4 * pi^2 * f(i)^2 * r(j)^2);
    end
end

%% Visualization - Reduced Mass

% 3D Surface plot: Reduced Mass vs Frequency and Radius
figure('Name', 'Reduced Mass Surface', 'NumberTitle', 'off');
surfc(F, R, reduced_mass);
xlabel('Frequency [Hz]');
ylabel('Radius [m]');
zlabel('Reduced Mass [kg]');
title('Reduced Mass as Function of Frequency and Radius');
grid on;
axis tight;
shading interp;
colorbar;
colormap(jet(100));

% 2D Plot: Reduced Mass vs Radius (for different frequencies)
figure('Name', 'Reduced Mass vs Radius', 'NumberTitle', 'off');
hold on;
plot(r, reduced_mass);
legend_entries = cell(1, NUM_POINTS);
for i = 1:NUM_POINTS
    legend_entries{i} = sprintf('f = %d Hz', i);
end
legend(legend_entries, 'Location', 'best');
xlabel('Radius [m]');
ylabel('Reduced Mass [kg]');
title('Reduced Mass vs Orbit Radius');
axis tight;
grid on;
hold off;

% 2D Plot: Reduced Mass vs Frequency (for different radii)
figure('Name', 'Reduced Mass vs Frequency', 'NumberTitle', 'off');
hold on;
plot(f, reduced_mass);
legend_entries = cell(1, NUM_POINTS);
legend_entries{1} = 'r = 0.3 mm';
legend_entries{NUM_POINTS} = 'r = 0.4 mm';
for i = 2:NUM_POINTS-1
    legend_entries{i} = sprintf('r%d', i);
end
legend(legend_entries, 'Location', 'best');
xlabel('Frequency [Hz]');
ylabel('Reduced Mass [kg]');
title('Reduced Mass vs Oscillation Frequency');
axis tight;
grid on;
hold off;

%% Visualization - Mass Isotropy Error

% 2D Plot: Isotropy Error vs Radius
figure('Name', 'Isotropy Error vs Radius', 'NumberTitle', 'off');
hold on;
plot(r, isotropy_error);
legend_entries = cell(1, NUM_POINTS);
for i = 1:NUM_POINTS
    legend_entries{i} = sprintf('f = %d Hz', i);
end
legend(legend_entries, 'Location', 'best');
xlabel('Radius [m]');
ylabel('Mean Center of Mass Error [%]');
title('Mass Isotropy Error vs Orbit Radius');
axis tight;
grid on;
hold off;

% 2D Plot: Isotropy Error vs Frequency
figure('Name', 'Isotropy Error vs Frequency', 'NumberTitle', 'off');
hold on;
plot(f, isotropy_error);
legend_entries = cell(1, NUM_POINTS);
legend_entries{1} = 'r = 0.3 mm';
legend_entries{NUM_POINTS} = 'r = 0.4 mm';
for i = 2:NUM_POINTS-1
    legend_entries{i} = sprintf('r%d', i);
end
legend(legend_entries, 'Location', 'best');
xlabel('Frequency [Hz]');
ylabel('Mean Center of Mass Error [%]');
title('Mass Isotropy Error vs Oscillation Frequency');
axis tight;
grid on;
hold off;

%% Functions

function E_kinetic = compute_kinetic_energy(f, r)
    % COMPUTE_KINETIC_ENERGY Calculate total kinetic energy of all masses
    %
    %   E_kinetic = compute_kinetic_energy(f, r)
    %
    %   Computes the average kinetic energy of the 8 masses in the four
    %   pantograph mechanisms over one complete oscillation cycle.

    % Geometric constants [m]
    L = 25.213e-3;
    l = 4e-3;
    a = 15e-3;

    % Dynamic parameters
    w = 2 * pi * f;     % Angular velocity [rad/s]
    m = 17.11e-3;       % Mass at each joint [kg]

    % Time discretization: 360 steps per cycle
    NUM_TIME_STEPS = 360;
    NUM_MASSES = 8;

    % Initialize position and velocity matrices
    % Dimensions: (coordinate x/y, mass index 1-8, time step)
    mass_positions = zeros(2, NUM_MASSES, NUM_TIME_STEPS);
    mass_velocities = zeros(2, NUM_MASSES, NUM_TIME_STEPS);

    % Compute mass positions over one complete cycle
    for i = 1:NUM_TIME_STEPS
        t = (i * 2 * pi) / (NUM_TIME_STEPS * w);

        % Center position (oscillating mass)
        center_x = r * cos(w * t);
        center_y = r * sin(w * t);

        % Compute pantograph angle
        alpha_val = compute_pantograph_angle(r, t, w, L, l, a);

        % Pantograph 1 (South) - masses 1 and 2
        mass_positions(1, 1, i) = center_x - 2*a * trig_x(pi + alpha_val);
        mass_positions(2, 1, i) = center_y - l - 2*a * trig_y(pi + alpha_val);
        mass_positions(1, 2, i) = center_x + 2*a * trig_x(pi + alpha_val);
        mass_positions(2, 2, i) = center_y - l - 2*a * trig_y(pi + alpha_val);

        % Pantograph 2 (West) - masses 3 and 4
        mass_positions(1, 3, i) = center_x - 2*a * trig_x(pi/2 + alpha_val);
        mass_positions(2, 3, i) = center_y - l - 2*a * trig_y(pi/2 + alpha_val);
        mass_positions(1, 4, i) = center_x + 2*a * trig_x(pi/2 + alpha_val);
        mass_positions(2, 4, i) = center_y - l - 2*a * trig_y(pi/2 + alpha_val);

        % Pantograph 3 (North) - masses 5 and 6
        mass_positions(1, 5, i) = center_x - 2*a * trig_x(-alpha_val);
        mass_positions(2, 5, i) = center_y - l - 2*a * trig_y(-alpha_val);
        mass_positions(1, 6, i) = center_x + 2*a * trig_x(-alpha_val);
        mass_positions(2, 6, i) = center_y - l - 2*a * trig_y(-alpha_val);

        % Pantograph 4 (East) - masses 7 and 8
        mass_positions(1, 7, i) = center_x - 2*a * trig_x(alpha_val);
        mass_positions(2, 7, i) = center_y - l - 2*a * trig_y(alpha_val);
        mass_positions(1, 8, i) = center_x + 2*a * trig_x(alpha_val);
        mass_positions(2, 8, i) = center_y - l - 2*a * trig_y(alpha_val);
    end

    % Compute velocities using finite differences
    time_step = (2 * pi) / (NUM_TIME_STEPS * w);

    for i = 1:NUM_TIME_STEPS-1
        for k = 1:NUM_MASSES
            mass_velocities(1, k, i) = (mass_positions(1, k, i+1) - mass_positions(1, k, i)) / time_step;
            mass_velocities(2, k, i) = (mass_positions(2, k, i+1) - mass_positions(2, k, i)) / time_step;
        end
    end

    % Compute total kinetic energy
    velocity_squared = zeros(NUM_MASSES, NUM_TIME_STEPS);
    total_kinetic_energy = 0;

    for i = 1:NUM_TIME_STEPS
        for k = 1:NUM_MASSES
            velocity_squared(k, i) = mass_velocities(1, k, i)^2 + mass_velocities(2, k, i)^2;
            total_kinetic_energy = total_kinetic_energy + 0.5 * m * velocity_squared(k, i);
        end
    end

    % Return average kinetic energy over the cycle
    E_kinetic = total_kinetic_energy / NUM_TIME_STEPS;
end


function mean_error = compute_mass_isotropy_error(f, r)
    % COMPUTE_MASS_ISOTROPY_ERROR Calculate center of mass deviation
    %
    %   mean_error = compute_mass_isotropy_error(f, r)
    %
    %   Evaluates the isotropy of the mass distribution by computing
    %   the mean deviation of the center of mass from the origin.

    % Geometric constants [m]
    L = 25.213e-3;
    l = 4e-3;
    a = 15e-3;
    w = 2 * pi * f;

    NUM_TIME_STEPS = 360;
    NUM_MASSES = 8;

    % Initialize position matrix
    mass_positions = zeros(2, NUM_MASSES, NUM_TIME_STEPS);
    mean_error = 0;

    for t_idx = 1:NUM_TIME_STEPS-1
        t = (t_idx * 2 * pi) / (NUM_TIME_STEPS * w);

        center_x = r * cos(w * t);
        center_y = r * sin(w * t);

        alpha_val = compute_pantograph_angle(r, t, w, L, l, a);

        % Pantograph 1
        mass_positions(1, 1, t_idx) = center_x - 2*a * trig_x(pi + alpha_val);
        mass_positions(2, 1, t_idx) = center_y - l - 2*a * trig_y(pi + alpha_val);
        mass_positions(1, 2, t_idx) = center_x + 2*a * trig_x(pi + alpha_val);
        mass_positions(2, 2, t_idx) = center_y - l - 2*a * trig_y(pi + alpha_val);

        % Pantograph 2
        mass_positions(1, 3, t_idx) = center_x - 2*a * trig_x(pi/2 + alpha_val);
        mass_positions(2, 3, t_idx) = center_y - l - 2*a * trig_y(pi/2 + alpha_val);
        mass_positions(1, 4, t_idx) = center_x + 2*a * trig_x(pi/2 + alpha_val);
        mass_positions(2, 4, t_idx) = center_y - l - 2*a * trig_y(pi/2 + alpha_val);

        % Pantograph 3
        mass_positions(1, 5, t_idx) = center_x - 2*a * trig_x(-alpha_val);
        mass_positions(2, 5, t_idx) = center_y - l - 2*a * trig_y(-alpha_val);
        mass_positions(1, 6, t_idx) = center_x + 2*a * trig_x(-alpha_val);
        mass_positions(2, 6, t_idx) = center_y - l - 2*a * trig_y(-alpha_val);

        % Pantograph 4
        mass_positions(1, 7, t_idx) = center_x - 2*a * trig_x(alpha_val);
        mass_positions(2, 7, t_idx) = center_y - l - 2*a * trig_y(alpha_val);
        mass_positions(1, 8, t_idx) = center_x + 2*a * trig_x(alpha_val);
        mass_positions(2, 8, t_idx) = center_y - l - 2*a * trig_y(alpha_val);

        % Compute mean distance from center for all masses at this time step
        mean_distance_at_t = 0;
        for i = 1:NUM_MASSES
            mean_distance_at_t = mean_distance_at_t + ...
                sqrt(mass_positions(1, i, t_idx)^2 + mass_positions(2, i, t_idx)^2);
        end
        mean_error = mean_error + mean_distance_at_t / NUM_MASSES;
    end

    mean_error = mean_error / NUM_TIME_STEPS;
end


function result = trig_x(alpha)
    % TRIG_X Compute x-component trigonometric factor
    result = cos(alpha / 2);
end


function result = trig_y(alpha)
    % TRIG_Y Compute y-component trigonometric factor
    result = sin(alpha / 2);
end


function alpha = compute_pantograph_angle(r, t, w, L, l, a)
    % COMPUTE_PANTOGRAPH_ANGLE Calculate internal pantograph angle
    %
    %   alpha = compute_pantograph_angle(r, t, w, L, l, a)
    %
    %   Computes the internal angle of the pantograph rhombus based on
    %   the current position of the oscillating center mass.

    numerator = sqrt(r^2 * cos(w*t)^2 + (r * sin(w*t) + L - l)^2);
    denominator = 2 * a;
    alpha = pi - 2 * asin(numerator / denominator);
end
