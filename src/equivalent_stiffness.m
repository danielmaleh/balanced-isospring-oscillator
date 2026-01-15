%% EQUIVALENT_STIFFNESS - Equivalent Spring Stiffness Analysis
%
% This script computes the equivalent linear stiffness (K_eq) of the
% Balanced Iso-Spring Oscillator (BISO) system as a function of the
% oscillating mass position in polar coordinates.
%
% The system consists of four pantograph mechanisms arranged symmetrically
% (North, South, East, West) that create an isotropic spring behavior.
%
% Author: Daniel Abraham Elmaleh
% Project: 2D Balanced Iso-Spring Oscillator (BISO)
%
%% Variable Definitions
%
%   Zp          - Equivalent spring stiffness matrix k_eq = f(rho, theta)
%   rho         - Effective orbit radius [m]
%   theta       - Angle of reduced mass in polar coordinates [rad]
%   rho_min     - Minimum orbit radius [m]
%   rho_max     - Maximum orbit radius [m]
%   k_beta      - Spring constant at beta pivot point [N*m/rad]
%   k_gamma     - Spring constant at gamma pivot points [N*m/rad]
%   a           - Pantograph arm length [m]
%   l           - Half-length of central platform [m]
%   L           - Reference length parameter [m]
%   d0, d1      - Distance variables during motion [m]

%% Clear Environment

clear;
clc;
close all;

%% Constants, Vectors & Matrices Generation

% Orbit radius range
rho_min = 0.3e-3;   % [m]
rho_max = 0.4e-3;   % [m]

% Generate coordinate vectors
rho   = linspace(rho_min, rho_max);
theta = linspace(0, 2*pi);

% Compute equivalent stiffness matrix
Zp = compute_equivalent_stiffness(rho, theta);

%% Visualization

figure('Name', 'Equivalent Stiffness K_eq = f(x, y)', 'NumberTitle', 'off');

polarplot3d(Zp, 'RadialRange', [rho_min, rho_max]);

xlabel('Position x [m]');
ylabel('Position y [m]');
zlabel('K_{eq} (Linear Stiffness) [N/m]');

grid on;
axis tight;
shading interp;
colorbar;
colormap('turbo');

title('Equivalent Spring Stiffness Distribution');

%% Functions

function Zp = compute_equivalent_stiffness(rho, theta)
    % COMPUTE_EQUIVALENT_STIFFNESS Calculate stiffness over position grid
    %
    %   Zp = compute_equivalent_stiffness(rho, theta)
    %
    %   Computes the equivalent spring stiffness for all four pantograph
    %   mechanisms (North, South, East, West) at each position (rho, theta).

    % Geometric constants [m]
    L = 25.213e-3;  % Reference length
    l = 4e-3;       % Half platform length
    a = 15e-3;      % Pantograph arm length

    % Spring constants [N*m/rad]
    k_gamma = 7.85e-3;
    k_beta  = 6.44e-3;

    % Initial configuration (equilibrium position)
    d00    = L - l;
    delta0 = asin((d00 / 2) / a);
    beta0  = 2 * delta0;
    gamma0 = delta0;

    % Initialize stiffness matrix
    n_rho   = numel(rho);
    n_theta = numel(theta);
    Zp = zeros(n_rho, n_theta);

    % Compute stiffness at each position
    for r_idx = 1:n_rho
        for t_idx = 1:n_theta

            % Convert polar to Cartesian coordinates
            x = rho(r_idx) * cos(theta(t_idx));
            y = rho(r_idx) * sin(theta(t_idx));

            % East pantograph contribution
            d0 = L - (x + l);
            d1 = y;
            Zp(r_idx, t_idx) = compute_stiffness_contribution(...
                a, d0, d1, k_gamma, k_beta, gamma0, beta0, ...
                r_idx, t_idx, Zp, rho);

            % West pantograph contribution
            d0 = L + (x - l);
            d1 = y;
            Zp(r_idx, t_idx) = compute_stiffness_contribution(...
                a, d0, d1, k_gamma, k_beta, gamma0, beta0, ...
                r_idx, t_idx, Zp, rho);

            % North pantograph contribution
            d0 = L - (y + l);
            d1 = x;
            Zp(r_idx, t_idx) = compute_stiffness_contribution(...
                a, d0, d1, k_gamma, k_beta, gamma0, beta0, ...
                r_idx, t_idx, Zp, rho);

            % South pantograph contribution
            d0 = L + (y - l);
            d1 = x;
            Zp(r_idx, t_idx) = compute_stiffness_contribution(...
                a, d0, d1, k_gamma, k_beta, gamma0, beta0, ...
                r_idx, t_idx, Zp, rho);
        end
    end
end


function k = compute_stiffness_contribution(a, d0, d1, k_gamma, k_beta, ...
                                            gamma0, beta0, r_idx, t_idx, Zp, rho)
    % COMPUTE_STIFFNESS_CONTRIBUTION Calculate single pantograph stiffness
    %
    %   Computes the elastic potential energy derivative to obtain the
    %   equivalent linear stiffness contribution from one pantograph.

    % Compute distance and angles
    d = sqrt(d0^2 + d1^2);

    phi   = acos(d0 / d);
    delta = asin((d / 2) / a);
    beta  = 2 * delta;

    % Gamma angles for the two pivot points
    gamma1 = delta - phi;
    gamma2 = delta + phi;

    % Elastic potential energy contributions
    e_gamma = 2 * k_gamma * ((gamma1 - gamma0)^2 + (gamma2 - gamma0)^2);
    e_beta  = 2 * k_beta * (beta - beta0)^2;

    % Accumulated stiffness (second derivative of potential energy)
    k = Zp(r_idx, t_idx) + (1 / rho(r_idx)^2) * (e_gamma + e_beta);
end
