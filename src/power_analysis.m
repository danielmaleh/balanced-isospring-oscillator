%% POWER_ANALYSIS - Power Dissipation and Reserve Margin Analysis
%
% This script analyzes the power requirements and energy reserve margins
% for the Balanced Iso-Spring Oscillator (BISO) system, considering the
% gear train efficiency and mainspring (barrel) characteristics.
%
% The analysis evaluates how the number of gear stages and barrel turns
% affects the available power reserve margin for maintaining oscillation.
%
% Author: Daniel Abraham Elmaleh
% Project: 2D Balanced Iso-Spring Oscillator (BISO)
%
%% Variable Definitions
%
%   k_r         - Equivalent spring constant of the system [N/m]
%   m_r         - Equivalent reduced mass of the system [kg]
%   r           - Effective orbit radius [m]
%   f           - Natural oscillation frequency [Hz]
%   P           - Power dissipated in the system [W]
%   Q           - Quality factor (dimensionless)
%   r_min       - Minimum orbit radius [m]
%   r_max       - Maximum orbit radius [m]
%   w           - Angular velocity of oscillating mass [rad/s]
%   n           - Number of gear stages in the train
%   N           - Number of barrel (mainspring) turns
%   eta         - Efficiency of the gear train
%   H           - Power reserve (autonomy) [s]
%   C_bar       - Barrel torque constant [N*m]

%% Clear Environment

clear;
clc;
close all;

%% System Constants

% Orbit radius range
r_min = 3e-4;       % Minimum orbit radius [m]
r_max = 3e-3;       % Maximum orbit radius [m]

% Oscillator parameters
Q = 900;            % Quality factor
f = 9;              % Operating frequency [Hz]
w = 2 * pi * f;     % Angular velocity [rad/s]

% Mechanical properties
m_r = 0.13;         % Reduced mass [kg]
k_r = 700;          % Spring constant [N/m]

% Minimum required power reserve
H_min = 691200;     % [s] (approximately 8 days)

%% Parameter Ranges

NUM_POINTS = 8;

% Number of gear stages (4 to 11)
n = linspace(4, 11, NUM_POINTS);

% Number of barrel turns (2 to 9)
N = [2, 3, 4, 5, 6, 7, 8, 9];

% Orbit radius range
r = linspace(r_min, r_max, NUM_POINTS);

%% Compute Barrel Torque Constants

C_bar = zeros(1, NUM_POINTS);
for i = 2:NUM_POINTS
    C_bar(i) = get_barrel_torque_constant(i);
end

%% Gear Train Efficiency

% Efficiency decreases with number of stages (98% per stage)
eta = (0.98) .^ n;

%% Power Calculations

% Power required by barrel accounting for efficiency losses
P_barrel = C_bar * w .* (1 ./ eta);

% Power dissipated by the oscillator: P = (8 * pi^3 * f^3 * m_r * r^2) / Q
P_dissipated = (8 * (pi * f)^3 * m_r * (r.^2)) / Q;

% Power reserve (autonomy) calculation
% H = (N * Q * eta * C_bar) / (2 * pi^2 * f^2 * m_r * r^2)
H = (N * Q .* eta .* C_bar) ./ (2 * pi^2 * f^2 * m_r * r.^2);

% Reserve margin above minimum requirement
margin = H - H_min;

%% Visualization

% Plot 1: Margin vs Number of Gear Stages
figure('Name', 'Reserve Margin vs Gear Stages', 'NumberTitle', 'off');
hold on;
plot(n, margin, 'LineWidth', 1.5);
xlabel('Number of Gear Stages');
ylabel('Reserve Margin [s]');
title('Power Reserve Margin vs Gear Train Length');
axis tight;
grid on;
hold off;

% Plot 2: Margin vs Orbit Radius
figure('Name', 'Reserve Margin vs Radius', 'NumberTitle', 'off');
hold on;
plot(r, margin, 'LineWidth', 1.5);
xlabel('Orbit Radius [m]');
ylabel('Reserve Margin [s]');
title('Power Reserve Margin vs Oscillation Amplitude');
axis tight;
grid on;
hold off;

% Plot 3: Margin vs Number of Barrel Turns
figure('Name', 'Reserve Margin vs Barrel Turns', 'NumberTitle', 'off');
hold on;
plot(N, margin, 'LineWidth', 1.5);
xlabel('Number of Barrel Turns');
ylabel('Reserve Margin [s]');
title('Power Reserve Margin vs Mainspring Capacity');
axis tight;
grid on;
hold off;

%% Functions

function C_b = get_barrel_torque_constant(N)
    % GET_BARREL_TORQUE_CONSTANT Return barrel torque for given turns
    %
    %   C_b = get_barrel_torque_constant(N)
    %
    %   Returns the torque constant [N*m] of the mainspring barrel
    %   as a function of the number of available turns.
    %
    %   These values are empirically determined for the specific
    %   barrel design used in this oscillator.

    switch N
        case 2
            C_b = 107;
        case 3
            C_b = 105;
        case 4
            C_b = 102;
        case 5
            C_b = 95;
        case 6
            C_b = 90;
        case 7
            C_b = 80;
        case 8
            C_b = 70;
        case 9
            C_b = 60;
        otherwise
            C_b = 0;
            warning('Barrel turn count %d not defined, using C_b = 0', N);
    end
end
