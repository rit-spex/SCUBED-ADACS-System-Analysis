%{
    Title: SCUBED ADACS Analysis
    Author: James E. Parkus, Amber Dubill
    Date: 10/29/2019
    Purpose: This script will calculate the characteristics of the SCUBED
    satellite's momentum wheels. The inputs are ____________, and the
    outputs are the required angular velocity of each wheels and the
    associated torques.

    Nomenclature:
    Earth Heliocentric Orbit - EHO

    Coordinate System:
        z - axis is the longitudinal axis of the CubeSat (axis that is always pointing towards the sun)
        y - axis is tangent to the orbital path
        x - axis completes the set

    Momentum Wheel Orientation:
        Wheel 1 -> principal axis is x-axis
        Wheel 2 -> principal axis is y-axis
        Wheel 3 -> principal axis is z-axis
%}

clc
clear
close all
format long

%% Memory Allocation
n = 10000;

%% Constants
T_EHO = 365.256363004*86400; % [s] - Orbital Period of EHO
orbit_altitude = 550; % [km]
earth_radius = 6378; % [km]
mu_earth = 398601.2; % [km^3/s^2] - Gravitational parameter for Earth
T_earth = 2*pi*(mu_earth)^(-1/2)*(orbit_altitude + earth_radius)^(3/2);

%% Simulink Simulation
% Initial Conditions
t = 2*T_earth; % simulation runtime
dt = 0.001;

Mgx = 10^-11; % [Nm] - Solar pressure
Mgy = 10^-11; % [Nm] - Solar pressure
Mgz = 10^-11; % [Nm] - Solar pressure
M_SRP = 10^-4; % [Nm] - Solar Sailing Pressure

I = 2*10^-3; % [kg m^2] - Spin moment of inertia
J = 10^-3;
A = 0.032;
B = 0.021;
C = 0.046;

J_I = A + I + 2*J;
J_II = B + I + 2*J;
J_III = C + I + 2*J;

% Wheel Initial Angular Velocities
omega1_i = 0; % Initial conditions of wheel 1 at t = 0
omega2_i = 0; % Initial conditions of wheel 2 at t = 0
omega3_i = 0; % Initial conditions of wheel 3 at t = 0

% Body Angular Velocity
omega_x = 2*pi/T_EHO;
omega_y = 2*pi/T_earth;
omega_z = Mgz/I;
alpha_x = 0;
alpha_y = 0;
alpha_z = 0;

% Run Simulation
Simulation = sim('Momentum_Wheel_Model');

% Extract results
tout = Simulation.omega_1.time;
omega_1 = Simulation.omega_1.signals.values;
omega_2 = Simulation.omega_2.signals.values;
omega_3 = Simulation.omega_3.signals.values;

% Analyse results
omega_magnitude = sqrt(omega_1.^2 + omega_2.^2 + omega_3.^2);

%% Plotting
rpm_conversion = 60/(2*pi);

hold on
plot(tout,omega_1.*rpm_conversion);
plot(tout,omega_2.*rpm_conversion);
plot(tout,omega_3.*rpm_conversion);
plot(tout,omega_magnitude.*rpm_conversion,'--');
hold off
grid on
xlabel('Time');
ylabel('\omega');
legend('\omega_1','\omega_2', '\omega_3','\omega_{RMS}');
title('Outputs vs. Time');
xlim([0 max(tout)]);

