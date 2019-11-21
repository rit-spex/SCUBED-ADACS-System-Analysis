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
% p = 1.49558*10^11; % [m] - Semilatus Rectum of Earth Heliocentric Orbit (EHO)
% a = 149.6*10^9; % [m] - Semimajor axis of EHO
% u = 1.3271544*10^20; % [m] - Heliocentric Gravitational Constant
% e = 0.0167086; % [-] - Eccentricity of EHO
% avg_r = 1.4959965*10^11; % [m] - Mean Distance from Earth to Sun (= 1 AU)
T_EHO = 365.256363004*86400; % [s] - Orbital Period of EHO
orbit_altitude = 550; % [km]
earth_radius = 6378; % [km]
mu_earth = 398601.2; % [km^3/s^2] - Gravitational parameter for Earth
T_earth = 2*pi*(mu_earth)^(-1/2)*(orbit_altitude + earth_radius)^(3/2);
% theta = (linspace(0,2*pi,n))'; % [rad] - True Anomaly of EHO

%% Simulink Simulation
% Initial Conditions
t = 2*T_earth; % simulation runtime
dt = 0.001;

Mgx = 10^-11; % [Nm] - Solar pressure
Mgy = 10^-11; % [Nm] - Solar pressure
Mgz = 10^-4; % [Nm] - Solar sailing
I = 2*10^-3; % [kg m^2] - Spin moment of inertia
J = 10^-3;
A = 0.032; 
B = 0.021;
C = 0.046;
MGX = Mgx/I;
J_I = (B-C)/I;
C_I = 1 + A/2 + 2*J/I;
MGY = Mgy/I;
J_II = (C-A)/I;
C_II = 1 + B/2 + 2*J/I;
MGZ = Mgz/I;
J_III = (A-B)/I;
C_III = 1 + C/2 + 2*J/I;

% Wheel Initial Angular Velocities
omega1_i = 0; % Initial conditions of wheel 1 at t = 0
omega2_i = 0; % Initial conditions of wheel 2 at t = 0
omega3_i = 0; % Initial conditions of wheel 3 at t = 0

% Body Angular Velocity
omega_x = MGX;
omega_y = MGY;
omega_z = MGY;
alpha_x = MGX;
alpha_y = MGY;
alpha_z = MGY;

% Run Simulation
Simulation = sim('Nonlinear_momentum_wheel_model');

% Extract results
tout = Simulation.omega_1.time;
omega_1 = Simulation.omega_1.signals.values;
omega_2 = Simulation.omega_2.signals.values;
omega_3 = Simulation.omega_3.signals.values;

% Analyse results
omega_magnitude = sqrt(omega_1.^2 + omega_2.^2 + omega_3.^2);
% 
% k = 2;
% for i = T_earth:T_earth:length(tout)
%     upper = floor(k*T_earth);
%     lower = floor((k-1)*T_earth);
%     max_omega_1(k,1) = max(omega_1(lower:upper,1));
%     max_omega_2(k,1) = max(omega_2(lower:upper,1));
%     max_omega_3(k,1) = max(omega_3(lower:upper,1));
%     max_time(k,1) = tout(upper,1);
%     k = k + 1;
% end

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

% hold on
% plot(max_time,max_omega_1.*rpm_conversion);
% plot(max_time,max_omega_2.*rpm_conversion);
% plot(max_time,max_omega_3.*rpm_conversion);
% hold off
% grid on
% xlabel('Time');
% ylabel('\omega');
% legend('\omega_1','\omega_2', '\omega_3');
% title('Outputs vs. Time');
% xlim([0 max(max_time)]);

% %% Analysis
% omega_0 = 2*pi/T; % [rad/s] - Required Angular Velocity to maintain sun-normal orientation throughout EHO
% Mgx = 10^-11; % [Nm] - Estimated magnitude of solar radiation pressure imposed S/C torque
% Mgy = Mgx;
% Mgz = 10^-8; % [Nm] - Estimated magnitude of SRP imposed torque due to diffractive sail
% I = 10^-3; % [kg m^2] - Spin-axis moment of inertia for momentum wheels
% 
% syms omega_1(t) omega_2(t) omega_3(t)
% 
% omega_1(t) = Mgx/I*t;                                                                  % [rad/s] - Angular Velocity of Momentum Wheel 1
% omega_2(t) = Mgy/(I*omega_0)*sin(omega_0*t) + Mgz/(I*omega_0)*(1 - cos(omega_0*t));    % [rad/s] - Angular Velocity of Momentum Wheel 2
% omega_3(t) = Mgy/(I*omega_0)*(cos(omega_0*t) - 1) + Mgz/(I*omega_0)*sin(omega_0*t);    % [rad/s] - Angular Velocity of Momentum Wheel 3
% 
% i = 1;
% for t = 0:T/1000:T
%     Omega(i,1:4) = [omega_1(t) omega_2(t) omega_3(t) t/86400];
%     Omega_Norm(i,1) = norm(Omega(i,1:3));
%     i = i + 1;
% end
% Omega_rpm(:,1:3) = Omega(:,1:3)*180/(pi);
% 
% max_W1 = max(Omega_rpm(:,1));
% max_W2 = max(Omega_rpm(:,2));
% max_W3 = max(Omega_rpm(:,3));
% 
% if max_W1 > max_W2 && max_W1 > max_W3
%     wheel = 1;
%     max_W = max_W1;
% end
% if max_W2 > max_W1 && max_W2 > max_W3
%     wheel = 2;
%     max_W = max_W2;
% end
% if max_W3 > max_W2 && max_W3 > max_W1
%     wheel = 3;
%     max_W = max_W3;
% end
% 
% fprintf('The maximum spin of reaction wheel %1.1g is %5.3f rpm.\n\n',wheel,max_W);
% 
% %% Plotting
% fig1 = figure();
% hold on
% plot(Omega(:,4),Omega_rpm(:,1:3));
% % plot(Omega(:,4),Omega_Norm,'--');
% ylabel('\omega [rad/s]');
% xlabel('Time [days]');
% title('Angular Velocity of Momentum Wheels through One Solar Orbit');
% legend('Wheel 1','Wheel 2','Wheel 3');
% xlim([0 max(t/86400)]);
% grid on
% hold off
% 
% % fig2 = figure();
% % hold on
% % plot(Omega(:,4),Omega_Norm);
% % ylabel('\omega [rad/s]');
% % xlabel('Time [s]');
% % title('Magnitude of Angular Velocity of Momentum Wheels through One Solar Orbit');
% % legend('Wheel 1','Wheel 2','Wheel 3');
% % xlim([0 max(t)]);
% % grid on
% % hold off
% 
% %}
% 
% 
