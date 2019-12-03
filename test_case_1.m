%{
    Title: Test Case #1 for SCUBED-ADACS model
    Author: James Parkus
    Date: 12/03/19
    Purpose: Test the trend of hand-calcs for Case #1 of the SCUBED-ADACS
    model. 
%}
clc
clear

%% Constants
T_EHO = 365.256363004*86400; % [s]

Mx = 10^-11;
My = Mx;
Mz = Mx;
M_SRP = 10^-4;

I = 2*10^-3;

omega_x = 2*pi/T_EHO;
i = 1;
for t = 0:10000:3.157*10^7

    omega_1(i) = Mx/I*t;
    omega_2(i) = Mz/(I*omega_x)*(1 - cos(omega_x*t)) - Mx/(I*omega_x)*sin(omega_x*t);
    omega_3(i) = Mx/(I*omega_x)*(1 - cos(omega_x*t)) + M_SRP*t/I + Mz/(I*omega_x)*sin(omega_x*t);
    time(i) = t;
    i = i + 1;
    
end

%% Plotting

figure()
hold on
plot(time,omega_1*(60/(2*pi)));
plot(time,omega_2*(60/(2*pi)));
plot(time,omega_3*(60/(2*pi)));
title('Case #1 Angular Velocity');
xlabel('Time [s]');
ylabel('Angular Velocity [rpm]');
grid on
hold off













