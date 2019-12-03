%{
    Title: Single Reaction Wheel Control Model
    Author: James Parkus
    Date: 12/2/19
    Purpose: Provide a reaction for a given wheel spin for the S-CUBED
    Cubesat preliminary ADACS system model.

    Input: Desired angular velocity of body in RPM.
    Output: Required flywheel angular velocity.
%}


function [w_flywheel] = reaction(w_body)

I_c = 0.05;
I_w = 2*10^-3;

w_flywheel = -w_body*(1 + I_c/I_w);

end