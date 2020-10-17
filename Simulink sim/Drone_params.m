% Motor Time Constants
Tau = 0.1;
Tau_1 = Tau;
Tau_2 = Tau;
Tau_3 = Tau;
Tau_4 = Tau;

% Steady State Gain of Motor+Prop Thrust
% Thrust / PWM
KM_GENERAL = 0.00509;
KM_1 = KM_GENERAL;
KM_2 = KM_GENERAL;
KM_3 = KM_GENERAL;
KM_4 = KM_GENERAL;

% Steady State Gain of Motor+Prop Torque
% Torque / PWM
KM_T_GENERAL = 0.000005;
KM_T_1 = KM_T_GENERAL;
KM_T_2 = KM_T_GENERAL;
KM_T_3 = KM_T_GENERAL;
KM_T_4 = KM_T_GENERAL;

% Inertia values
Ixx = 0.002100;
Iyy = 0.002417;
Izz = 0.004633;

% Mass of Drone
m = 0.475;

% PWM cutoff, usually 1000
PWM_CUT = 1000;

% Prop to Prop axial length
X_length = 0.155;
Y_length = 0.200;
