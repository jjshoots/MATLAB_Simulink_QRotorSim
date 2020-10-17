function [xn, uin, uctl, Itg] = Qrotor_inner_sim( xn, xc, ref, Itg, Ts, wnd )
%%========================================================================%
% Coded by Phang Swee King                                                %
% Maintained by Phang Swee King (swee_king@hotmail.com)                   %
% Copyright 2018 Taylor's Robotic Center, Taylor's University             %
% Updated: 9 Aug 2018           Version 1.0.3                             %
%%========================================================================%
% Single Step Simulation of Quadrotor with Inner Loop Controller Only!
% Outputs:
% xn  = [pos; v_b; att; omg] ---------------- updated actual system state
% uin = [u1; u2; u3; u4] -------------------- current input to the motors
% uctl= [ail; ele; rud; thr] ---------------- current input calculated by control
% Itg = [Itg1; Itg2; Itg3] ------------------ integration terms saved for next loop
%
% Inputs
% xn  = [pos; v_b; att; omg] ---------------- current actual system state
%        pos = [x; y; z] -------------------- current ground-axis position, m
%        vel = [u; v; w] -------------------- current body-axis velocity, m/s
%        att = [phi; tht; psi] -------------- current attitude angles, rad
%        omg = [p; q; r] -------------------- current angular rates, rad/s
% xc  = [pos; v_b; att; omg] ---------------- current measured system state for control 
%        pos = [x; y; z] -------------------- current ground-axis position, m
%        vel = [u; v; w] -------------------- current body-axis velocity, m/s
%        att = [phi; tht; psi] -------------- current attitude angles, rad
%        omg = [p; q; r] -------------------- current angular rates, rad/s
% ref = [phi_r; tht_r; psi_r] --------------- current reference input (angle reference), rad
% Itg = [Itg1; Itg2; Itg3] ------------------ integration terms saved from last loop
% Ts ---------------------------------------- sampling interval of control, sec
% wnd = [wnd_p; wnd_q; wnd_r] --------------- disturbance at 3 rotating axes, rad/s^2
%

% controller parameter, Ki & Kp for pitch/roll, Kir & Kpr for yaw
global K_P K_I K_Ir K_Pr iTs SSI_or_PI;
global Ixx Iyy Izz Ixm Iym Izm;

% Model parameters (these are values use in control)
g   = 9.781;        % gravitational acceleration, m/s^2
m   = 4.027;        % total mass of helicopter, kg

Kt1 = 0.00012551;       % Thrust coefficient, counter clock wise - Kt_ccw
Kt2 = 0.00013137;       % Thrust coefficient, clock wise         - Kt_cw
Kq  = 0.0000033;        % Torque coefficient                     - average Kq_cw, Kq_ccw

Kw1 = [566.6193 57.4600];     % the gain of PWM to rotation speed, counter clock wise. - pwm2omg_1st_ccw
Kw2 = [557.1556 57.9755];     % the gain of PWM to rotation speed, clock wise          - pwm2omg_1st_cw

Lm  = 0.303;        % the length of the arm of the platform, m

%% Inputs
% ref 1/2/3 is the phi/tht/psi reference
phi_r = ref(1);
tht_r = ref(2);
psi_r = ref(3);

psi_r = wrapToPi(psi_r); % convert to [-pi pi]


%% Single Control Step, Multiple Simulation Steps
ts = iTs;     % step size for inner-loop, sec
Ts = max( ts, Ts );

%% current states
phi = wrapToPi(xc(7));  % convert to [-pi pi]
tht = wrapToPi(xc(8));
psi = wrapToPi(xc(9));
omg = xc(10:12);

R = eye(3,3); % transformation from body coordinates to ground coordinates
R(1,1) = cos(psi)*cos(tht);
R(1,2) = cos(psi)*sin(tht)*sin(phi) - sin(psi)*cos(phi);
R(1,3) = cos(psi)*sin(tht)*cos(phi) + sin(psi)*sin(phi);
R(2,1) = sin(psi)*cos(tht);
R(2,2) = sin(psi)*sin(tht)*sin(phi) + cos(psi)*cos(phi);
R(2,3) = sin(psi)*sin(tht)*cos(phi) - cos(psi)*sin(phi);
R(3,1) = -sin(tht);
R(3,2) = cos(tht)*sin(phi);
R(3,3) = cos(tht)*cos(phi);

%% control parameter calculation for Jet Controller (integral and derivative)

e_phi = phi_r - phi;
e_tht = tht_r - tht;
e_psi = psi_r - psi;

e_omgp = e_phi - omg(1);
e_omgq = e_tht - omg(2);
e_omgr = e_psi - omg(3);

ie_omgp = min(max(Itg(1) + Ts * e_omgp, -1), 1);
ie_omgq = min(max(Itg(2) + Ts * e_omgq, -1), 1);
ie_omgr = min(max(Itg(3) + Ts * e_omgr, -1), 1);

Itg(1:3) = [ie_omgp; ie_omgq; ie_omgr];


%% For inner loop control
if SSI_or_PI == 1
    u_p = K_P * e_omgp + K_I * K_P * ie_omgp - K_I * Ixm * omg(1);
    u_q = K_P * e_omgq + K_I * K_P * ie_omgq - K_I * Iym * omg(2);
    u_r = K_Pr * e_omgr + K_Ir * K_Pr * ie_omgr - K_Ir * Izm * omg(3);
else
    u_p = K_P * e_omgp + K_I * K_P * ie_omgp + K_I * Ixm * e_omgp;
    u_q = K_P * e_omgq + K_I * K_P * ie_omgq + K_I * Iym * e_omgq;
    u_r = K_Pr * e_omgr + K_Ir * K_Pr * ie_omgr + K_Ir * Izm * e_omgr;
end

u_z = -m* 9.781;

uctl = [u_p; u_q; u_r; u_z];

%% ESC control
omega2_u = [ -Lm*Kt1/sqrt(2)  Lm*Kt1/sqrt(2)  Lm*Kt2/sqrt(2)  -Lm*Kt2/sqrt(2);
              Lm*Kt1/sqrt(2) -Lm*Kt1/sqrt(2)  Lm*Kt2/sqrt(2)  -Lm*Kt2/sqrt(2);
                   Kq              Kq              -Kq              -Kq      ;
                  -Kt1            -Kt1             -Kt2             -Kt2      ];      

% A\b is more efficient way to inv(A)*b
omega2 = omega2_u \ [u_p; u_q; u_r; u_z];   

% lower limit to 0
omega2(omega2<0) = 0;

omega  = sqrt(omega2);

%% calculate the PWM signal for the motors based on first order approximation
uin(1) = ( omega(1) - Kw1(2) ) / Kw1(1); 
uin(2) = ( omega(2) - Kw1(2) ) / Kw1(1);
uin(3) = ( omega(3) - Kw2(2) ) / Kw2(1);
uin(4) = ( omega(4) - Kw2(2) ) / Kw2(1);

%% check the range
uin_max = 1;
uin_min = 0;

% limit input signals
uin(uin>uin_max) = uin_max;
uin(uin<uin_min) = uin_min;

%% state update (for simulation only, actual drone will not run this)
% approximate states using 4th order Runge-Kutta (RK4)
for k = 1:Ts/ts
%     if k == 1
%         [k1, a_b] = QrotorDynMod( xn(1:12), uin, wnd );
%         %k1 = QrotorDynMod( xn, uin , [0;0;0] );
%     else
%         k1 = QrotorDynMod( xn(1:12), uin , wnd );
%     end
    k1 = QrotorDynMod( xn(1:12), uin , wnd );
    k2 = QrotorDynMod( xn(1:12)+ts/2*k1, uin , wnd );
    k3 = QrotorDynMod( xn(1:12)+ts/2*k2, uin , wnd );
    k4 = QrotorDynMod( xn(1:12)+ts  *k3, uin , wnd );
    xn(1:12) = xn(1:12) + ts/6 * ( k1 + 2*k2 + 2*k3 + k4 );
    
    % calculate omg_r, since psi is modelled using a first order system 
    xn(12) =  ( k1(9) + 2*k2(9) + 2*k3(9) + k4(9) ) / 6;
end

%% calculate a_b
% xn(7) = wrapToPi(xn(7));
% xn(8) = wrapToPi(xn(8));
% xn(9) = wrapToPi(xn(9));
% phi = xn(7);
% tht = xn(8);
% psi = xn(9); 
% 
% R(1,1) = cos(psi)*cos(tht);
% R(1,2) = cos(psi)*sin(tht)*sin(phi) - sin(psi)*cos(phi);
% R(1,3) = cos(psi)*sin(tht)*cos(phi) + sin(psi)*sin(phi);
% R(2,1) = sin(psi)*cos(tht);
% R(2,2) = sin(psi)*sin(tht)*sin(phi) + cos(psi)*cos(phi);
% R(2,3) = sin(psi)*sin(tht)*cos(phi) - cos(psi)*sin(phi);
% R(3,1) = -sin(tht);
% R(3,2) = cos(tht)*sin(phi);
% R(3,3) = cos(tht)*cos(phi);
% 
% a_b = R' * ( R * xn(4:6) - v_g_true ) / Ts; % differentiate ground velocity

end % end of qrotor_inner_sim function


function varargout = QrotorDynMod( x, u, wnd )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%        Linear Dynamic Model of Quadrotor UAV                           %%%%
%%%%        X-Config Quadrotor with attitude stability augmentation         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [dx,a_b] = QrotorDynMod( x, u, wnd )
% x   = [pos; vel; att; omg] ------------ current system state
% pos = [x; y; z] ----------------------- current ground-axis position, m
% vel = [u; v; w] ----------------------- current body-axis velocity, m/s%
% att = [phi; theta; psi] --------------- current attitude angles, rad
% datt= [dphi; dtht, dpsi] -------------- current derivatives of attitude angles, rad/s
% u   = [u1; u2; u3; u4] ---------------- PWM input in [0 1] without units
% wnd = [wnd_p; wnd_q; wnd_r] ----------- body-axis disturbance at rotating axis, rad/s^2
% dx  = dot of state x

global iTs;
ts = iTs;     % step size for model simulation, sec

persistent g m Lm b;

%% Model parameters

%if isempty(g)
g   = 9.781;        % gravitational acceleration, m/s^2
m   = 4.027;        % total mass of helicopter, kg

Ixx = 0.0707;       % rolling moment of inertia, kg*m^2  (real measured value)
Iyy = 0.0705;       % pitching moment of inertia, kg*m^2 (scaled based on CIFER)
Izz = 0.2322;       % yawing moment of inertia, kg*m^2   (real measured value)

Lm  = 0.303;        % the length of the arm of the platform

b   = -0.544;       % body drag-force coefficient for x and y direction, kg/s (wind tunnel data)

Jr  = 0;            % rotor inertia
OMG_r = 0;          % rotation speed of the motor


%% Parse Function Inputs
u1 = u(1);          % input to the esc
u2 = u(2);
u3 = u(3);
u4 = u(4);

vel_u = x(4);       % m/s
vel_v = x(5);       % m/s
vel_w = x(6);       % m/s

phi   = x(7);       % rad
tht   = x(8);       % rad
psi   = x(9);       % rad
omg_p = x(10);      % rad/s
omg_q = x(11);      % rad/s
omg_r = x(12);      % rad/s

% disturbances
wnd_p = wnd(1);
wnd_q = wnd(2);
wnd_r = wnd(3);

R = eye(3,3); % rotation matrix transforming body-axis vector to ned-axis vector
R(1,1) = cos(psi)*cos(tht);
R(1,2) = cos(psi)*sin(tht)*sin(phi) - sin(psi)*cos(phi);
R(1,3) = cos(psi)*sin(tht)*cos(phi) + sin(psi)*sin(phi);
R(2,1) = sin(psi)*cos(tht);
R(2,2) = sin(psi)*sin(tht)*sin(phi) + cos(psi)*cos(phi);
R(2,3) = sin(psi)*sin(tht)*cos(phi) - cos(psi)*sin(phi);
R(3,1) = -sin(tht);
R(3,2) = cos(tht)*sin(phi);
R(3,3) = cos(tht)*cos(phi);

%% the thrust of the rotors
%  
%  configuration of the uav   
% 
%   3(cw)   1(ccw)
%       \  / 
%        \/
%        /\   
%       /  \
%  2(ccw)   4(cw)

pwm2omg_ccw = [ 449.5054, -923.2888,  1058.7,    8.1576];
omg2Tz_ccw  = [ 0.00012551,   -0.0069,    0.1698];
omg2Rz_ccw  = [ 0.0000037357,  -0.00030354,  0.0095];

pwm2omg_cw  = [ 498.5041, -993.2777,  1075.2,    6.9775];
omg2Tz_cw   = [ 0.00013137,   -0.0077,    0.2126];
omg2Rz_cw   = [ 0.0000040001,  -0.00036572,  0.0112];

Omg1 = pwm2omg_ccw * [ u1^3; u1^2; u1; 1]; 
Omg2 = pwm2omg_ccw * [ u2^3; u2^2; u2; 1]; 
Omg3 = pwm2omg_cw  * [ u3^3; u3^2; u3; 1]; 
Omg4 = pwm2omg_cw  * [ u4^3; u4^2; u4; 1]; 

T_rt1 = omg2Tz_ccw * [ Omg1^2; Omg1; 1]; 
T_rt2 = omg2Tz_ccw * [ Omg2^2; Omg2; 1]; 
T_rt3 = omg2Tz_cw  * [ Omg3^2; Omg3; 1]; 
T_rt4 = omg2Tz_cw  * [ Omg4^2; Omg4; 1]; 

Q_rt1 = omg2Rz_ccw * [ Omg1^2; Omg1; 1];  
Q_rt2 = omg2Rz_ccw * [ Omg2^2; Omg2; 1];  
Q_rt3 = omg2Rz_cw  * [ Omg3^2; Omg3; 1];  
Q_rt4 = omg2Rz_cw  * [ Omg4^2; Omg4; 1];  


%% rotor force and moments
Xrt = 0;
Yrt = 0; 
Zrt = -(T_rt1 + T_rt2 + T_rt3 + T_rt4);

Lrt = Lm / sqrt(2) * ( -T_rt1 + T_rt2 + T_rt3 - T_rt4 );
Mrt = Lm / sqrt(2) * (  T_rt1 - T_rt2 + T_rt3 - T_rt4 );
Nrt = Q_rt1 + Q_rt2 - Q_rt3 - Q_rt4;


%% resulted forces and moments along body axes (actuator action)
% x and y linear axes include wind drag force
F_u = Xrt + b*vel_u;
F_v = Yrt + b*vel_v;
F_w = Zrt;
M_p = Lrt + wnd_p;
M_q = Mrt + wnd_q;
M_r = Nrt + wnd_r;


%% Rigid body dynamics
% body gyro effect + propeller gyro effect + roll/pitch/yaw actuator action
domg_p = ( ( Iyy - Izz ) * omg_q * omg_r + Jr * omg_q * OMG_r + M_p ) / Ixx; 
domg_q = ( ( Ixx - Izz ) * omg_p * omg_r + Jr * omg_p * OMG_r + M_q ) / Iyy;
domg_r = ( ( Ixx - Iyy ) * omg_p * omg_q                      + M_r ) / Izz;

dphi = [1  sin(phi)*tan(tht)  cos(phi)*tan(tht)] * [omg_p; omg_q; omg_r];
dtht = [0       cos(phi)          -sin(phi)    ] * [omg_p; omg_q; omg_r];
dpsi = [0  sin(phi)*sec(tht)  cos(phi)*sec(tht)] * [omg_p; omg_q; omg_r];

agx = R(1,:) * [ F_u; F_v; F_w; ] / m;
agy = R(2,:) * [ F_u; F_v; F_w; ] / m;
agz = R(3,:) * [ F_u; F_v; F_w; ] / m + g;

a_u = F_u / m - g * sin(tht);
a_v = F_v / m + g * sin(phi)*cos(tht);
a_w = F_w / m + g * cos(phi)*cos(tht);

dvel_u = a_u - omg_q*vel_w + omg_r*vel_v;
dvel_v = a_v - omg_r*vel_u + omg_p*vel_w;
dvel_w = a_w - omg_p*vel_v + omg_q*vel_u;

v_g  = R * [vel_u; vel_v; vel_w];
dpos = v_g;

  
%% Function Outputs
dx = [dpos;             dvel_u; dvel_v; dvel_w;
      dphi; dtht; dpsi; domg_p; domg_q; domg_r;];
if nargout == 1
    varargout{1} = dx;
elseif nargout == 2
    varargout{1} = dx;
    varargout{2} = [F_u; F_v; F_w]/m;
end

end % end of qrotordynmdl function


