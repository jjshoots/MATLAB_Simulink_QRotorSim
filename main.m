%%========================================================================%
% Coded by Phang Swee King                                                %
% Maintained by Phang Swee King (swee_king@hotmail.com)                   %
% Copyright 2018 Taylor's Robotic Center, Taylor's University             %
% Updated: 9 Aug 2018             Version 1.0.3                           %
%%========================================================================%

% close all;
clc;

% User-input Variable
% Global to be used in functions/script
global K_P K_I K_Ir K_Pr iTs SSI_or_PI;
global Ixx Iyy Izz Ixm Iym Izm;
K_P  = 3;
K_I  = 3;
K_Pr = 3;
K_Ir = 50;
SSI_or_PI = 1; % 1 for SSI, 0 for PI

% actual J
Ixx = 0.0707;       % rolling moment of inertia,  kg*m^2 (real measured value)
Iyy = 0.0705;       % pitching moment of inertia, kg*m^2 (real measured value)
Izz = 0.2322;       % yawing moment of inertia,   kg*m^2 (real measured value)

% measured J
Ixm = Ixx*2;
Iym = Iyy*2;
Izm = Izz*2;

Td = 150;   % total simulation time, s
Ts = 0.02;  % update rate for simulation, s
iTs = 0.005; % Runge Kutta update rate, s

wind_x_start = 20;    % x axis wind start time (s)
wind_x_end   = 140;    % x axis wind end time (s)
wind_x_vm    = 0.5;     % x axis wind max acceleration (rad/s^2)
wind_y_start = 0;    % y axis wind start time (s)
wind_y_end   = Td;   % y axis wind end time (s)
wind_y_vm    = 0;     % y axis wind max acceleration (rad/s^2)
wind_z_start = 0;    % z axis wind start time (s)
wind_z_end   = Td;    % z axis wind end time (s)
wind_z_vm    = 0;     % z axis wind max acceleration (rad/s^2)

% Reference generation (angular reference)
t   = 0:Ts:Td;
n   = length(t);

% Input signal generation
ref_angle = step_gen(10, 1, n, Ts);

phi_r = ref_angle;
%phi_r = zeros(n, 1);
tht_r = zeros(n, 1);
psi_r = zeros(n, 1);

% State
xn  = zeros(12, 1);   % 12 states in total, pos, vel, att, omg
xc  = xn;            % current states, initialized as xn

% Integral terms
Itg = zeros(3, 1);    % integral states for next input

% Wind disturbance
wnd_x = zeros(n, 1);
wnd_y = zeros(n, 1);
wnd_z = zeros(n, 1);

wnd_x(wind_x_start/Ts+1:wind_x_end/Ts) = GustGeneration(Ts,wind_x_end-wind_x_start,wind_x_vm);
wnd_y(wind_y_start/Ts+1:wind_y_end/Ts) = GustGeneration(Ts,wind_y_end-wind_y_start,wind_y_vm);
wnd_z(wind_z_start/Ts+1:wind_z_end/Ts) = GustGeneration(Ts,wind_z_end-wind_z_start,wind_z_vm);

wnd = [wnd_x(1:Td/Ts+1) wnd_y(1:Td/Ts+1) wnd_z(1:Td/Ts+1)];
% wnd = zeros(n, 3);    % currently no wind

%alpha holders
alpha = zeros(3, 1);
alpha_predix = zeros(3, 1);
delta_alpha = zeros(3, 1);

% Data logging
xn_all   = zeros(n, 12); 
uin_all  = zeros(n, 4); 
uctl_all = zeros(n, 4); 
Itg_all  = zeros(n, 3);

alpha_all = zeros(n, 3);
alpha_predix_all = zeros(n, 3);
delta_alpha_all = zeros(n, 3);

k_learn = 0.001;

%inertia holder
Ixx_all = zeros(n, 3);

display('Simulation Starts...');
for i = 1:n
    if rem(i, 100) == 1
        fprintf('%f / %d seconds total simulated. \n',Ts*(i-1), Td);
    end
    
    % setup reference
    ref   = [phi_r(i), tht_r(i), psi_r(i)];
    
    % simulation starts
    [xn, uin, uctl, Itg] = Qrotor_inner_sim( xn, xc, ref, Itg, Ts, wnd(i,:));
    
    % updated xn is new current state xc
    xc = xn;
    
    % make alpha prediction
    alpha_predix(1) = abs((Iym - Izm) * xn(11) * xn(12) + uctl(1)) / Ixm;
    alpha_predix(2) = abs((Ixm - Izm) * xn(10) * xn(12) + uctl(2)) / Iym;
    alpha_predix(3) = abs((Ixm - Iym) * xn(10) * xn(11) + uctl(3)) / Izm;
    
    if i ~= 1
        alpha = [abs(xn(10:12)' - xn_all(i-1, 10:12)) ./ Ts]';
        delta_alpha = alpha - alpha_predix;
    end
    
    % saving all in _all variables
    xn_all(i,:)   = xn';
    Itg_all(i,:)  = Itg';
    uin_all(i,:)  = uin';
    uctl_all(i,:) = uctl';
    
    alpha_predix_all(i, :) = alpha_predix';
    alpha_all(i, :) = alpha';
    delta_alpha_all(i, :) = delta_alpha';
    
    Ixm = Ixm - k_learn * delta_alpha(1);
    Iym = Iym - k_learn * delta_alpha(2);
    Izm = Izm - k_learn * delta_alpha(3);
    Ixx_all(i, 1) = Ixm;
    Iyy_all(i, 2) = Iym;
    Izz_all(i, 3) = Izm;
end

display('Simulation Ends!!');


% Parse function outputs
x_m    = xn_all(:,1);
y_m    = xn_all(:,2);
z_m    = xn_all(:,3);
u_m    = xn_all(:,4);
v_m    = xn_all(:,5);
w_m    = xn_all(:,6);
phi_m  = xn_all(:,7);
tht_m  = xn_all(:,8);
psi_m  = xn_all(:,9);
omgp_m = xn_all(:,10);
omgq_m = xn_all(:,11);
omgr_m = xn_all(:,12);

ail_m  = uin_all(:,1);
ele_m  = uin_all(:,2);
thr_m  = uin_all(:,3);
rud_m  = uin_all(:,4);

u_p    = uctl_all(:,1);
u_q    = uctl_all(:,2);
u_r    = uctl_all(:,3);

%% Plot graphs
figure(1)
plot(t, wnd(:,1), 'k:', t, xn_all(:,7), 'g', t, phi_m, 'b', t, phi_r, 'r--', t, Itg_all(:, 1), 'g--');
ylim([-2,2])
legend('Disturbance (in rad/s^2)', 'Angle response PI (in rad)', 'Angle response SIPIC (in rad)', 'Angle reference (in rad)', 'asdsa')
xlabel('Time (s)');
ylabel('Angle (rad) / Disturbance (rad/s^2)');

% plot position
figure(2)
subplot(3,1,1)
plot(t, x_m);
title('x')
subplot(3,1,2)
plot(t, y_m);
title('y')
subplot(3,1,3)
plot(t, z_m);
title('z')

% plot velocity
figure(3)
subplot(3,1,1)
plot(t, u_m);
title('u')
subplot(3,1,2)
plot(t, v_m);
title('v')
subplot(3,1,3)
plot(t, w_m);
title('w')

% plot error
figure(4)
subplot(3,1,1)
plot(t, phi_r-phi_m);
title('Phi Error')
ylim([-0.05 0.05])
subplot(3,1,2)
plot(t, tht_r-tht_m);
title('Theta Error')
subplot(3,1,3)
plot(t, psi_r-psi_m);
title('Psi Error')


% plot angle
figure(5)
subplot(3,1,1)
plot(t, phi_m, 'b', t, phi_r, 'r--');
title('Phi')
subplot(3,1,2)
plot(t, tht_m, 'b', t, tht_r, 'r--');
title('Theta')
subplot(3,1,3)
plot(t, psi_m, 'b', t, psi_r, 'r--');
title('Psi')

% plot angular rate
% figure(6)
% subplot(3,1,1)
% plot(t, omgp_m);
% title('Omega p')
% subplot(3,1,2)
% plot(t, omgq_m);
% title('Omega q')
% subplot(3,1,3)
% plot(t, omgr_m);
% title('Omega r')

% figure(7)
% plot(t, alpha_all(:, 1), 'b', t, alpha_predix_all(:, 1), 'r--');

% figure(10)
% plot(t, delta_alpha_all(:, 1), 'b');

figure(11)
plot(t, (Ixx - Ixx_all(:, 1)) ./ Ixx, 'r--', t, (Iyy - Ixx_all(:, 2)) ./ Iyy, 'g--', t, (Izz - Ixx_all(:, 3)) ./ Izz, 'b--');

% figure(12)
% plot(t(1:n-1), diff(xn_all(:, 10))./Ts);