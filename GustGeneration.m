%%========================================================================%
% Coded by Phang Swee King                                                %
% Maintained by Phang Swee King (swee_king@hotmail.com)                   %
% Copyright 2017 research UAV team, National University of Singapore      %
% Updated: 6 June 2017                                                    %
%%========================================================================%
function vel_out = GustGeneration(Ts,time_period,v_max)
% The Discrete Wind Gust Model block implements a wind gust of the standard
% "1-cosine" shape. This block implements the mathematical representation
% in the Military Specification MIL-F-8785C. The gust is applied to
% each axis individually, or to all three axes at once. You specify the
% gust amplitude (the increase in wind speed generated by the gust), the
% gust length (length, in meters, over which the gust builds up) and the
% gust start time.

t = time_period;
g_l = floor(t/3); % for simplicity, gust_length assumed to be 1/3 of the total wind period
v_m = v_max;

length_t = t/Ts;
vel = zeros(1,length_t);

for i = 1:g_l/Ts
    vel(i) = (v_m/2)*(1-cos(pi*i*Ts/g_l));
end
for i = g_l/Ts:length_t-g_l/Ts
    vel(i) = v_m;
end
for i = length_t-g_l/Ts:length_t
    vel(i) = (v_m/2)*(1-cos(pi*(length_t-i)*Ts/g_l));
end

vel_out = vel;

end