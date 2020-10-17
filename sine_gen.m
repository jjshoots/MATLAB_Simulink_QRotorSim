%Generates sine signal

function a_ref = sine_gen( cycles, amp, n, Ts )

ref_angle = cycles*2*pi*(0:(n-1))/n;
a_ref  = amp*sin( ref_angle )';

end