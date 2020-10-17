%Generates step signal

function a_ref = step_gen( step_time, step_amp, n, Ts )

a_ref = zeros(n, 1);
a_ref(step_time/Ts:n) = step_amp;

end


