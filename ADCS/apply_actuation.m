function T = apply_actuation(x,u,actuators,cpu_speed)
wG_true = x(14:16);
wG_commanded = u;
wGd_commanded = (wG_commanded-wG_true)/cpu_speed;
T = actuators.JWHEEL*wGd_commanded; % commanded torque
if (T > actuators.RWA_MAXTORQUE)
    T = actuators.RWA_MAXTORQUE;
end
end

