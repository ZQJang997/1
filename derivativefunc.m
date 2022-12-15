%% derivative function that returns nonlinear xdot when given x_ and u_ 

function xdot_nonlinear = derivativefunc(x, u, t)
    % CONTROL GAINS
    Kpsibyx = 0.011;   % deg/(ft ⋅sec)
    Kpsi = 0.0667;     % deg/(deg⋅sec)
    Kv = 1;            % kt/kt/sec

    % Now, assembling the derivatives of each of the states.
    % NOTE: shown for readability, not computational efficiency

    % derivative of x, with V converted from knots to ft per sec
    xderivative = x(3) * sind(x(4)) * 6076/3600;

    % derivative of y, with V converted from knots to ft per sec
    yderivative = x(3) * cosd(x(4)) * 6076/3600;

    % derivative of V - a simple proportional control law on error in V
    Vderivative = Kv * (u(2) - x(3));

    % derivative of psi - a simple proportional control law on error in x & psi
    psiderivative = Kpsibyx * (u(1) - x(1)) + Kpsi * (u(3) - x(4));

    % put limits on the YAW rate
    psiderivative = min(max(psiderivative, -3), 3);

    xdot_nonlinear = [xderivative; yderivative; Vderivative; psiderivative];
end