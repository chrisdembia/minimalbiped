function minimalbiped
    % MINIMALBIPED Implementation of the minimal biped model in Srinivasan,
    % 2006), by Chris Dembia.

    global counter;
    counter = 0;

    p.lmax = 1; % m

    p.V = 0.5;
    p.D = 0.5;

    p.m = 1; % kg
    p.g = 9.8; % m/s^2

    p.N = 10;

    % Constraints.
    Fbarmax = 1; % (-)
    Fbarmin = 0;
    tauSmin = p.D / p.V; % TODO 1e-8;
    tauSmax = p.D / p.V;

    % Guesses.
    l0 = p.lmax;
    X00 = -p.D/2;
    Xd00 = 0.5;
    Y00 = sqrt(-(0.5*p.D)^2 + l0^2);
    Yd00 = -0.15;
    Fbar = [0; 0 * ones(p.N - 2, 1); 1];
    tauS = tauSmax;
    arg0 = [X00; Xd00; Y00; Yd00; Fbar; tauS];
    
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [-inf;
          -inf;
          -inf;
          -inf;
          Fbarmin * ones(p.N, 1),
          tauSmin];
    ub = [inf;
          inf;
          inf;
          inf;
          [Fbarmin; zeros(p.N - 2, 1); Fbarmax]; % TODO Fbarmax * ones(p.N, 1);
          tauSmax];

    argLast = [];
    myf = [];
    myc = [];
    myceq = [];

    fun = @(arg) objfun(arg, p);
    cfun = @(arg) constr(arg, p);

    hf = figure;
    set(hf, 'Position', [500, 400, 1000, 600]);
    options = optimset('Algorithm', 'sqp');
    [optArg Cmin] = fmincon(fun, arg0, A, b, Aeq, beq, lb, ub, cfun, options);

    [t, future] = simulation(optArg, p);
    x = future(:, 1);
    xd = future(:, 2);
    y = future(:, 3);
    yd = future(:, 4);

    figure;
    for i = 1:length(t)
        plot(x(i), y(i), 'o');
        axis equal;
        hold on;
    end

    X0 = optArg(1);
    Xd0 = optArg(2);
    Y0 = optArg(3);
    Yd0 = optArg(4);
    Fbar = optArg(4+1:4+p.N);
    tauS = optArg(end);

    figure;
    plot(Fbar);

    function obj = objfun(arg, p)
        if ~isequal(arg, argLast)
            [myf, myc, myceq] = computeall(arg, p);
            argLast = arg;
        end
        obj = myf;
    end

    function [c, ceq] = constr(arg, p)
        if ~isequal(arg, argLast)
            [myf, myc, myceq] = computeall(arg, p);
        end
        c = myc;
        ceq = myceq;
    end
end

function [d, v] =  dimizeDV(p)

    d = p.D * p.lmax;
    v = p.V * sqrt(p.g * p.lmax);

end

function [C, c, ceq] = computeall(arguments, p)

    global counter;
    counter = counter + 1;

    % Simulate.
    [t, future] = simulation(arguments, p);

    % Unpack results.
    x = future(:, 1);
    xd = future(:, 2);
    y = future(:, 3);
    yd = future(:, 4);

    W_runningSum = future(:, 5);

    dActual = x(end) - x(1);
    vActual = trapz(t, x) / (t(end) - t(1));

    % Compute objective function.
    W = W_runningSum(end);
    C = W / p.m / p.g / dActual;

    % Compute constraint values.
    % Leg length throughout the simulation.
    l = sqrt(x.^2 + y.^2);

    [dDesired, vDesired] = dimizeDV(p);

    % Enforced at the grid points.
    [d, v] = dimizeDV(p);
    tStep = d / v;
    %c = [linterp(t/tStep, l, p.tauGrid) - p.lmax];
    c = [l - p.lmax];

    % TODO do not need to enforce the velocity constraint.
    ceq = [dActual - dDesired;
           y(end) - y(1);
           xd(end) - xd(1);
           yd(end) - yd(1)];

    % See intermediate results.
    clf;
    subplot(2, 4, 1);
    for i = 1:length(t)
        plot(x(i), y(i), 'o', 'Color', [i/length(t), 0 0]);
        axis equal;
        axis([min(x) max(x) 0 max(y)]);
        hold on;
    end

    subplot(2, 4, 2);
    Fbar = arguments(4+1:4+p.N);
    plot(Fbar);
    title('fbar');

    subplot(2, 4, 3);
    plot(t, l);
    title('l(tau)');

    subplot(2, 4, 4);
    plot(t, W_runningSum);
    title('W running sum');

    subplot(2, 4, 5);
    axis off;
    text(0, 0.5, ...
        {sprintf('iteration: %i', counter), ...
        sprintf('cost: %f', C), ...
        sprintf('xi -> xf: %f -> %f', x(1), x(end)), ...
        sprintf('(dDes, dAct: %f , %f)', dDesired, dActual), ...
        sprintf('yi -> yf: %f -> %f', y(1), y(end)), ...
        sprintf('xdi -> xdf: %f -> %f', xd(1), xd(end)), ...
        sprintf('ydi -> ydf: %f -> %f', yd(1), yd(end)), ...
        sprintf('tauS: %f', arguments(end))});
    pause(.5);

end

function [t, future] = simulation(arguments, p)
    X0 = arguments(1);
    Xd0 = arguments(2);
    Y0 = arguments(3);
    Yd0 = arguments(4);
    Fbar = arguments(4+1:4+p.N);
    tauS = arguments(end);

    nonDimSpeed = sqrt(p.g * p.lmax);

    x0 = X0 * p.lmax;
    xd0 = Xd0 * nonDimSpeed;
    y0 = Y0 * p.lmax;
    yd0 = Yd0 * nonDimSpeed;
    p.Fgrid = Fbar * p.m * p.g;
    p.tS = tauS * sqrt(p.lmax / p.g);

    W0 = 0;

    initConditions = [x0; xd0; y0; yd0; W0];

    [d, v] = dimizeDV(p);

    p.tGrid = linspace(0, p.tS, p.N);

    tStep = d / v;
    tSpan = [0 tStep];

    options = odeset('Events', @(t, z) startSwing(t, z, p));
    [t, future] = ode45(@(t, z) eqnsOfMotion(t, z, p), ...
        tSpan, initConditions, options);
end

function [value, isterminal, direction] = startSwing(t, z, p)

    global isStance;

    isStance = t <= p.tS;

    value =  t - p.tS;
    isterminal = false;
    direction = 1;

end

function derivatives = eqnsOfMotion(t, z, p)

    global isStance;

    x = z(1);
    xd = z(2);
    y = z(3);
    yd = z(4);

    xc = 0;

    l = sqrt((x - xc)^2 + y^2);
    ldot = (2*(x - xc)*xd + 2*y*yd) / (2 * sqrt((x - xc)^2 + y^2));

    isStance = t <= p.tS;

    if isStance
        FatT = linterp(p.tGrid, p.Fgrid, t);
        xdd = 1 / p.m * FatT * (x - xc) / l;
        ydd = -p.g + 1 / p.m * FatT * y / l;
        Wdot = ((FatT * ldot) > 0) * FatT * ldot;
    else
        xdd = 0;
        ydd = -p.g;
        Wdot = 0;
    end

    derivatives = [xd; xdd; yd; ydd; Wdot];

end
