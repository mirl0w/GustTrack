function mppt_wind()
% MPPT comparison for a small wind turbine (simulation).
% Controllers: Fixed resistive load (baseline), Optimal-TSR torque control, Perturb & Observe (P&O)
% Replace cpCurve() with YOUR blade's Cp(lambda) from BEM/QBlade for the final version.
% Save as mppt_wind.m and run:  >> mppt_wind

    % ---------------- Parameters ----------------
    P.rho = 1.225;          % air density [kg/m^3]
    P.R   = 0.6;            % rotor radius [m]
    P.A   = pi * P.R^2;     % swept area [m^2]
    P.J   = 0.15;           % rotor + generator inertia [kg m^2]
    P.eta = 0.85;           % generator + converter efficiency
    P.dt  = 0.01;           % time step [s]
    P.Tend = 600;           % duration [s]
    P.vmean = 7.0;          % mean wind speed [m/s]

    lams = linspace(0.5, 14, 2000);
    [P.CPMAX, idx] = max(cpCurve(lams));
    P.LAMOPT = lams(idx);
    P.Kopt = 0.5 * P.rho * P.A * P.R^3 * P.CPMAX / P.LAMOPT^3;   % Te = Kopt*w^2

    % ---------------- Main run (seed 1) ----------------
    [t, v] = windProfile(P, 1);
    names  = {'fixed', 'tsr', 'po'};
    labels = {'Fixed load', 'Optimal TSR', 'P&O'};
    res = struct();
    for i = 1:3
        res.(names{i}) = runSim(names{i}, v, P);
    end

    pwind  = 0.5 * P.rho * P.A * v.^3;
    Eideal = sum(P.eta * P.CPMAX * pwind) * P.dt / 3600;          % [Wh]
    fprintf('LAMBDA_OPT=%.2f  CP_MAX=%.3f  Kopt=%.5f\n', P.LAMOPT, P.CPMAX, P.Kopt);
    fprintf('Ideal (always at Cpmax, eta applied): %.2f Wh\n', Eideal);
    E = struct();
    for i = 1:3
        r = res.(names{i});
        E.(names{i}) = sum(r.P) * P.dt / 3600;
    end
    for i = 1:3
        n = names{i}; r = res.(n);
        fprintf('%-6s energy = %7.2f Wh | %5.1f%% of ideal | %+5.1f%% vs fixed | mean lambda=%.2f\n', ...
            n, E.(n), 100*E.(n)/Eideal, 100*(E.(n)/E.fixed - 1), mean(r.lam));
    end

    % ---------------- Robustness: 5 wind seeds ----------------
    fprintf('\nMean over 5 random wind seeds (%% of ideal energy):\n');
    acc = zeros(5, 3);
    for sd = 1:5
        [~, vv] = windProfile(P, sd);
        ideal = sum(P.eta * P.CPMAX * 0.5 * P.rho * P.A * vv.^3);
        for i = 1:3
            r = runSim(names{i}, vv, P);
            acc(sd, i) = 100 * sum(r.P) / ideal;
        end
    end
    for i = 1:3
        fprintf('  %-6s %5.1f +/- %.1f\n', names{i}, mean(acc(:, i)), std(acc(:, i)));
    end

    % ---------------- Plots ----------------
    figure('Name', 'MPPT results', 'Position', [100 100 1200 750]);
    subplot(2,2,1);
    plot(lams, cpCurve(lams), 'LineWidth', 1.5); hold on;
    xline(P.LAMOPT, '--r');
    title(sprintf('Cp vs TSR (max %.2f @ %.1f)', P.CPMAX, P.LAMOPT));
    xlabel('TSR \lambda'); ylabel('C_p'); grid on;

    subplot(2,2,2);
    plot(t, v); title('Wind speed'); xlabel('s'); ylabel('m/s'); grid on;

    subplot(2,2,3); hold on;
    for i = 1:3
        plot(t, smooth200(res.(names{i}).P), 'DisplayName', labels{i});
    end
    title('Electrical power (smoothed)'); xlabel('s'); ylabel('W'); legend; grid on;

    subplot(2,2,4); hold on;
    for i = 1:3
        plot(t, smooth200(res.(names{i}).cp), 'DisplayName', labels{i});
    end
    title('C_p achieved (smoothed)'); xlabel('s'); ylabel('C_p'); legend; grid on;

    % ---------------- Export log for the telemetry/dashboard part ----------------
    log = [t(:), v(:), res.tsr.w(:), res.tsr.P(:)];
    writematrix(log, 'sim_log.csv');   % columns: time, wind, rpm_rad_s, power_W
    fprintf('\nSaved sim_log.csv\n');

    % ---------------- Mean wind speed study: 5, 7, 9 m/s ----------------
    fprintf('\nMean wind speed study (%% of ideal energy, mean over 5 seeds):\n');
    for vm = [5 7 9]
        P.vmean = vm;
        acc = zeros(5, 3);
        for sd = 1:5
            [~, vv] = windProfile(P, sd);
            ideal = sum(P.eta * P.CPMAX * 0.5 * P.rho * P.A * vv.^3);
            for i = 1:3
                r = runSim(names{i}, vv, P);
                acc(sd, i) = 100 * sum(r.P) / ideal;
            end
        end
        fprintf('vmean=%d m/s: fixed %.1f | tsr %.1f | po %.1f\n', vm, mean(acc));
    end
end

% =====================================================================
function cp = cpCurve(lam)
% Generic Cp(lambda), beta = 0 (Heier model). Max ~0.48 at lambda ~ 8.1
    lam = max(lam, 1e-3);
    li  = 1 ./ (1 ./ lam - 0.035);
    li(li <= 0) = 1e3;
    cp = 0.5176 * (116 ./ li - 5) .* exp(-21 ./ li) + 0.0068 * lam;
    cp = max(cp, 0);
end

% =====================================================================
function [t, v] = windProfile(P, seed)
    rng(seed);
    n = round(P.Tend / P.dt);
    t = (0:n-1) * P.dt;
    slow = 1.5 * sin(2*pi*t/120) + 1.0 * sin(2*pi*t/37);
    tau = 3.0; sig = 0.9;                       % Ornstein-Uhlenbeck turbulence
    ou = zeros(1, n);
    for i = 2:n
        ou(i) = ou(i-1) - ou(i-1) * P.dt / tau + sig * sqrt(2 * P.dt / tau) * randn();
    end
    v = min(max(P.vmean + slow + ou, 2.5), 13.0);
end

% =====================================================================
function out = runSim(controller, v, P)
    n = numel(v);
    w = P.LAMOPT * P.vmean / P.R;               % start near optimum for the mean wind [rad/s]
    out.P = zeros(1, n); out.w = zeros(1, n);
    out.lam = zeros(1, n); out.cp = zeros(1, n); out.k = zeros(1, n);

    cFixed = P.Kopt * (P.LAMOPT * P.vmean / P.R);   % fixed resistive load: Te = c*w (tuned to the mean wind)

    % P&O state: perturb torque gain k in Te = k*w^2
    k = 0.7 * P.Kopt;  dk = 0.02 * P.Kopt;  TsPO = 8.0;
    pAcc = 0; pPrev = 0; dirn = 1; cnt = 0;
    nHalf = floor(TsPO / P.dt / 2);

    for i = 1:n
        lam = w * P.R / v(i);
        cp  = cpCurve(lam);
        Tm  = 0.5 * P.rho * P.A * P.R * v(i)^2 * cp / max(lam, 1e-3);   % aerodynamic torque

        switch controller
            case 'fixed'
                Te = cFixed * w;
            case 'tsr'
                Te = P.Kopt * w^2;
            case 'po'
                if cnt >= nHalf
                    pAvg = pAcc / cnt;
                    if pAvg < pPrev
                        dirn = -dirn;
                    end
                    k = min(max(k + dirn * dk, 0.3 * P.Kopt), 2.0 * P.Kopt);
                    pPrev = pAvg; pAcc = 0; cnt = 0;
                end
                Te = k * w^2;
        end
        Te = max(Te, 0);

        w  = max(w + (Tm - Te) / P.J * P.dt, 0.5);     % rotor dynamics (Euler)
        Pe = P.eta * Te * w;

        if strcmp(controller, 'po')
            pAcc = pAcc + Pe; cnt = cnt + 1;
        end
        out.P(i) = Pe; out.w(i) = w; out.lam(i) = lam; out.cp(i) = cp; out.k(i) = k;
    end
end

% =====================================================================
function y = smooth200(x)
    y = conv(x, ones(1, 200) / 200, 'same');
end