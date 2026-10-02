function telemetry_link()
% Telemetry over an unreliable link: packet loss + delay.
    D = readmatrix('sim_log.csv');
    t = D(:,1);  v = D(:,2);  w = D(:,3);  P = D(:,4);
    n  = numel(t);
    dt = t(2) - t(1);                       % simulation step [s] (0.01)

    Ts = 0.1;                               % telemetry period [s] -> 10 packets per second
    scenarios = [ 0    0;                   % [loss rate, one-way delay in s]
                  0    0.5;
                  0.05 0.5;
                  0.20 0.5;
                  0.20 2.0 ];

    Etrue = sum(P) * dt / 3600;             % true energy [Wh]
    fprintf('True energy: %.2f Wh\n\n', Etrue);
    fprintf('loss   delay   energy error   NRMSE of power\n');

    Pm_all = cell(size(scenarios, 1), 1);
    for s = 1:size(scenarios, 1)
        loss  = scenarios(s, 1);
        delay = scenarios(s, 2);
        Pm = monitorSignal(P, dt, Ts, loss, delay, 1);   % what the monitor sees
        Pm_all{s} = Pm;

        Emon  = sum(Pm) * dt / 3600;                     % energy the monitor would report
        errE  = 100 * (Emon - Etrue) / Etrue;            % energy error [%]
        nrmse = 100 * sqrt(mean((Pm - P).^2)) / mean(P); % instantaneous error [% of mean power]
        fprintf('%4.0f%%  %4.1fs   %+8.2f %%     %6.1f %%\n', 100*loss, delay, errE, nrmse);
    end

    % ---------------- Plot: true vs monitor (20% loss, 2 s delay) ----------------
    win = (t >= 100 & t <= 160);
    figure('Name', 'Telemetry link', 'Position', [100 100 1100 450]);
    plot(t(win), P(win), 'LineWidth', 1.2, 'DisplayName', 'True power'); hold on;
    plot(t(win), Pm_all{5}(win), 'LineWidth', 1.2, 'DisplayName', 'Monitor (20% loss, 2 s delay)');
    plot(t(win), Pm_all{3}(win), '--', 'LineWidth', 1.0, 'DisplayName', 'Monitor (5% loss, 0.5 s delay)');
    xlabel('Time [s]'); ylabel('Power [W]'); legend; grid on;
    title('True turbine power vs what the remote monitor sees');
    saveas(gcf, 'telemetry_results.png');
end

% =====================================================================
function Pm = monitorSignal(P, dt, Ts, loss, delay, seed)
% Sample P every Ts seconds, drop packets randomly, delay the rest,
% and hold the last received value (zero-order hold) between arrivals.
    rng(seed);
    n    = numel(P);
    idx  = 1:round(Ts / dt):n;                          % sampling instants (indices)
    ok   = rand(1, numel(idx)) >= loss;                 % true = packet delivered
    arr  = idx + round(delay / dt);                     % arrival index of each packet
    good = ok & (arr <= n);

    sig = nan(n, 1);
    sig(arr(good)) = P(idx(good));                      % place delivered values at arrival time

    Pm   = zeros(n, 1);
    last = 0;                                           % monitor shows 0 until first packet arrives
    for i = 1:n
        if ~isnan(sig(i))
            last = sig(i);
        end
        Pm(i) = last;
    end
end