function [cost_fps, cost_std, cost_inc] = CostStep(N_s, N_ip, N_g)
    % 1. FPS Steps (k = 1 to N_ip)
    k_fps = 1:N_ip;
    cost_fps = 3 .* (N_s - k_fps);

    % 2. Greedy Steps (k = N_ip + 1 to N_ip + N_g)
    k_greedy = (N_ip + 1):(N_ip + N_g);

    % Standard Greedy step cost
    cost_std = (1/3) .* k_greedy.^3 + 2.5 .* k_greedy.^2 + 2 .* k_greedy .* (N_s - k_greedy);

    % Incremental Greedy step cost
    cost_inc = 3 .* k_greedy.^2 + 2 .* (N_s - k_greedy);

    % 3. ADD INITIAL CHOLESKY FACTORIZATION TO THE FIRST GREEDY STEP
    % This represents the setup cost right as the greedy loop starts
    cost_chol_init = (1/3) * N_ip^3;
    cost_std(1) = cost_std(1) + cost_chol_init;
    cost_inc(1) = cost_inc(1) + cost_chol_init;
end