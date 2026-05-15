function [P, A, q, l, u, scaling] = modified_ruiz_equilibrate(P, A, q, l, u, scalingIterations, options)
    %MODIFIED_RUIZ_EQUILIBRATE MATLAB translation of OSQP's scale_data routine.
    %   [P, A, q, l, u, scaling] = modified_ruiz_equilibrate(P, A, q, l, u, scaling_iter)
    %   applies Ruiz equilibration and cost normalization similarly to
    %   src/scaling.c::scale_data in OSQP.
    %
    %   Inputs:
    %     P, A         Problem matrices (dense or sparse)
    %     q, l, u      Problem vectors
    %     scaling_iter Number of scaling iterations (same meaning as OSQP setting)
    %
    %   Outputs:
    %     P, A, q, l, u  Scaled problem data
    %     scaling        Struct with fields c, cinv, D, Dinv, E, Einv
    %
    %   Notes:
    %     - OSQP constants are mirrored: OSQP_MIN_SCALING=1e-4, OSQP_MAX_SCALING=1e4.
    %     - Infinity norms are column/row max abs values, matching OSQP behavior.

    %% Input validation
    arguments
        P (:,:) double;
        A (:,:) double;
        q (:,1) double;
        l (:,1) double;
        u (:,1) double;
        scalingIterations (1, 1) {mustBePositive, mustBeInteger} = 10;
        options.minScaling (1, 1) {mustBeReal, mustBeFinite, mustBePositive} = 1e-4;
        options.maxScaling (1, 1) {mustBeReal, mustBeFinite, mustBePositive} = 1e4;
    end

    n = size(P, 1);
    m = size(A, 1);
    minScaling = options.minScaling;
    maxScaling = options.maxScaling;

    assert(n == size(P, 2), 'modified_ruiz_equilibrate:DimensionMismatch', 'P must be square, but got size %d-by-%d.', size(P, 1), size(P, 2));
    assert(size(A, 2) == n, 'modified_ruiz_equilibrate:DimensionMismatch', 'A must have %d columns to match the size of P (%d-by-%d).', n, n, n);
    assert(minScaling <= maxScaling, 'modified_ruiz_equilibrate:InvalidScalingBounds', 'options.minScaling must be less than or equal to options.maxScaling.');
    assert(numel(q) == n, 'modified_ruiz_equilibrate:DimensionMismatch', 'q must have %d elements to match the number of columns in P and A.', n);
    assert(numel(l) == m && numel(u) == m, 'modified_ruiz_equilibrate:DimensionMismatch', 'l and u must each have %d elements to match the number of rows in A.', m);

    %% 
    scaling.c = 1.0;
    scaling.D = ones(n, 1);
    scaling.E = ones(m, 1);

    for k = 1:scalingIterations
        % Compute column infinity norms of KKT blocks without forming KKT.
        col_norm_P = vecnorm(P, inf); % 1 x n
        col_norm_A = vecnorm(A, inf); % 1 x n
        D_temp = max(col_norm_P, col_norm_A); % 1 x n

        % Row infinity norms for A map to the second KKT block.
        E_temp = vecnorm(A, inf, 2); % m x 1

        % Clamp scaling vectors.
        D_temp = clampScaling(D_temp, minScaling, maxScaling);
        E_temp = clampScaling(E_temp, minScaling, maxScaling);

        % Inverse square-root scaling.
        D_temp = 1 ./ sqrt(D_temp);
        E_temp = 1 ./ sqrt(E_temp);

        % Equilibrate matrices and vector.
        %    n x 1  n x n   1 x n
        P = D_temp' .* P .* D_temp; % n x n
        %   m x 1   m x n  1 x n
        A = E_temp .* A .* D_temp; % m x n
        %  n x 1 n x 1
        q = q .* D_temp'; % n x 1

        % Accumulate total scaling.
        %               n x 1              n x 1
        scaling.D = scaling.D .* D_temp';
        %               m x 1              m x 1 
        scaling.E = scaling.E .* E_temp;

        % Cost normalization.
        col_norm_P = vecnorm(P, inf);
        c_temp = mean(col_norm_P);
        inf_norm_q = norm(q, inf);
        inf_norm_q = clampScaling(inf_norm_q, minScaling, maxScaling);
        c_temp = max(c_temp, inf_norm_q);
        c_temp = clampScaling(c_temp, minScaling, maxScaling);
        c_temp = 1.0 / c_temp;
        P = c_temp * P;
        q = c_temp * q;
        scaling.c = scaling.c * c_temp;
    end

    scaling.cinv = 1.0 / scaling.c;
    scaling.Dinv = 1 ./ scaling.D;
    scaling.Einv = 1 ./ scaling.E;

    l = l .* scaling.E;
    u = u .* scaling.E;

end

function v = clampScaling(v, minValue, maxValue)
    v = min(v, maxValue);
    v(v < minValue) = 1.0;
end
