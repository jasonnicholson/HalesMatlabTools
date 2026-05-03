classdef newton_realStressTest < matlab.unittest.TestCase
    %newton_realStressTest  Stress-test newton_real with randomly generated polynomials.
    %
    % Generates N_TESTS monic polynomials with i.i.d. N(0,1) non-leading
    % coefficients and degrees drawn uniformly from [MIN_DEG, MAX_DEG].
    % All polynomial evaluations and assertions are computed in
    % TestMethodSetup so that each Test method contains only assertions.
    %
    % The pass criterion is a *scaled* residual:
    %
    %     max|p(r_k)| / max(1, max|r_k|^N) < RESIDUAL_TOL
    %
    % Scaling by max(1, max|r|^N) is necessary because |p(r)| can grow as
    % |r|^N when an evaluated root is even slightly off in modulus, even
    % though that root is correct to full machine precision. The scaled
    % residual stays ~eps for accurate roots regardless of |r|.
    %
    % See also: newton_real, newton_realTest

    % =====================================================================
    %  Configuration constants
    % =====================================================================
    properties (Constant)
        N_TESTS             = 200      % number of random polynomials
        MIN_DEG             = 4        % smallest degree
        MAX_DEG             = 25       % largest degree
        RESIDUAL_TOL        = 1e-6     % max scaled |p(root)| threshold for pass
        RANDOM_SEED         = 42       % fixed seed for reproducibility
        OVERALL_PASS_RATE   = 0.98     % required overall success rate
        LOW_DEG_PASS_RATE   = 1.00     % required rate for degrees  4–10
        MID_DEG_PASS_RATE   = 0.97     % required rate for degrees 11–18
        HIGH_DEG_PASS_RATE  = 0.97     % required rate for degrees 19–25
        MEDIAN_RESIDUAL_MAX = 1e-12    % required median scaled |p(root)|
        CONVERGE_RATE       = 0.98     % required err==0 rate
    end

    % =====================================================================
    %  Per-test result arrays — populated in TestMethodSetup
    % =====================================================================
    properties
        Degrees         % (1 × N_TESTS) degree of each random polynomial
        ScaledResiduals % (1 × N_TESTS) max|p(root_k)| / max(1,|r|^N)
        ErrFlags        % (1 × N_TESTS) err return value from newton_real
        RootCounts      % (1 × N_TESTS) numel(roots_out) from newton_real
        Passed          % (1 × N_TESTS) logical — scaled residual < RESIDUAL_TOL
    end

    % =====================================================================
    %  Shared setup — runs once before each test method
    % =====================================================================
    methods (TestMethodSetup)
        function generateAndRunPolynomials(testCase)
            rng(testCase.RANDOM_SEED, 'twister');

            n = testCase.N_TESTS;
            testCase.Degrees         = randi([testCase.MIN_DEG, testCase.MAX_DEG], 1, n);
            testCase.ScaledResiduals = inf(1, n);
            testCase.ErrFlags        = zeros(1, n);
            testCase.RootCounts      = zeros(1, n);
            testCase.Passed          = false(1, n);

            for k = 1 : n
                deg    = testCase.Degrees(k);
                coeffs = [1.0, randn(1, deg)];   % monic, random remaining terms

                [r, err] = newton_real(coeffs);

                res    = abs(polyval(coeffs, r));
                scale  = max(1, max(abs(r))^deg);
                testCase.ScaledResiduals(k) = max(res) / scale;
                testCase.ErrFlags(k)        = err;
                testCase.RootCounts(k)      = numel(r);
                testCase.Passed(k)          = (testCase.ScaledResiduals(k) < testCase.RESIDUAL_TOL);
            end
        end
    end

    % =====================================================================
    %  Test methods — assertions only, no logic
    % =====================================================================
    methods (Test)

        % -----------------------------------------------------------------
        %  Root-count sanity: newton_real must return exactly deg roots
        % -----------------------------------------------------------------
        function testRootCountMatchesDegree(testCase)
            testCase.verifyEqual(testCase.RootCounts, testCase.Degrees, ...
                'newton_real must return exactly degree(p) roots for every polynomial.');
        end

        % -----------------------------------------------------------------
        %  Overall success rate across all degrees
        % -----------------------------------------------------------------
        function testOverallSuccessRate(testCase)
            rate = mean(testCase.Passed);
            testCase.verifyGreaterThanOrEqual(rate, testCase.OVERALL_PASS_RATE, ...
                sprintf('Overall success rate %.1f%% is below the %.0f%% requirement.', ...
                100*rate, 100*testCase.OVERALL_PASS_RATE));
        end

        % -----------------------------------------------------------------
        %  Per-degree-band success rates
        % -----------------------------------------------------------------
        function testLowDegreeSuccessRate(testCase)
            mask = (testCase.Degrees >= 4) & (testCase.Degrees <= 10);
            rate = mean(testCase.Passed(mask));
            testCase.verifyGreaterThanOrEqual(rate, testCase.LOW_DEG_PASS_RATE, ...
                sprintf('Low-degree (4–10) success rate %.1f%% is below %.0f%%.', ...
                100*rate, 100*testCase.LOW_DEG_PASS_RATE));
        end

        function testMidDegreeSuccessRate(testCase)
            mask = (testCase.Degrees >= 11) & (testCase.Degrees <= 18);
            rate = mean(testCase.Passed(mask));
            testCase.verifyGreaterThanOrEqual(rate, testCase.MID_DEG_PASS_RATE, ...
                sprintf('Mid-degree (11–18) success rate %.1f%% is below %.0f%%.', ...
                100*rate, 100*testCase.MID_DEG_PASS_RATE));
        end

        function testHighDegreeSuccessRate(testCase)
            mask = (testCase.Degrees >= 19) & (testCase.Degrees <= 25);
            rate = mean(testCase.Passed(mask));
            testCase.verifyGreaterThanOrEqual(rate, testCase.HIGH_DEG_PASS_RATE, ...
                sprintf('High-degree (19–25) success rate %.1f%% is below %.0f%%.', ...
                100*rate, 100*testCase.HIGH_DEG_PASS_RATE));
        end

        % -----------------------------------------------------------------
        %  Residual quality for passing tests
        % -----------------------------------------------------------------
        function testPassingResidualsAreBelowTolerance(testCase)
            passingResiduals = testCase.ScaledResiduals(testCase.Passed);
            testCase.verifyTrue(all(passingResiduals < testCase.RESIDUAL_TOL), ...
                'Every test flagged as passing must have scaled max|p(root)| < RESIDUAL_TOL.');
        end

        function testMedianResidualIsSmall(testCase)
            passingResiduals = testCase.ScaledResiduals(testCase.Passed);
            testCase.verifyLessThan(median(passingResiduals), testCase.MEDIAN_RESIDUAL_MAX, ...
                sprintf('Median scaled residual of passing tests exceeds %.0e.', ...
                testCase.MEDIAN_RESIDUAL_MAX));
        end

        % -----------------------------------------------------------------
        %  Convergence flag: most polynomials should reach err = 0
        % -----------------------------------------------------------------
        function testMostPolynomialsConverge(testCase)
            rate = mean(testCase.ErrFlags == 0);
            testCase.verifyGreaterThanOrEqual(rate, testCase.CONVERGE_RATE, ...
                sprintf('Convergence rate %.1f%% is below the %.0f%% requirement.', ...
                100*rate, 100*testCase.CONVERGE_RATE));
        end

    end
end
