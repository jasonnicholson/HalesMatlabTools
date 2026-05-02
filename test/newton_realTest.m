classdef newton_realTest < matlab.unittest.TestCase

  methods (TestClassSetup)
    % Shared setup for the entire test class
  end

  methods (TestMethodSetup)
    % Setup for each test
  end

  methods (Test)
    % Test linear equation: x - 2 = 0 -> root = 2
    function testLinear(testCase)
      coeff = [1, -2];
      [res, err] = newton_real(coeff);
      testCase.verifyEqual(err, 0, "Error should be 0");
      testCase.verifyEqual(numel(res), 1, "Should have 1 root");
      testCase.verifyEqual(real(res(1)), 2.0, 'AbsTol', 1e-14, "Root should be 2");
    end

    % Test quadratic with real roots: x^2 - 5x + 6 = 0 -> roots = [2, 3]
    function testQuadraticRealRoots(testCase)
      coeff = [1, -5, 6];
      [res, err] = newton_real(coeff);
      testCase.verifyEqual(err, 0, "Error should be 0");
      testCase.verifyEqual(numel(res), 2, "Should have 2 roots");
      
      % Sort roots for consistent comparison
      roots_sorted = sort(real(res));
      testCase.verifyEqual(roots_sorted(1), 2.0, 'AbsTol', 1e-14, "First root should be 2");
      testCase.verifyEqual(roots_sorted(2), 3.0, 'AbsTol', 1e-14, "Second root should be 3");
    end

    % Test quadratic with complex roots: x^2 + 1 = 0 -> roots = [i, -i]
    function testQuadraticComplexRoots(testCase)
      coeff = [1, 0, 1];
      [res, err] = newton_real(coeff);
      testCase.verifyEqual(err, 0, "Error should be 0");
      testCase.verifyEqual(numel(res), 2, "Should have 2 roots");
      
      % Roots should be ±i
      testCase.verifyEqual(abs(imag(res(1))), 1.0, 'AbsTol', 1e-14, "Root should have imag part = ±1");
      testCase.verifyEqual(abs(imag(res(2))), 1.0, 'AbsTol', 1e-14, "Root should have imag part = ±1");
      testCase.verifyEqual(real(res(1)), 0.0, 'AbsTol', 1e-14, "Real part should be 0");
      testCase.verifyEqual(real(res(2)), 0.0, 'AbsTol', 1e-14, "Real part should be 0");
    end

    % Test cubic with one real root: x^3 - 1 = 0 -> root = 1 (and 2 complex roots)
    function testCubic(testCase)
      coeff = [1, 0, 0, -1];
      [res, err] = newton_real(coeff);
      testCase.verifyEqual(err, 0, "Error should be 0");
      testCase.verifyEqual(numel(res), 3, "Should have 3 roots");
      
      % Find the real root (should be 1)
      real_roots = res(abs(imag(res)) < 1e-14);
      testCase.verifyGreaterThan(length(real_roots), 0, "Should have at least one real root");
      testCase.verifyEqual(real(real_roots(1)), 1.0, 'AbsTol', 1e-12, "Real root should be 1");
    end

    % Test polynomial with zero roots: x^2(x - 1) = 0 -> roots = [0, 0, 1]
    function testPolynomialWithZeroRoots(testCase)
      coeff = [1, -1, 0, 0];
      [res, err] = newton_real(coeff);
      testCase.verifyEqual(err, 0, "Error should be 0");
      testCase.verifyEqual(length(res), 3, "Should have 3 roots");
      
      % Count zeros
      zero_count = sum(abs(res) < 1e-10);
      testCase.verifyGreaterThanOrEqual(zero_count, 2, "Should have at least 2 zero roots");
    end

    % Test simple polynomial: (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6
    function testCubicWithRealRoots(testCase)
      coeff = [1, -6, 11, -6];
      [res, err] = newton_real(coeff);
      testCase.verifyEqual(err, 0, "Error should be 0");
      testCase.verifyEqual(length(res), 3, "Should have 3 roots");
      
      % All roots should be real
      imag_parts = abs(imag(res));
      testCase.verifyLessThan(max(imag_parts), 1e-10, "All roots should be real");
      
      % Sort and verify
      roots_sorted = sort(real(res));
      testCase.verifyEqual(roots_sorted(1), 1.0, 'AbsTol', 1e-10, "First root should be 1");
      testCase.verifyEqual(roots_sorted(2), 2.0, 'AbsTol', 1e-10, "Second root should be 2");
      testCase.verifyEqual(roots_sorted(3), 3.0, 'AbsTol', 1e-10, "Third root should be 3");
    end

    % Test quartic with all complex roots: x^4 + 1 = 0
    % Triggers complex deflation since algorithm finds complex conjugate pairs
    function testQuarticAllComplex(testCase)
      coeff = [1, 0, 0, 0, 1];
      [res, err] = newton_real(coeff);
      testCase.verifyEqual(err, 0, "Error should be 0");
      testCase.verifyEqual(length(res), 4, "Should have 4 roots");
      
      % All roots should be complex (non-zero imaginary parts)
      complex_count = sum(abs(imag(res)) > 1e-8);
      testCase.verifyEqual(complex_count, 4, "All 4 roots should be complex");
      
      % Verify roots come in conjugate pairs and satisfy equation
      for k = 1:4
        z = res(k);
        val = z^4 + 1;
        testCase.verifyEqual(abs(val), 0.0, 'AbsTol', 1e-12, ...
          sprintf("Root %d should satisfy z^4 + 1 = 0", k));
      end
    end

    % Test quartic with mixed roots: (x-1)(x-2)(x^2+1) = x^4 - 3x^3 + 3x^2 - 3x + 2
    % Forces algorithm to find both real and complex conjugate pair roots
    function testQuarticMixedRoots(testCase)
      coeff = [1, -3, 3, -3, 2];
      [res, err] = newton_real(coeff);
      testCase.verifyEqual(err, 0, "Error should be 0");
      testCase.verifyEqual(length(res), 4, "Should have 4 roots");
      
      % Should have 2 real roots and 2 complex roots
      real_roots = res(abs(imag(res)) < 1e-10);
      complex_roots = res(abs(imag(res)) > 1e-10);
      testCase.verifyEqual(length(real_roots), 2, "Should have 2 real roots");
      testCase.verifyEqual(length(complex_roots), 2, "Should have 2 complex roots");
      
      % Verify real roots are 1 and 2
      real_sorted = sort(real(real_roots));
      testCase.verifyEqual(real_sorted(1), 1.0, 'AbsTol', 1e-12);
      testCase.verifyEqual(real_sorted(2), 2.0, 'AbsTol', 1e-12);
    end

    % Test quintic: x^5 - 2x^4 + x - 2 = (x-2)(x^4+1)
    % This forces the algorithm through complex deflation multiple times
    function testQuinticComplex(testCase)
      % Coefficients for (x-2)(x^4+1) = x^5 - 2x^4 + x - 2
      coeff = [1, -2, 0, 0, 1, -2];
      [res, err] = newton_real(coeff);
      testCase.verifyEqual(err, 0, "Error should be 0");
      testCase.verifyEqual(length(res), 5, "Should have 5 roots");
      
      % Verify one real root (x=2)
      real_roots = res(abs(imag(res)) < 1e-10);
      testCase.verifyEqual(length(real_roots), 1, "Should have 1 real root");
      testCase.verifyEqual(real(real_roots(1)), 2.0, 'AbsTol', 1e-12);
      
      % Verify roots satisfy polynomial equation
      for k = 1:5
        z = res(k);
        val = z^5 - 2*z^4 + z - 2;
        testCase.verifyEqual(abs(val), 0.0, 'AbsTol', 1e-10, ...
          sprintf("Root %d should satisfy polynomial", k));
      end
    end

    % Test sextic with all complex roots except one real: x^6 - 1 = (x-1)(x+1)(x^4+x^2+1)
    % Root x=1 is real, x=-1 is real, others are complex
    function testSextic(testCase)
      coeff = [1, 0, 0, 0, 0, 0, -1];
      [res, err] = newton_real(coeff);
      testCase.verifyEqual(err, 0, "Error should be 0");
      testCase.verifyEqual(length(res), 6, "Should have 6 roots");
      
      % Verify all roots satisfy z^6 - 1 = 0
      for k = 1:6
        z = res(k);
        val = z^6 - 1;
        testCase.verifyEqual(abs(val), 0.0, 'AbsTol', 1e-12, ...
          sprintf("Root %d should satisfy z^6 - 1 = 0", k));
      end
      
      % Should have at least the real roots 1 and -1
      real_roots = res(abs(imag(res)) < 1e-10);
      testCase.verifyEqual(length(real_roots), 2, "Should have 2 real roots (1 and -1)");
    end

    % --- err output -------------------------------------------------

    % err output is 0 for a well-conditioned cubic
    function testErrZeroOnNormalConvergence(testCase)
      [~, err] = newton_real([1, -6, 11, -6]);
      testCase.verifyEqual(err, 0, "err should be 0 for well-conditioned input");
    end

    % Single-output call must work without capturing err
    function testErrOutputPresentSingleOutput(testCase)
      res = newton_real([1, -5, 6]);
      testCase.verifyEqual(numel(res), 2, "Single-output call should return 2 roots");
    end

    % --- Edge cases -------------------------------------------------

    % Constant polynomial p(x) = 5 has no roots
    function testConstantPolynomial(testCase)
      res = newton_real([5]);
      testCase.verifyEmpty(res, "Constant polynomial should return zero roots");
    end

    % p(x) = x^3 — three roots all at zero
    function testAllZeroRoots(testCase)
      res = newton_real([1, 0, 0, 0]);
      testCase.verifyEqual(numel(res), 3, "Should have 3 roots");
      testCase.verifyLessThan(max(abs(res)), 1e-7, "All roots should be at zero");
    end

    % --- Linear -----------------------------------------------------

    % p(x) = x + 4 -> root = -4
    function testLinearNegativeRoot(testCase)
      [res, err] = newton_real([1, 4]);
      testCase.verifyEqual(err, 0, "Error should be 0");
      testCase.verifyEqual(real(res(1)), -4.0, 'AbsTol', 1e-14, "Root should be -4");
      testCase.verifyEqual(imag(res(1)),  0.0, 'AbsTol', 1e-14, "Root should be real");
    end

    % --- Quadratic --------------------------------------------------

    % p(x) = (x-3)^2 = x^2 - 6x + 9 — repeated root at 3
    function testQuadraticRepeatedRoot(testCase)
      [res, err] = newton_real([1, -6, 9]);
      testCase.verifyEqual(err, 0, "Error should be 0");
      testCase.verifyEqual(numel(res), 2, "Should have 2 roots");
      testCase.verifyEqual(real(res(1)), 3.0, 'AbsTol', 1e-5, "Root should be 3");
      testCase.verifyEqual(real(res(2)), 3.0, 'AbsTol', 1e-5, "Root should be 3");
    end

    % p(x) = x^2 - 4 -> roots ±2 (no linear term — depressed quadratic)
    function testQuadraticDepressed(testCase)
      [res, err] = newton_real([1, 0, -4]);
      testCase.verifyEqual(err, 0, "Error should be 0");
      roots_sorted = sort(real(res));
      testCase.verifyEqual(roots_sorted(1), -2.0, 'AbsTol', 1e-14, "Root should be -2");
      testCase.verifyEqual(roots_sorted(2),  2.0, 'AbsTol', 1e-14, "Root should be 2");
    end

    % p(x) = -x^2 + 4 -> roots ±2 (negative leading coefficient)
    function testNegativeLeadingCoefficient(testCase)
      [res, err] = newton_real([-1, 0, 4]);
      testCase.verifyEqual(err, 0, "Error should be 0");
      roots_sorted = sort(real(res));
      testCase.verifyEqual(roots_sorted(1), -2.0, 'AbsTol', 1e-14, "Root should be -2");
      testCase.verifyEqual(roots_sorted(2),  2.0, 'AbsTol', 1e-14, "Root should be 2");
    end

    % --- Cubic ------------------------------------------------------

    % p(x) = (x-1)(x^2+1) -> one real root (1) and two complex roots (±i)
    function testCubicOneRealTwoComplex(testCase)
      coeff = [1, -1, 1, -1];
      [res, err] = newton_real(coeff);
      testCase.verifyEqual(err, 0, "Error should be 0");
      testCase.verifyEqual(numel(res), 3, "Should have 3 roots");
      is_real = abs(imag(res)) < 1e-10;
      testCase.verifyEqual(sum(is_real), 1, "Expected exactly one real root");
      testCase.verifyEqual(real(res(is_real)), 1.0, 'AbsTol', 1e-12, "Real root should be 1");
      complex_roots = res(~is_real);
      testCase.verifyEqual(abs(complex_roots(1)), 1.0, 'AbsTol', 1e-12, "Complex root modulus should be 1");
    end

    % --- Quartic ----------------------------------------------------

    % p(x) = (x^2+1)(x^2+4) -> roots ±i, ±2i (two complex conjugate pairs)
    function testQuarticTwoComplexPairs(testCase)
      coeff = [1, 0, 5, 0, 4];
      [res, err] = newton_real(coeff);
      testCase.verifyEqual(err, 0, "Error should be 0");
      testCase.verifyEqual(numel(res), 4, "Should have 4 roots");
      moduli = sort(abs(res));
      testCase.verifyEqual(moduli(1), 1.0, 'AbsTol', 1e-12, "Two roots should have modulus 1");
      testCase.verifyEqual(moduli(2), 1.0, 'AbsTol', 1e-12);
      testCase.verifyEqual(moduli(3), 2.0, 'AbsTol', 1e-12, "Two roots should have modulus 2");
      testCase.verifyEqual(moduli(4), 2.0, 'AbsTol', 1e-12);
      testCase.verifyLessThan(max(abs(real(res))), 1e-10, "All roots should be purely imaginary");
    end

    % p(x) = (x-1)(x-2)(x-3)(x-4) -> four distinct real roots
    function testQuarticAllRealRoots(testCase)
      coeff = [1, -10, 35, -50, 24];
      [res, err] = newton_real(coeff);
      testCase.verifyEqual(err, 0, "Error should be 0");
      roots_sorted = sort(real(res));
      testCase.verifyEqual(roots_sorted(1), 1.0, 'AbsTol', 1e-12);
      testCase.verifyEqual(roots_sorted(2), 2.0, 'AbsTol', 1e-12);
      testCase.verifyEqual(roots_sorted(3), 3.0, 'AbsTol', 1e-12);
      testCase.verifyEqual(roots_sorted(4), 4.0, 'AbsTol', 1e-12);
    end

    % --- Roots at zero ----------------------------------------------

    % p(x) = x^2(x-5) -> roots 0, 0, 5
    function testSingleRootAtZero(testCase)
      coeff = [1, -5, 0, 0];
      [res, err] = newton_real(coeff);
      testCase.verifyEqual(err, 0, "Error should be 0");
      roots_sorted = sort(real(res));
      testCase.verifyEqual(roots_sorted(1), 0.0, 'AbsTol', 1e-10);
      testCase.verifyEqual(roots_sorted(2), 0.0, 'AbsTol', 1e-10);
      testCase.verifyEqual(roots_sorted(3), 5.0, 'AbsTol', 1e-12);
    end

    % --- Higher degree ----------------------------------------------

    % p(x) = (x-1)(x-2)(x-3)(x-4)(x-5) -> five distinct real roots
    function testDegree5AllReal(testCase)
      true_roots = (1:5).';
      coeff = poly(true_roots);
      [res, err] = newton_real(coeff);
      testCase.verifyEqual(err, 0, "Error should be 0");
      roots_sorted = sort(real(res));
      testCase.verifyEqual(roots_sorted, true_roots, 'AbsTol', 1e-6);
    end

    % p(x) = (x-1)(x-3)(x^2+1)(x^2+4) -> 2 real + 4 complex roots
    function testDegree6MixedRoots(testCase)
      coeff = conv(poly([1, 3]), conv([1, 0, 1], [1, 0, 4]));
      [res, err] = newton_real(coeff);
      testCase.verifyEqual(err, 0, "Error should be 0");
      testCase.verifyEqual(numel(res), 6, "Should have 6 roots");
      for k = 1:6
        testCase.verifyEqual(abs(polyval(coeff, res(k))), 0.0, 'AbsTol', 1e-8, ...
          sprintf("Root %d residual too large", k));
      end
    end

    % --- Residual checks --------------------------------------------

    % p(x) = (x-1)^4 — near-multiple root; verify by residual
    function testResidualRepeatedRoot(testCase)
      coeff = [1, -4, 6, -4, 1];
      [res, err] = newton_real(coeff);
      testCase.verifyLessThanOrEqual(err, 0, "err should be non-positive");
      for k = 1:4
        testCase.verifyEqual(abs(polyval(coeff, res(k))), 0.0, 'AbsTol', 1e-5, ...
          sprintf("Root %d residual too large", k));
      end
    end

    % 8 Chebyshev-spaced roots in (-1,1) — tests robustness to clustered roots
    function testResidualChebyshevSpaced(testCase)
      true_roots = cos((1:2:15) * pi / 16).';
      coeff = real(poly(true_roots));
      [res, err] = newton_real(coeff);
      testCase.verifyLessThanOrEqual(err, 0, "err should be non-positive");
      for k = 1:8
        testCase.verifyEqual(abs(polyval(coeff, res(k))), 0.0, 'AbsTol', 1e-5, ...
          sprintf("Root %d residual too large", k));
      end
    end

    % --- Consistency with built-in roots() -------------------------

    % Compare sorted real parts of quartic against MATLAB companion-matrix roots
    function testConsistencyBuiltinRootsQuartic(testCase)
      coeff   = [1, -10, 35, -50, 24];   % (x-1)(x-2)(x-3)(x-4)
      r_ref   = sort(real(roots(coeff)));
      [r_mine, err] = newton_real(coeff);
      testCase.verifyEqual(err, 0, "Error should be 0");
      testCase.verifyEqual(sort(real(r_mine)), r_ref, 'AbsTol', 1e-8);
    end

    % Compare root moduli of degree-4 complex-pair polynomial
    function testConsistencyBuiltinRootsComplexPair(testCase)
      coeff  = [1, 0, 5, 0, 4];
      r_ref  = sort(abs(roots(coeff)));
      [r_mine, err] = newton_real(coeff);
      testCase.verifyEqual(err, 0, "Error should be 0");
      testCase.verifyEqual(sort(abs(r_mine)), r_ref, 'AbsTol', 1e-7);
    end

  end

end