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
  end

end