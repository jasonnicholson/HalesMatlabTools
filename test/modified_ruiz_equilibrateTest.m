classdef modified_ruiz_equilibrateTest < matlab.unittest.TestCase

  methods (Test)

    function testAcceptsRowVectorsForQlu(testCase)
      P = speye(2);
      A = sparse([1, 0; 0, 1; 1, 1]);
      q = [1, 2];
      l = [-1, 0, 0];
      u = [2, 3, 4];

      [~, ~, qOut, lOut, uOut, scaling] = modified_ruiz_equilibrate(P, A, q, l, u, 1);

      testCase.verifySize(qOut, [2, 1]);
      testCase.verifySize(lOut, [3, 1]);
      testCase.verifySize(uOut, [3, 1]);
      testCase.verifyEqual(numel(scaling.D), 2);
      testCase.verifyEqual(numel(scaling.E), 3);
    end

    function testRejectsAMatrixWithWrongColumnCount(testCase)
      P = speye(2);
      A = sparse([1, 0, 0; 0, 1, 1]);
      q = [1; 2];
      l = [-1; -1];
      u = [1; 1];

      testCase.verifyError(@() modified_ruiz_equilibrate(P, A, q, l, u), ...
        'modified_ruiz_equilibrate:DimensionMismatch');
    end

    function testRejectsMismatchedVectorLengths(testCase)
      P = speye(2);
      A = sparse([1, 0; 0, 1; 1, 1]);
      q = [1; 2; 3];
      l = [-1; -1; -1];
      u = [1; 1; 1];

      testCase.verifyError(@() modified_ruiz_equilibrate(P, A, q, l, u), ...
        'modified_ruiz_equilibrate:DimensionMismatch');
    end

    function testRejectsInvalidScalingOptionOrder(testCase)
      P = speye(2);
      A = sparse([1, 0; 0, 1]);
      q = [1; 2];
      l = [0; 0];
      u = [1; 1];

      testCase.verifyError(@() modified_ruiz_equilibrate(P, A, q, l, u, 1, ...
        minScaling=10, maxScaling=1), ...
        'modified_ruiz_equilibrate:InvalidScalingBounds');
    end

  end
end
