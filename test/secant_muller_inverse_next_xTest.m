classdef secant_muller_inverse_next_xTest < matlab.unittest.TestCase

  methods (TestClassSetup)
    function addSourcePath(~)
      thisFile = mfilename('fullpath');
      testDir = fileparts(thisFile);
      sourceDir = fullfile(fileparts(testDir), 'source');
      addpath(sourceDir);
      addpath(testDir);
    end
  end

  methods (Test)
    function testVectorsParity(testCase)
      vectors = secant_muller_inverse_test_vectors();

      for k = 1:numel(vectors)
        v = vectors(k);
        out = secant_muller_inverse_next_x( ...
          v.points, v.targetY, v.stepSize, v.method, ...
          struct('denominatorTol', 1e-12, 'duplicateXTol', 1e-12));

        switch v.name
          case "secant-basic"
            testCase.verifyTrue(out.ok);
            testCase.verifyEqual(out.xNext, 3.157162792479945, 'AbsTol', 1e-12);
            testCase.verifyEqual(string(out.methodUsed), "secant");
          case "secant-next"
            testCase.verifyTrue(out.ok);
            testCase.verifyEqual(out.xNext, 3.14154625560721, 'AbsTol', 1e-12);
            testCase.verifyEqual(string(out.methodUsed), "secant");
          case "muller-near-pi"
            testCase.verifyTrue(out.ok);
            testCase.verifyEqual(out.xNext, 3.141592670639796, 'AbsTol', 1e-12);
            testCase.verifyEqual(string(out.methodUsed), "muller");
          case "step-size-half"
            testCase.verifyTrue(out.ok);
            testCase.verifyEqual(out.xNext, 3.1493545240536047, 'AbsTol', 1e-12);
            testCase.verifyEqual(string(out.methodUsed), "secant");
          case "zero-denominator"
            testCase.verifyFalse(out.ok);
            testCase.verifyEqual(string(out.blockedReason), "zero_denominator");
            testCase.verifyEqual(string(out.methodUsed), "secant");
          case "repeated-x"
            testCase.verifyFalse(out.ok);
            testCase.verifyEqual(string(out.blockedReason), "repeated_x");
            testCase.verifyEqual(string(out.methodUsed), "muller");
          otherwise
            testCase.assertFail("Unhandled vector case: " + v.name);
        end
      end
    end

    function testYEditTriggersNextRowCompute(testCase)
      data = baseData();
      data{1, 4} = 4.0;
      data{1, 5} = -0.513604990615856;
      data{2, 4} = 3.0;
      data{2, 5} = 1.28224001611973;
      data{2, 2} = 'secant';
      data{2, 3} = 1.0;

      [outData, eventOut] = secant_muller_inverse_apply_edit(data, 2, 5, 1.0);

      testCase.verifyTrue(eventOut.triggered);
      testCase.verifyTrue(eventOut.ok);
      testCase.verifyEqual(eventOut.triggerTargetRow, 3);
      testCase.verifyEqual(outData{3, 4}, 3.157162792479945, 'AbsTol', 1e-12);
    end

    function testStepSizeEditTriggersNextRowCompute(testCase)
      data = baseData();
      data{1, 4} = 3.0;
      data{1, 5} = 1.28224001611973;
      data{2, 4} = 3.1571627925;
      data{2, 5} = 0.968860980423343;
      data{2, 2} = 'secant';
      data{2, 3} = 0.5;

      [outData, eventOut] = secant_muller_inverse_apply_edit(data, 2, 3, 1.0);

      testCase.verifyTrue(eventOut.triggered);
      testCase.verifyTrue(eventOut.ok);
      testCase.verifyEqual(outData{3, 4}, 3.1493545240536047, 'AbsTol', 1e-12);
    end

    function testMethodEditTriggersNextRowCompute(testCase)
      data = baseData();
      data{1, 4} = 3.0;
      data{1, 5} = 1.28224001611973;
      data{2, 4} = 3.1571627925;
      data{2, 5} = 0.968860980423343;
      data{3, 4} = 3.1415462556;
      data{3, 5} = 1.00009279600125;
      data{3, 2} = 'muller';
      data{3, 3} = 1.0;

      [outData, eventOut] = secant_muller_inverse_apply_edit(data, 3, 2, 1.0);

      testCase.verifyTrue(eventOut.triggered);
      testCase.verifyTrue(eventOut.ok);
      testCase.verifyEqual(outData{4, 4}, 3.141592670639796, 'AbsTol', 1e-12);
      testCase.verifyEqual(string(eventOut.methodUsed), "muller");
    end

    function testNonTriggerEditDoesNotCompute(testCase)
      data = baseData();
      [outData, eventOut] = secant_muller_inverse_apply_edit(data, 2, 4, 1.0);

      testCase.verifyFalse(eventOut.triggered);
      testCase.verifyEqual(string(eventOut.blockedReason), "non_trigger_field");
      testCase.verifyEqual(outData, data);
    end

    function testOutOfRangeEditIsBlocked(testCase)
      data = baseData();
      [~, eventOut] = secant_muller_inverse_apply_edit(data, 99, 5, 1.0);

      testCase.verifyFalse(eventOut.triggered);
      testCase.verifyEqual(string(eventOut.blockedReason), "out_of_range");
    end

    function testLastRowTriggerHasMissingTargetRow(testCase)
      data = baseData();
      lastRow = size(data, 1);
      [~, eventOut] = secant_muller_inverse_apply_edit(data, lastRow, 5, 1.0);

      testCase.verifyTrue(eventOut.triggered);
      testCase.verifyEqual(string(eventOut.blockedReason), "target_row_missing");
      testCase.verifyFalse(eventOut.ok);
    end

    function testTriggeredEditWithInsufficientPointsPropagatesBlock(testCase)
      data = baseData();
      data{1, 4} = 3.0;
      data{1, 5} = 1.2;

      [~, eventOut] = secant_muller_inverse_apply_edit(data, 1, 5, 1.0);

      testCase.verifyTrue(eventOut.triggered);
      testCase.verifyEqual(string(eventOut.blockedReason), "insufficient_points");
      testCase.verifyFalse(eventOut.ok);
    end

    function testAutoMethodSecantWithTwoPoints(testCase)
      pts = [4.0, -0.513604990615856; 3.0, 1.28224001611973];
      out = secant_muller_inverse_next_x(pts, 1.0, 1.0, "", struct());

      testCase.verifyTrue(out.ok);
      testCase.verifyEqual(string(out.methodUsed), "secant");
    end

    function testAutoMethodMullerWithThreePoints(testCase)
      pts = [
        3.0, 1.28224001611973;
        3.1571627925, 0.968860980423343;
        3.1415462556, 1.00009279600125
      ];
      out = secant_muller_inverse_next_x(pts, 1.0, 1.0, "", struct());

      testCase.verifyTrue(out.ok);
      testCase.verifyEqual(string(out.methodUsed), "muller");
    end

    function testExplicitMullerWithOnlyTwoPointsBlocked(testCase)
      pts = [4.0, -0.5; 3.0, 1.2];
      out = secant_muller_inverse_next_x(pts, 1.0, 1.0, "muller", struct());

      testCase.verifyFalse(out.ok);
      testCase.verifyEqual(string(out.blockedReason), "insufficient_points");
    end

    function testInvalidTargetYBlocked(testCase)
      pts = [4.0, -0.5; 3.0, 1.2];
      out = secant_muller_inverse_next_x(pts, NaN, 1.0, "secant", struct());

      testCase.verifyFalse(out.ok);
      testCase.verifyEqual(string(out.blockedReason), "invalid_number");
    end

    function testInvalidStepSizeBlocked(testCase)
      pts = [4.0, -0.5; 3.0, 1.2];
      out = secant_muller_inverse_next_x(pts, 1.0, Inf, "secant", struct());

      testCase.verifyFalse(out.ok);
      testCase.verifyEqual(string(out.blockedReason), "invalid_number");
    end

    function testComplexMullerStepBlocked(testCase)
      pts = [
        -1.0, 1.0;
        0.0, 0.0;
        1.0, 1.0
      ];
      out = secant_muller_inverse_next_x(pts, -1.0, 1.0, "muller", struct());

      testCase.verifyFalse(out.ok);
      testCase.verifyEqual(string(out.blockedReason), "complex_muller_step");
    end

    function testSecantRepeatedXBlocked(testCase)
      pts = [
        2.0, 0.1;
        2.0, 0.2
      ];
      out = secant_muller_inverse_next_x(pts, 0.0, 1.0, "secant", struct());

      testCase.verifyFalse(out.ok);
      testCase.verifyEqual(string(out.blockedReason), "repeated_x");
    end

    function testMullerSmallSpacingBlockedByDenominatorTol(testCase)
      pts = [
        0.0, 1.0;
        1e-13, 2.0;
        2e-13, 3.0
      ];
      out = secant_muller_inverse_next_x( ...
        pts, 0.0, 1.0, "muller", ...
        struct('denominatorTol', 1e-12, 'duplicateXTol', 0.0));

      testCase.verifyFalse(out.ok);
      testCase.verifyEqual(string(out.blockedReason), "zero_denominator");
    end

    function testMullerEEqualsZeroBlocked(testCase)
      pts = [
        0.0, 0.0;
        1.0, 0.0;
        2.0, 0.0
      ];
      out = secant_muller_inverse_next_x(pts, 0.0, 1.0, "muller", struct());

      testCase.verifyFalse(out.ok);
      testCase.verifyEqual(string(out.blockedReason), "zero_denominator");
    end

    function testDefaultsPathWithMinimalArguments(testCase)
      pts = [4.0, -0.513604990615856; 3.0, 1.28224001611973];
      out = secant_muller_inverse_next_x(pts);

      testCase.verifyTrue(out.ok);
      testCase.verifyEqual(string(out.methodUsed), "secant");
    end

    function testEmptyPointsBlocked(testCase)
      out = secant_muller_inverse_next_x(zeros(0, 2), 1.0, 1.0, "secant", struct());
      testCase.verifyFalse(out.ok);
      testCase.verifyEqual(string(out.blockedReason), "insufficient_points");
    end

    function testApplyEditDefaultsAndOptionPassthrough(testCase)
      data = baseData();
      data{1, 4} = 4.0;
      data{1, 5} = -0.513604990615856;
      data{2, 4} = 3.0;
      data{2, 5} = 1.28224001611973;
      data{2, 2} = 'secant';
      data{2, 3} = 1.0;

      [outData, eventOut] = secant_muller_inverse_apply_edit( ...
        data, 2, 5, 1.0, struct('denominatorTol', 1e-12, 'duplicateXTol', 1e-12));

      testCase.verifyTrue(eventOut.ok);
      testCase.verifyTrue(isfinite(outData{3, 4}));

      % Call with omitted targetY/opts to hit argument defaults.
      [~, eventOutDefault] = secant_muller_inverse_apply_edit(data, 2, 5);
      testCase.verifyTrue(eventOutDefault.triggered);
    end

    function testApplyEditStringAndFallbackParsingBranches(testCase)
      data = baseData();
      data{1, 4} = 4.0;
      data{1, 5} = -0.513604990615856;
      data{2, 4} = 3.0;
      data{2, 5} = 1.28224001611973;
      data{2, 2} = 7;      % non-char -> method fallback path in cell_to_default
      data{2, 3} = "0.5";  % string parse success in scalar_from_cell

      [outData1, eventOut1] = secant_muller_inverse_apply_edit(data, 2, 3, 1.0);
      testCase.verifyTrue(eventOut1.ok);
      testCase.verifyEqual(outData1{3, 4}, 3.0785813962399725, 'AbsTol', 1e-12);

      data{2, 3} = "not-a-number"; % string parse failure -> fallback branch
      [outData2, eventOut2] = secant_muller_inverse_apply_edit(data, 2, 3, 1.0);
      testCase.verifyTrue(eventOut2.ok);
      testCase.verifyEqual(outData2{3, 4}, 3.157162792479945, 'AbsTol', 1e-12);

      data{2, 3} = []; % isempty branch
      [outData3, eventOut3] = secant_muller_inverse_apply_edit(data, 2, 3, 1.0);
      testCase.verifyTrue(eventOut3.ok);
      testCase.verifyEqual(outData3{3, 4}, 3.157162792479945, 'AbsTol', 1e-12);
    end

    function testAppBuildAndCallbackExecution(testCase)
      app = secant_muller_inverse_app('NumRows', 5, 'Visible', 'off');
      cleanup = onCleanup(@() delete(app.Figure)); %#ok<NASGU>

      testCase.verifySize(app.Table.Data, [5, 5]);
      testCase.verifyEqual(app.TargetYField.Value, 1.0);

      % Seed two points then invoke table callback for row-2 y edit.
      data = app.Table.Data;
      data{1, 4} = 4.0;
      data{1, 5} = -0.513604990615856;
      data{2, 4} = 3.0;
      data{2, 5} = 1.28224001611973;
      app.Table.Data = data;

      evt = struct('Indices', [2, 5]);
      cb = app.Table.CellEditCallback;
      cb(app.Table, evt);

      testCase.verifyTrue(isfinite(app.Table.Data{3, 4}));
      testCase.verifyNotEmpty(app.StatusLabel.Text);
    end

    function testAppCallbackNonTriggerStatus(testCase)
      app = secant_muller_inverse_app('NumRows', 5, 'Visible', 'off');
      cleanup = onCleanup(@() delete(app.Figure)); %#ok<NASGU>

      evt = struct('Indices', [2, 4]); % x edit is non-trigger
      cb = app.Table.CellEditCallback;
      cb(app.Table, evt);

      testCase.verifyEqual(string(app.StatusLabel.Text), "Edit accepted (no recompute trigger).");
    end

    function testAppCallbackBlockedStatus(testCase)
      app = secant_muller_inverse_app('NumRows', 5, 'Visible', 'off');
      cleanup = onCleanup(@() delete(app.Figure)); %#ok<NASGU>

      % Trigger edit with insufficient points so next-x is blocked.
      data = app.Table.Data;
      data{1, 4} = 3.0;
      data{1, 5} = 1.2;
      app.Table.Data = data;

      evt = struct('Indices', [1, 5]);
      cb = app.Table.CellEditCallback;
      cb(app.Table, evt);

      testCase.verifyTrue(contains(string(app.StatusLabel.Text), "compute blocked"));
    end
  end
end

function data = baseData()
cfg = secant_muller_inverse_contract();
nRows = 6;
data = cell(nRows, 5);
for i = 1:nRows
  data{i, 1} = i;
  if i >= 3
    data{i, 2} = char(cfg.defaultMethodByPointCount.threePoints);
  else
    data{i, 2} = char(cfg.defaultMethodByPointCount.twoPoints);
  end
  data{i, 3} = 1.0;
  data{i, 4} = NaN;
  data{i, 5} = NaN;
end
end
