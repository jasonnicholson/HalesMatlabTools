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
