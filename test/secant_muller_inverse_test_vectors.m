function vectors = secant_muller_inverse_test_vectors()
%SECANT_MULLER_INVERSE_TEST_VECTORS Deterministic vectors for parity checks.

vectors = struct();

vectors(1).name = "secant-basic";
vectors(1).targetY = 1.0;
vectors(1).stepSize = 1.0;
vectors(1).method = "secant";
vectors(1).points = [
  4.0, -0.513604990615856;
  3.0, 1.28224001611973
];

vectors(2).name = "secant-next";
vectors(2).targetY = 1.0;
vectors(2).stepSize = 1.0;
vectors(2).method = "secant";
vectors(2).points = [
  3.0, 1.28224001611973;
  3.1571627925, 0.968860980423343
];

vectors(3).name = "muller-near-pi";
vectors(3).targetY = 1.0;
vectors(3).stepSize = 1.0;
vectors(3).method = "muller";
vectors(3).points = [
  3.0, 1.28224001611973;
  3.1571627925, 0.968860980423343;
  3.1415462556, 1.00009279600125
];

vectors(4).name = "step-size-half";
vectors(4).targetY = 1.0;
vectors(4).stepSize = 0.5;
vectors(4).method = "secant";
vectors(4).points = [
  3.0, 1.28224001611973;
  3.1571627925, 0.968860980423343
];

vectors(5).name = "zero-denominator";
vectors(5).targetY = 1.0;
vectors(5).stepSize = 1.0;
vectors(5).method = "secant";
vectors(5).points = [
  2.0, 1.25;
  3.0, 1.25
];

vectors(6).name = "repeated-x";
vectors(6).targetY = 1.0;
vectors(6).stepSize = 1.0;
vectors(6).method = "muller";
vectors(6).points = [
  3.0, 1.28224001611973;
  3.0, 1.28224001611973;
  3.1571627925, 0.968860980423343
];
end
