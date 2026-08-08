function cfg = secant_muller_inverse_contract()
%SECANT_MULLER_INVERSE_CONTRACT Shared behavior contract for inverse solver UI.
% Returns constants and rules used by both MATLAB and web implementations.

cfg = struct();

% Table schema and defaults
cfg.columns = ["Iteration", "Method", "Step Size", "x", "y"];
cfg.methodOptions = ["secant", "muller"];
cfg.defaultMethodByPointCount = struct("twoPoints", "secant", "threePoints", "muller");
cfg.defaultStepSize = 1.0;
cfg.stepSizeMin = 0.0;

% Numeric behavior
cfg.targetYDefault = 1.0;
cfg.displaySignificantDigits = 8;
cfg.tolerance = struct( ...
  "denominator", 1e-12, ...
  "duplicateX", 1e-12);

% UX behavior gates
cfg.dropdownPolicy = "disable-invalid-choices";
cfg.stepSizeScope = "per-row";
cfg.numericDisplayPolicy = "fixed-significant-digits";

% Blocked compute reasons (stable IDs for callbacks/UI messaging)
cfg.blockedReason = struct( ...
  "insufficient_points", "insufficient_points", ...
  "repeated_x", "repeated_x", ...
  "zero_denominator", "zero_denominator", ...
  "complex_muller_step", "complex_muller_step", ...
  "invalid_number", "invalid_number");

% Event contract
cfg.event = struct( ...
  "trigger", "y-or-step-size-or-method-edit-complete", ...
  "triggerFields", ["y", "Step Size", "Method"], ...
  "triggerSourceRow", "edited-row", ...
  "triggerTargetRow", "edited-row-plus-one", ...
  "action", "compute-next-x");
end
