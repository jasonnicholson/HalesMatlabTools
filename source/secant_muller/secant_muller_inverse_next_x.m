function out = secant_muller_inverse_next_x(points, targetY, stepSize, method, opts)
  %SECANT_MULLER_INVERSE_NEXT_X Compute next x from recent (x,y) points.
  % points: Nx2 numeric matrix [x,y]. Uses latest rows.
  % targetY: scalar desired y.
  % stepSize: scalar multiplier on base delta x.
  % method: "secant" or "muller"; if empty, auto-select by available points.
  % opts: struct with optional fields:
  %   - denominatorTol (default 1e-12)
  %   - duplicateXTol (default 1e-12)
  %
  % Returns struct:
  %   ok, methodUsed, xNext, deltaBase, blockedReason

  arguments
    points (:,2) double
    targetY (1,1) double = 1.0
    stepSize (1,1) double = 1.0
    method string = ""
    opts struct = struct()
  end
  
  denomTol = get_opt(opts, "denominatorTol", 1e-12);
  dupTol = get_opt(opts, "duplicateXTol", 1e-12);

  out = struct( ...
    "ok", false, ...
    "methodUsed", "", ...
    "xNext", NaN, ...
    "deltaBase", NaN, ...
    "blockedReason", "");

  % Basic input validation
  if ~isnumeric(points) || size(points,2) ~= 2 || isempty(points)
    out.blockedReason = "insufficient_points";
    return;
  end
  if ~isscalar(targetY) || ~isfinite(targetY)
    out.blockedReason = "invalid_number";
    return;
  end
  if ~isscalar(stepSize) || ~isfinite(stepSize)
    out.blockedReason = "invalid_number";
    return;
  end

  finiteMask = all(isfinite(points), 2);
  valid = points(finiteMask, :);
  n = size(valid, 1);

  if n < 2
    out.blockedReason = "insufficient_points";
    return;
  end

  methodUsed = normalize_method(method, n);
  if methodUsed == ""
    out.blockedReason = "insufficient_points";
    return;
  end

  out.methodUsed = methodUsed;

  if methodUsed == "secant"
    p = valid(max(1, n-1):n, :);
    [deltaBase, reason] = secant_delta(p, targetY, denomTol, dupTol);
  else
    p = valid(max(1, n-2):n, :);
    [deltaBase, reason] = muller_delta(p, targetY, denomTol, dupTol);
  end

  if reason ~= ""
    out.blockedReason = reason;
    return;
  end

  xCurrent = p(end, 1);
  out.deltaBase = deltaBase;
  out.xNext = xCurrent + stepSize * deltaBase;
  out.ok = true;
end

function [delta, reason] = secant_delta(p, targetY, denomTol, dupTol)
  reason = "";
  delta = NaN;

  x1 = p(1,1); y1 = p(1,2);
  x2 = p(2,1); y2 = p(2,2);

  if abs(x2 - x1) <= dupTol
    reason = "repeated_x";
    return;
  end

  g1 = y1 - targetY;
  g2 = y2 - targetY;
  denom = g2 - g1;

  if abs(denom) <= denomTol
    reason = "zero_denominator";
    return;
  end

  delta = -g2 * (x2 - x1) / denom;
end

function [delta, reason] = muller_delta(p, targetY, denomTol, dupTol)
  reason = "";
  delta = NaN;

  x0 = p(1,1); y0 = p(1,2);
  x1 = p(2,1); y1 = p(2,2);
  x2 = p(3,1); y2 = p(3,2);

  if abs(x1 - x0) <= dupTol || abs(x2 - x1) <= dupTol || abs(x2 - x0) <= dupTol
    reason = "repeated_x";
    return;
  end

  g0 = y0 - targetY;
  g1 = y1 - targetY;
  g2 = y2 - targetY;

  h0 = x1 - x0;
  h1 = x2 - x1;

  if abs(h0) <= denomTol || abs(h1) <= denomTol || abs(h0 + h1) <= denomTol
    reason = "zero_denominator";
    return;
  end

  d0 = (g1 - g0) / h0;
  d1 = (g2 - g1) / h1;
  a = (d1 - d0) / (h1 + h0);
  b = d1 + h1 * a;
  c = g2;

  D = sqrt(b^2 - 4 * a * c);
  Eplus = b + D;
  Eminus = b - D;

  if abs(Eplus) >= abs(Eminus)
    E = Eplus;
  else
    E = Eminus;
  end

  if abs(E) <= denomTol
    reason = "zero_denominator";
    return;
  end

  deltaComplex = -2 * c / E;
  if abs(imag(deltaComplex)) > denomTol
    reason = "complex_muller_step";
    return;
  end

  delta = real(deltaComplex);
end

function value = get_opt(opts, key, defaultValue)
  if isfield(opts, key)
    value = opts.(key);
  else
    value = defaultValue;
  end
end

function methodUsed = normalize_method(method, nPoints)
  m = lower(string(method));
  if m == ""
    if nPoints >= 3
      methodUsed = "muller";
    else
      methodUsed = "secant";
    end
    return;
  end

  if m == "secant" && nPoints >= 2
    methodUsed = "secant";
    return;
  end

  if m == "muller" && nPoints >= 3
    methodUsed = "muller";
    return;
  end

  methodUsed = "";
end
