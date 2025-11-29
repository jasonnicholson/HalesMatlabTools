function [res, err] = newton_real(coeff)
  % NEWTON_REAL Solve n degree polynomial using Newton's (Madsen) method ::
  %
  %   [res, err] = newton_real(coeff)
  %
  % Inputs:
  %   coeff: Vector of real coefficients [a_n, a_{n-1}, ..., a_0]
  %          corresponding to      P(x) = a_n*x^n + ...    + a_0.
  %
  % Outputs:
  %   res: Vector of complex roots.
  %   err: Error code (0 if successful, < 0 if max iterations reached).
  %
  % Description:
  %   This function implements Newton's method (as described by Madsen and Reid)
  %   to find all roots of a real-coefficient polynomial. It iteratively refines
  %   guesses for the roots and performs deflation to reduce the polynomial degree
  %   after each root is found. The method handles both real and complex roots.
  %
  % Example:
  %
  %
  % References:
  %   - Madsen, Kaj, and John K. Reid. "Fortran subroutines for finding polynomial zeros." (1975).

  arguments
    coeff (:,1) double
  end

  % Constants
  MAX_ITERATIONS = 50;
  DBL_MANT_DIG = 53; % Standard double precision

  % Determine degree n
  n = length(coeff) - 1;

  a = coeff;

  % Initialize result array (roots)
  res = zeros(1, n);

  err = 0;

  % Remove trailing zeros (roots at 0).
  % TODO vectorize this loop
  while n >= 1 && a(n+1) == 0.0
    res(n) = 0;
    n = n - 1;
  end

  % Resize a to match current n
  a = a(1:n+1);

  % Allocate derivative array
  a1 = zeros(1, n);

  while n > 2
    % Calculate coefficients of f'(x)
    for i = 0:(n-1)
      a1(i+1) = a(i+1) * (n - i);
    end

    u = startpoint(n, a);
    z0 = complex(0, 0);
    f0 = 2.0 * a(n+1) * a(n+1);
    ff = f0;
    fz0 = complex(a(n)); % a[n-1]

    % Initial guess z
    if a(n) == 0.0
      z_real = 1.0;
    else
      z_real = -a(n+1) / a(n);
    end

    z = complex(z_real, 0);

    if abs(real(z)) == 0
      factor = 0;
    else
      factor = real(z) / abs(real(z));
    end
    z = complex(factor * u * 0.5, 0);

    dz = z;
    [f, fz] = feval_poly(n, a, z);

    r0 = 2.5 * u;
    r = abs(dz);
    eps_val = 4 * n * n * f0 * pow2(-DBL_MANT_DIG * 2);

    % Start iteration
    itercnt = 0;
    stage1 = true;

    % Main Iteration loop
    while ( (real(z) + real(dz) ~= real(z)) || (imag(z) + imag(dz) ~= imag(z)) ) && ...
        (f > eps_val) && (itercnt < MAX_ITERATIONS)

      itercnt = itercnt + 1;

      [u_val, fz1] = feval_poly(n - 1, a1, z);

      if u_val == 0.0 % True saddle point
        dz = alterdirection(dz, 5.0);
      else
        % Newton step
        % dz = fz / fz1
        dz = fz / fz1;

        % Which stage are we on
        fwz = fz0 - fz1;
        wz = z0 - z;

        if abs(wz) == 0
          f2 = 0;
        else
          f2 = (real(fwz)^2 + imag(fwz)^2) / (real(wz)^2 + imag(wz)^2);
        end

        stage1 = (f2/u_val > u_val/f/4) || (f ~= ff);
        r = abs(dz);
        if r > r0
          dz = alterdirection(dz, r0 / r);
        end
        r0 = r * 5.0;
      end

      z0 = z;
      f0 = f;
      fz0 = fz;

      % Inner retry loop (equivalent to goto iter2 logic)
      retry_iter2 = true;
      while retry_iter2
        retry_iter2 = false; % Default to not retrying unless condition met

        z = z0 - dz;
        [f, fz] = feval_poly(n, a, z);
        ff = f;

        if stage1
          wz = z;
          div2 = (f > f0);

          for i = 1:n
            if div2
              dz = dz * 0.5;
              wz = z0 - dz;
            else
              wz = wz - dz;
            end

            [fw, fwz] = feval_poly(n, a, wz);

            if fw >= f
              break;
            end

            f = fw;
            fz = fwz;
            z = wz;

            if div2 && i == 2
              dz = alterdirection(dz, 0.5);
              z = z0 - dz;
              [f, fz] = feval_poly(n, a, z);
              break;
            end
          end
        else
          % Calculate upper bound of errors using Adam's test
          eps_val = upperbound(n, a, z, DBL_MANT_DIG);
        end

        % Domain rounding errors check
        if ( r < abs(z) * pow2(-DBL_MANT_DIG/2) ) && (f >= f0)
          z = z0;
          dz = alterdirection(dz, 0.5);
          if z + dz ~= z
            retry_iter2 = true; % Re-run inner loop
          end
        end
      end
    end

    if itercnt >= MAX_ITERATIONS
      err = err - 1;
    end

    z0_real = complex(real(z), 0.0);
    [f_z0, ~] = feval_poly(n, a, z0_real);

    if f_z0 <= f
      % Real root
      res(n) = complex(real(z), 0);
      [n, a] = realdeflation(n, a, real(z));
    else
      % Complex root
      res(n) = z;
      res(n - 1) = complex(real(z), -imag(z));
      [n, a] = complexdeflation(n, a, z);
    end
  end

  % Solve remaining linear or quadratic
  res = quadratic(n, a, res);
end

% --------------------------------------------------------------------------
% Helper Functions
% --------------------------------------------------------------------------

function res = quadratic(n, a, res)
  % Solve linear or quadratic equation

  if n == 2
    if a(2) == 0 % a[1] == 0
      r = -a(3) / a(1); % -a[2]/a[0]
      if r < 0
        res(2) = complex(0, sqrt(-r));
        res(3) = complex(0, -imag(res(2)));
      else
        res(2) = complex(sqrt(r), 0);
        res(3) = complex(-real(res(2)), 0);
      end
    else
      r = 1 - 4 * a(1) * a(3) / (a(2) * a(2));
      if r < 0
        val1 = complex( -a(2)/(2*a(1)), a(2)*sqrt(-r)/(2*a(1)) );
        res(2) = val1;
        res(3) = complex(real(val1), -imag(val1));
      else
        val1 = complex( (-1 - sqrt(r)) * a(2) / (2 * a(1)), 0 );
        res(2) = val1;
        res(3) = complex( a(3) / (a(1) * real(val1)), 0 );
      end
    end
  elseif n == 1
    res(2) = complex( -a(2) / a(1), 0 );
  end
end

function [val, fz] = feval_poly(n, a, z)
  % Performed function evaluation. Horners algorithm.
  % Returns |f(z)|^2 as val, and f(z) as fz.

  p = -2.0 * real(z);
  q = real(z)^2 + imag(z)^2;
  s = 0;
  r = a(1); % a[0]

  for i = 1:(n-1)
    t = a(i+1) - p * r - q * s;
    s = r;
    r = t;
  end

  fz_real = a(n+1) + real(z) * r - q * s;
  fz_imag = imag(z) * r;
  fz = complex(fz_real, fz_imag);

  val = real(fz)^2 + imag(fz)^2;
end

function min_val = startpoint(n, a)
  % Determine starting point
  r = log( abs( a(n+1) ) ); % a[n]
  min_val = exp( ( r - log( abs( a(1) ) ) ) / n ); % a[0]

  for i = 1:(n-1)
    if a(i+1) ~= 0
      u = exp( ( r - log( abs( a(i+1) ) ) ) / ( n - i ) );
      if u < min_val
        min_val = u;
      end
    end
  end
end

function e_sq = upperbound(n, a, z, DBL_MANT_DIG)
  % Calculate a upper bound for the rounding errors (Adam's test)

  p = -2.0 * real(z);
  q = real(z)^2 + imag(z)^2;
  u = sqrt(q);
  s = 0.0;
  r = a(1);
  e = abs(r) * (3.5 / 4.5);

  for i = 1:(n-1)
    t = a(i+1) - p * r - q * s;
    s = r;
    r = t;
    e = u * e + abs(t);
  end

  t = a(n+1) + real(z) * r - q * s;
  e = u * e + abs(t);
  e = ( 9.0 * e - 7.0 * ( abs(t) + abs(r) * u ) + ...
    abs(real(z)) * abs(r) * 2.0 ) * pow2(-DBL_MANT_DIG+1);

  e_sq = e * e;
end

function dz = alterdirection(dz, m)
  x = ( real(dz) * 0.6 - imag(dz) * 0.8 ) * m;
  y = ( real(dz) * 0.8 + imag(dz) * 0.6 ) * m;
  dz = complex(x, y);
end

function [n_new, a] = realdeflation(n, a, x)
  % Real root forward deflation
  r = 0;
  for i = 0:(n-1)
    r = r * x + a(i+1);
    a(i+1) = r;
  end
  n_new = n - 1;
end

function [n_new, a] = complexdeflation(n, a, z)
  % Complex root forward deflation
  r = -2.0 * real(z);
  u = real(z)^2 + imag(z)^2;

  a(2) = a(2) - r * a(1);
  for i = 2:(n-2)
    a(i+1) = a(i+1) - r * a(i) - u * a(i-1);
  end
  n_new = n - 2;
end