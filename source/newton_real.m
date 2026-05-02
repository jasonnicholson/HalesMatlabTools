function [roots_out, err] = newton_real(coeff)
%NEWTON_REAL  Find all roots of a real-coefficient polynomial.
%
%   Implements Madsen's two-stage modified Newton method (1973), translated
%   from the C++ reference implementation by Future Team Aps (2002) and
%   rewritten in idiomatic MATLAB. Uses native complex arithmetic, MATLAB's
%   polyval for Horner evaluation, deconv for deflation, and vectorised
%   helper computations.
%
% Syntax
%   roots_out        = newton_real(coeff)
%   [roots_out, err] = newton_real(coeff)
%
% Inputs
%   coeff  - (1×N+1) or (N+1×1) real coefficient vector, leading-coefficient
%            first — the same convention used by polyval and roots:
%                p(x) = coeff(1)·x^N + coeff(2)·x^(N-1) + … + coeff(N+1)
%
% Outputs
%   roots_out - N×1 complex vector of all polynomial roots.
%   err       - Integer convergence flag (mirrors the C++ return value):
%                 0        all roots converged within MAXITER iterations.
%               < 0        |err| roots hit the MAXITER limit; the
%                          corresponding entries in roots_out are the best
%                          approximations available but may be inaccurate.
%
% Algorithm
%   Roots are extracted one (or two, for complex conjugate pairs) at a time
%   via forward deflation.  For each root the two-stage Newton method is
%   applied:
%
%     Stage 1 – Search mode
%       The Newton step is accepted only if it reduces |p(z)|.  If the
%       step is too large it is halved; if reducing it still helps the
%       iterate is extended further along the Newton direction, giving
%       super-linear convergence to multiple roots.
%
%     Stage 2 – Newton mode
%       Plain Newton iteration, entered when a Kantorowitz-type heuristic
%       indicates quadratic convergence has begun.  Any failed step causes
%       an immediate return to Stage 1.
%
%     Convergence
%       Iteration stops when the step size is below machine resolution,
%       |p(z)| ≈ 0 (Adam's rounding-error bound in Stage 2), or MAXITER
%       (50) is reached.
%
%     Root classification (real-coefficient optimisation)
%       After converging to z, the residual at real(z) is compared to that
%       at z.  If real(z) achieves a lower residual it is accepted as a
%       real root and the polynomial is deflated by (x − real(z)).
%       Otherwise z and conj(z) are recorded as a conjugate pair and the
%       polynomial is deflated by the quadratic factor (x−z)(x−conj(z)).
%
%   After deflation reduces the degree to ≤ 2 the remaining factor is
%   solved by exact formulas.
%
% Notes
%   • Coefficient convention is leading-first (same as polyval/roots),
%     which differs from the Madsen-Reid Fortran report (constant-first).
%   • polyval (Horner's rule) and deconv use MATLAB's optimised BLAS.
%   • deconv performs polynomial long division; tiny remainders due to
%     finite-precision roots are expected and discarded.
%   • Maximum 50 Newton iterations per root extraction (MAXITER).
%   • err mirrors the int return value of the C++ newton_real(): starts at
%     0 and is decremented by 1 for each root that exhausts MAXITER.
%
% Example
%   % Cubic with three real roots
%   p = [1, -6, 11, -6];                    % (x−1)(x−2)(x−3)
%   [r, err] = newton_real(p)         % r ≈ [1; 2; 3], err = 0
%
%   % Quadratic with complex roots
%   p = [1, 0, 1];                   % x^2 + 1
%   r = newton_real(p)        % ≈ [1i; −1i]
%
% References
%   Madsen, K. (1973). A root-finding algorithm based on Newton's method.
%   BIT Numerical Mathematics, 13(1), 71–75.
%
%   Madsen, K. & Reid, J. (1975). Fortran subroutines for finding
%   polynomial zeros. AERE-R 7986, Harwell.
%
% See also: roots, polyval, deconv

 arguments
      coeff (1,:) double {mustBeReal, mustBeNonempty, mustBeFinite}
    end

    % --- Algorithm constants -------------------------------------------
    MAXITER  = 50;
    RADIX    = 2;
    MANT_DIG = 53;                       % IEEE 754 double: 53 mantissa bits
    SQRT_EPS = RADIX^(-MANT_DIG / 2);   % ≈ 1.49e-8; step-size guard

    % --- Input normalisation ------------------------------------------

    a  = double(coeff(:).');             % row vector; a(1) = leading coeff
    N  = numel(a) - 1;                   % polynomial degree

    roots_out = complex(zeros(N, 1));
    err = 0;                             % convergence flag (0 = all OK)
    rp = N;                              % root pointer — filled top-down

    % Strip zero constant terms; each zero constant means a root at x = 0
    while N > 0 && a(N + 1) == 0
        roots_out(rp) = 0;
        N  = N  - 1;
        rp = rp - 1;
    end
    a = a(1:N + 1);

    if N <= 0, return; end

    % ==================================================================
    % Main deflation loop: extract roots until degree ≤ 2
    % ==================================================================
    while N > 2

        % Derivative coefficient vector (vectorised).
        % For p(x) = Σ a(k+1)·x^(N-k), k=0..N, the derivative is
        % p'(x) = Σ (N-k)·a(k+1)·x^(N-k-1), k=0..N-1.
        % In leading-first storage: da(k) = (N-k+1)·a(k), k=1..N.
        da = a(1:N) .* (N:-1:1);       % leading-first; length N = degree of p'

        % --- Initial guess: Cauchy-type lower bound on root modulus -----
        u_start = cauchy_bound(a, N);

        % Seed on the real axis, direction from −a(N+1)/a(N).
        % This mimics the C++ starting-direction heuristic.
        if a(N) ~= 0
            seed_sign = sign(-a(N + 1) / a(N));
        else
            seed_sign = 1;
        end
        if seed_sign == 0, seed_sign = 1; end

        z   = seed_sign * u_start * 0.5;  % real scalar; MATLAB promotes later
        dz  = z;
        z0  = 0;
        fz0 = a(N);                        % proxy: p'(0) / (N-1)! ≈ a(N)
        [f,  fz]  = peval(a,  z);
        f0  = 2 * abs(a(N + 1))^2;         % loose initial bound for f0
        ff  = f0;
        r0  = 2.5 * u_start;
        r   = abs(dz);
        % Convergence threshold: stop when |p(z)|² falls below the machine
        % rounding-noise floor estimated from f0 — mirrors the C++ formula
        %   eps = 4·n²·f0·RADIX^(−2·MANT_DIG)
        % Using f > 0 instead (as the old code did) causes the loop to grind
        % through all MAXITER iterations even when the root is already at full
        % machine precision, producing spurious err < 0 flags.
        eps_iter = 4 * N^2 * f0 * RADIX^(-2*MANT_DIG);

        stage1  = true;
        itercnt = 0;

        % --- Newton iteration loop ------------------------------------
        while itercnt < MAXITER && f > eps_iter && ...
              (real(z) + real(dz) ~= real(z) || imag(z) + imag(dz) ~= imag(z))

            itercnt = itercnt + 1;

            [u2, fz1] = peval(da, z);      % evaluate p'(z); u2 = |p'(z)|²

            if u2 == 0
                % True saddle point: rotate and enlarge step to escape
                dz = rotate_step(dz, 5.0);
            else
                % Newton step: p(z)/p'(z), computed using conjugate division
                % to stay close to the real-arithmetic form of the C++ code.
                dz = (fz * conj(fz1)) / u2;  % equiv to fz / fz1

                % Kantorowitz-type stage decision ----------------------
                %   Stage 1 needed when a finite-difference proxy f2 exceeds
                %   the Kantorowitz threshold  2|p||p''| ≤ |p'|²|Δz|.
                %   The heuristic mirrors the C++: f2/u2 > u2/(4f) ≡ stage1.
                fwz   = fz0 - fz1;   % Δp' ≈ change in derivative
                wdiff = z0  - z;     % Δz
                if abs(wdiff) > 0
                    f2     = abs(fwz)^2 / abs(wdiff)^2;
                    stage1 = (f2 / u2 > u2 / (4 * f)) || (f ~= ff);
                else
                    stage1 = true;
                end

                r = abs(dz);
                if r > r0
                    dz = rotate_step(dz, r0 / r);   % cap step magnitude
                end
                r0 = r * 5.0;
            end

            % Save current iterate for the next stage decision
            z0 = z;  f0 = f;  fz0 = fz;

            % --- Candidate step ----------------------------------------
            z = z0 - dz;
            [f, fz] = peval(a, z);
            ff = f;

            if stage1
                % Line search: halve or extend dz to guarantee |p| decreases.
                % A halving strategy is used when the initial step overshoots
                % (f > f0); an extension strategy is used when it undershoots.
                wz   = z;
                div2 = (f > f0);    % true → step was too large; halve it
                for i = 1:N
                    if div2
                        dz = dz * 0.5;
                        wz = z0 - dz;
                    else
                        wz = wz - dz;        % extend one more Newton length
                    end
                    [fw, fwz_val] = peval(a, wz);
                    if fw >= f, break; end   % no further improvement; stop
                    f = fw;  fz = fwz_val;  z = wz;
                    if div2 && i == 2
                        % Two halvings insufficient: rotate once and retry
                        dz = rotate_step(dz, 0.5);
                        z  = z0 - dz;
                        [f, fz] = peval(a, z);
                        break;
                    end
                end
            else
                % Stage 2: tighten convergence threshold to Adam's local
                % rounding-error bound (C++ upperbound() call).  Iteration
                % stops as soon as |p(z)|² is within the evaluation noise.
                eps_iter = adam_bound(a, N, z);
            end

            % Domain rounding-error guard: if the step has shrunk below
            % machine resolution without improving the residual, alter the
            % direction; if still stuck, break (accept current z0).
            if r < abs(z) * SQRT_EPS && f >= f0
                z  = z0;
                dz = rotate_step(dz, 0.5);
                if real(z) + real(dz) == real(z) && ...
                   imag(z) + imag(dz) == imag(z)
                    break;
                end
            end

        end % Newton iteration loop

        if itercnt >= MAXITER
            err = err - 1;   % flag that this root did not converge
        end

        % --- Classify root: real or complex conjugate pair ---------------
        % Exploit real-coefficient symmetry: if projecting z to the real
        % axis reduces the residual, the root is real.
        [f_re, ~] = peval(a, real(z));
        if f_re <= f
            % Real root — deflate by linear factor (x − z_real)
            roots_out(rp) = real(z);
            a  = deconv(a, [1, -real(z)]);
            N  = N - 1;
            rp = rp - 1;
        else
            % Complex conjugate pair — deflate by (x−z)(x−conj(z))
            roots_out(rp)     = z;
            roots_out(rp - 1) = conj(z);
            a  = deconv(a, [1, -2*real(z), abs(z)^2]);
            N  = N - 2;
            rp = rp - 2;
        end

    end % deflation while loop

    % --- Solve the remaining linear or quadratic factor ------------------
    if N > 0
        roots_out(1:N) = solve_low_degree(a(1:N + 1), N);
    end

end


% =========================================================================
%  Local helper functions
% =========================================================================

function u = cauchy_bound(a, n)
%CAUCHY_BOUND  Cauchy lower bound on the smallest root modulus.
%
%   Computes  u = min_{0≤k<n, a(k+1)≠0}  (|a(n+1)| / |a(k+1)|)^{1/(n-k)}
%
%   This is a standard lower bound for the modulus of the smallest root
%   (Cauchy, 1829).  The starting iterate for Newton's method is placed at
%   u/2 so that the initial point lies inside the smallest root's modulus.
%   Vectorised: all valid coefficient ratios are evaluated simultaneously.
%
%   a is in leading-first order; a(n+1) is the constant term.

    an = abs(a(n + 1));                % |constant term|
    if an == 0, u = 1; return; end     % zero constant → root at origin already stripped

    r    = log(an);
    k    = 1:n;                        % non-constant coefficient indices
    pows = n + 1 - k;                  % Cauchy exponents: n, n-1, …, 1

    % Always include leading coefficient (index 1); skip exact zeros
    valid    = (a(k) ~= 0);
    valid(1) = true;                   % leading coeff must be nonzero for a valid poly

    sel = k(valid);                    % selected indices into a
    u   = min(exp((r - log(abs(a(sel)))) ./ pows(valid)));
end


function [mod2, val] = peval(c, z)
%PEVAL  Evaluate polynomial c at z; return (|p(z)|², p(z)).
%
%   c is in leading-first order (same as MATLAB polyval).
%   MATLAB's polyval applies Horner's rule internally using optimised BLAS.
%   For complex z the evaluation is performed in complex arithmetic.

    val  = polyval(c, z);
    mod2 = real(val)^2 + imag(val)^2;   % |val|² without a redundant sqrt
end


function dz = rotate_step(dz, m)
%ROTATE_STEP  Rotate step vector by ~53.13° and scale by m.
%
%   Applies the same fixed Givens-like rotation as C++ alterdirection():
%       [ 0.6  -0.8 ] [ Re(dz) ]
%       [ 0.8   0.6 ] [ Im(dz) ] * m
%
%   cos(53.13°) ≈ 0.6, sin(53.13°) ≈ 0.8.  Rotating the Newton step
%   escapes saddle points and avoids oscillation on quadratic-factor axes.
%   Scaling by m < 1 reduces step magnitude; m > 1 enlarges it.

    re = real(dz);  im = imag(dz);
    dz = complex(re * 0.6 - im * 0.8, ...
                 re * 0.8 + im * 0.6) * m;
end


function e2 = adam_bound(a, n, z)
%ADAM_BOUND  Adam's rounding-error upper bound on |p(z)|², squared.
%
%   Computes a local estimate of the squared rounding error in evaluating
%   p at z via the real-arithmetic Horner recurrence (Adam, 1983; used
%   in the Madsen-Reid Fortran library as the Stage-2 convergence guard).
%   Returns e² so it can be compared directly against |p(z)|² = f.
%
%   a is leading-first; n = degree (length(a)−1).
%   Mirrors C++ upperbound() in Newton_Real.cpp.

    p = -2 * real(z);
    q = real(z)^2 + imag(z)^2;
    u = sqrt(q);
    s = 0;  r = a(1);  e = abs(r) * (3.5 / 4.5);
    for i = 2:n
        t = a(i) - p*r - q*s;
        s = r;  r = t;
        e = u*e + abs(t);
    end
    t = a(n+1) + real(z)*r - q*s;
    e = u*e + abs(t);
    % eps(1) = RADIX^(−MANT_DIG+1) = 2^(−52); matches C++ _DBL_RADIX^(−DBL_MANT_DIG+1)
    e = (9*e - 7*(abs(t) + abs(r)*u) + abs(real(z))*abs(r)*2) * eps(1);
    e2 = e^2;
end


function res = solve_low_degree(a, n)
%SOLVE_LOW_DEGREE  Exact roots of a linear (n=1) or quadratic (n=2).
%
%   Uses the numerically stable formulas ported from the C++ reference:
%     • Linear:    root = −a(2)/a(1)
%     • Quadratic (a(2)=0):   roots = ±sqrt(−a(3)/a(1))
%     • Quadratic (general):  Vieta's formula avoids catastrophic cancellation
%                             in the discriminant form.
%
%   a is in leading-first order: a(1)·x^n + a(2)·x^(n-1) + … + a(n+1).

    res = complex(zeros(n, 1));

    if n == 1
        % Linear: a(1)·x + a(2) = 0
        res(1) = -a(2) / a(1);

    else  % n == 2 quadratic: a(1)·x² + a(2)·x + a(3) = 0
        if a(2) == 0
            % Depressed quadratic: x = ±sqrt(−a(3)/a(1))
            r = -a(3) / a(1);
            if r < 0
                res(1) =  1i * sqrt(-r);
                res(2) = -1i * sqrt(-r);
            else
                res(1) =  sqrt(r);
                res(2) = -sqrt(r);
            end
        else
            % General quadratic via scaled discriminant  d = 1 − 4a₁a₃/a₂²
            % This form avoids loss of significance when a₂² >> 4a₁a₃.
            disc = 1 - 4 * a(1) * a(3) / a(2)^2;
            if disc < 0
                % Complex conjugate pair
                re_part = -a(2) / (2 * a(1));
                im_part =  a(2) * sqrt(-disc) / (2 * a(1));
                res(1) = complex(re_part,  im_part);
                res(2) = complex(re_part, -im_part);
            else
                % Two real roots — use Vieta's product to avoid cancellation
                res(1) = (-1 - sqrt(disc)) * a(2) / (2 * a(1));
                res(2) = a(3) / (a(1) * real(res(1)));  % Vieta: r1·r2 = a3/a1
            end
        end
    end
end
