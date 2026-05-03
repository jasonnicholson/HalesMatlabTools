function [roots_out, err] = newton_real(coeff)
%NEWTON_REAL  Find all roots of a real-coefficient polynomial.
%
%   Implements Madsen's two-stage modified Newton method (1973). Roots are
%   extracted one (real) or two (complex conjugate pair) at a time via
%   forward deflation. The polynomial is evaluated with a paired-real
%   Horner scheme that exploits real coefficients at complex arguments,
%   and deflation uses O(n) synthetic division.
%
% Syntax
%   roots_out        = newton_real(coeff)
%   [roots_out, err] = newton_real(coeff)
%
% Inputs
%   coeff  - real coefficient vector, leading-coefficient first — the same
%            convention used by polyval and roots:
%                p(x) = coeff(1)·x^N + coeff(2)·x^(N-1) + … + coeff(N+1)
%
% Outputs
%   roots_out - N×1 complex vector of all polynomial roots.
%   err       - Integer convergence flag (mirrors the original C++ return).
%                 0   all roots converged within MAXITER iterations.
%               < 0   |err| roots hit the MAXITER limit; the corresponding
%                     entries in roots_out are best approximations.
%
% Algorithm
%   Each root is found by Madsen's two-stage Newton iteration:
%
%     Stage 1 — Search mode
%       The Newton step is accepted only if it reduces |p(z)|. If the
%       step is too large it is halved; if reducing helps, the iterate is
%       extended further along the Newton direction (super-linear
%       convergence to multiple roots).
%
%     Stage 2 — Newton mode
%       Plain Newton iteration, entered when a Kantorowitz-type heuristic
%       indicates quadratic convergence has begun. A failed step returns
%       to Stage 1.
%
%     Convergence
%       Iteration stops when the step size is below machine resolution,
%       |p(z)|^2 falls below Adam's local rounding-error bound, or
%       MAXITER (50) is reached.
%
%     Root classification
%       After converging to z, the residual at real(z) is compared to that
%       at z. If real(z) achieves a no-worse residual it is accepted as a
%       real root and the polynomial is deflated by (x − real(z)).
%       Otherwise z and conj(z) are recorded and the polynomial is deflated
%       by the quadratic factor (x−z)(x−conj(z)).
%
%   After deflation reduces the degree to ≤ 2 the remaining factor is
%   solved with closed-form formulas.
%
% Implementation notes
%   • Coefficient convention is leading-first (same as polyval/roots),
%     which differs from the Madsen-Reid Fortran report (constant-first).
%   • peval() uses a paired-real Horner recurrence, which performs only
%     real-real multiplications until the final p(z) reconstruction —
%     about half the multiplications of a generic complex Horner.
%   • Real deflation uses synthetic division by (x − r); complex deflation
%     uses synthetic division by (x^2 − 2·Re(z)·x + |z|^2). Both run in
%     O(n) and overwrite the working coefficient buffer in place.
%   • Maximum 50 Newton iterations per root extraction (MAXITER).
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
    SQRT_EPS = RADIX^(-MANT_DIG / 2);    % ≈ 1.49e-8; step-size guard

    % --- Input normalisation ------------------------------------------
    a  = double(coeff(:).');             % row vector; a(1) = leading coeff
    N  = numel(a) - 1;                   % polynomial degree

    roots_out = complex(zeros(N, 1));
    err = 0;                             % convergence flag (0 = all OK)
    rp = N;                              % root pointer — filled top-down

    if N <= 0, return; end

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
        da = a(1:N) .* (N:-1:1);       % length N = degree of p'

        % --- Initial guess: Cauchy-type lower bound on root modulus -----
        u_start = cauchy_bound(a, N);

        % Seed on the real axis, direction from −a(N+1)/a(N).
        if a(N) ~= 0
            seed_sign = sign(-a(N + 1) / a(N));
        else
            seed_sign = 1;
        end
        if seed_sign == 0, seed_sign = 1; end

        z   = complex(seed_sign * u_start * 0.5, 0);
        dz  = z;
        z0  = complex(0);
        fz0 = complex(a(N));               % proxy: p'(0) / (N-1)! ≈ a(N)
        [f,  fz]  = peval(N, a, z);
        f0  = 2 * abs(a(N + 1))^2;
        ff  = f0;
        r0  = 2.5 * u_start;

        % Convergence threshold: stop when |p(z)|² falls below the machine
        % rounding-noise floor estimated from f0 — mirrors the C++ formula
        %   eps = 4·n²·f0·RADIX^(−2·MANT_DIG)
        eps_iter = 4 * N^2 * f0 * RADIX^(-2*MANT_DIG);

        stage1  = true;
        itercnt = 0;

        % --- Newton iteration loop ------------------------------------
        while itercnt < MAXITER && f > eps_iter && ...
              (real(z) + real(dz) ~= real(z) || imag(z) + imag(dz) ~= imag(z))

            % Evaluate p'(z) using paired-real Horner on da (degree N-1)
            [u_deriv, fz1] = peval(N - 1, da, z);

            if u_deriv == 0
                % True saddle point: rotate and enlarge step to escape
                dz = rotate_step(dz, 5.0);
            else
                % Newton step dz = p(z)/p'(z) via conjugate division
                dz = complex( ...
                    (real(fz)*real(fz1) + imag(fz)*imag(fz1)) / u_deriv, ...
                    (imag(fz)*real(fz1) - real(fz)*imag(fz1)) / u_deriv);

                % Kantorowitz-type stage decision ----------------------
                fwz   = fz0 - fz1;          % Δp' ≈ change in derivative
                wdiff = z0  - z;            % Δz
                if abs(wdiff) > 0
                    f2     = (real(fwz)^2 + imag(fwz)^2) / ...
                             (real(wdiff)^2 + imag(wdiff)^2);
                    stage1 = (f2 / u_deriv > u_deriv / (4 * f)) || (f ~= ff);
                else
                    stage1 = true;
                end

                r = abs(dz);
                if r > r0
                    dz = rotate_step(dz, r0 / r);   % cap step magnitude
                end
                r0 = abs(dz) * 5.0;
            end

            % Save current iterate for the next stage decision
            z0 = z;  f0 = f;  fz0 = fz;

            % --- Step / line-search / rounding-guard inner loop --------
            % This mirrors the C++ "goto iter2" structure: when the step
            % collapses below machine resolution without improving |p|,
            % the direction is altered and the step retried.
            retry = true;
            while retry
                retry = false;
                z       = z0 - dz;
                [f, fz] = peval(N, a, z);
                ff      = f;

                if stage1
                    % Line search: halve or extend dz to guarantee |p|
                    % decreases. Halve when initial step overshoots
                    % (f > f0); extend when it undershoots.
                    wz   = z;
                    div2 = (f > f0);
                    for i = 1:N
                        if div2
                            dz = dz * 0.5;
                            wz = z0 - dz;
                        else
                            wz = wz - dz;        % extend one Newton length
                        end
                        [fw, fwz_val] = peval(N, a, wz);
                        if fw >= f, break; end
                        f = fw;  fz = fwz_val;  z = wz;
                        if div2 && i == 2
                            % Two halvings insufficient: rotate and retry
                            dz = rotate_step(dz, 0.5);
                            z  = z0 - dz;
                            [f, fz] = peval(N, a, z);
                            break;
                        end
                    end
                else
                    % Stage 2: tighten convergence threshold to Adam's
                    % local rounding-error bound.
                    eps_iter = adam_bound(N, a, z);
                end

                rstep = abs(dz);
                if rstep < abs(z) * SQRT_EPS && f >= f0
                    % Domain rounding-error guard: alter direction and
                    % retry; if the iterate is still unchanged, accept.
                    z  = z0;
                    dz = rotate_step(dz, 0.5);
                    if real(z) + real(dz) ~= real(z) || ...
                       imag(z) + imag(dz) ~= imag(z)
                        retry = true;
                    end
                end
            end

            itercnt = itercnt + 1;
        end % Newton iteration loop

        if itercnt >= MAXITER
            err = err - 1;   % flag that this root did not converge
        end

        % --- Classify root: real or complex conjugate pair ---------------
        % Exploit real-coefficient symmetry: if projecting z to the real
        % axis does not worsen the residual, the root is real.
        [f_re, ~] = peval(N, a, complex(real(z), 0));
        if f_re <= f
            roots_out(rp) = real(z);
            a  = real_deflate(N, a, real(z));
            N  = N - 1;
            rp = rp - 1;
        else
            roots_out(rp)     = z;
            roots_out(rp - 1) = conj(z);
            a  = complex_deflate(N, a, z);
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
%   Standard lower bound for the modulus of the smallest root (Cauchy,
%   1829). The Newton starting iterate is placed at u/2.

    an = abs(a(n + 1));
    if an == 0, u = 1; return; end

    rLog = log(an);
    u    = exp((rLog - log(abs(a(1)))) / n);
    for i = 2:n
        if a(i) ~= 0
            cand = exp((rLog - log(abs(a(i)))) / (n - i + 1));
            if cand < u
                u = cand;
            end
        end
    end
end


function [mod2, val] = peval(n, a, z)
%PEVAL  Evaluate degree-n real-coefficient polynomial at complex z.
%   Returns (|p(z)|², p(z)) using the paired-real Horner recurrence:
%
%       p = -2·Re(z),  q = |z|²
%       r_k = a_k − p·r_{k-1} − q·r_{k-2}
%
%   The recurrence operates entirely in real arithmetic until the final
%   reconstruction p(z) = (a_n + Re(z)·r − q·s) + i·Im(z)·r, halving the
%   multiplications of a generic complex Horner.

    p = -2.0 * real(z);
    q = real(z)^2 + imag(z)^2;
    s = 0;
    r = a(1);
    for i = 2:n
        t = a(i) - p*r - q*s;
        s = r;
        r = t;
    end
    val  = complex(a(n+1) + real(z)*r - q*s, imag(z)*r);
    mod2 = real(val)^2 + imag(val)^2;
end


function dz = rotate_step(dz, m)
%ROTATE_STEP  Rotate step vector by ~53.13° and scale by m.
%
%   Applies a fixed Givens-like rotation:
%       [ 0.6  -0.8 ] [ Re(dz) ]
%       [ 0.8   0.6 ] [ Im(dz) ] * m
%
%   Rotating the Newton step escapes saddle points and avoids oscillation
%   on quadratic-factor axes. m < 1 reduces step magnitude; m > 1 enlarges.

    re = real(dz);  im = imag(dz);
    dz = complex(re * 0.6 - im * 0.8, ...
                 re * 0.8 + im * 0.6) * m;
end


function e2 = adam_bound(n, a, z)
%ADAM_BOUND  Adam's rounding-error upper bound on |p(z)|², squared.
%
%   Returns e² so it can be compared directly against |p(z)|² = f.
%   Mirrors C++ upperbound() in Newton_Real.cpp.

    p = -2.0 * real(z);
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
    e = (9*e - 7*(abs(t) + abs(r)*u) + abs(real(z))*abs(r)*2) * eps(1);
    e2 = e^2;
end


function a = real_deflate(n, a, x)
%REAL_DEFLATE  Synthetic division of a (degree n) by (z − x).
%   Result is degree n−1 stored in-place in a(1:n).

    r = 0;
    for i = 1:n
        r = r * x + a(i);
        a(i) = r;
    end
end


function a = complex_deflate(n, a, z)
%COMPLEX_DEFLATE  Synthetic division of a (degree n) by
%   (x² − 2·Re(z)·x + |z|²). Result is degree n−2 stored in a(1:n−1).

    rr = -2.0 * real(z);
    u  = real(z)^2 + imag(z)^2;
    a(2) = a(2) - rr * a(1);
    for i = 3:n-1
        a(i) = a(i) - rr*a(i-1) - u*a(i-2);
    end
end


function res = solve_low_degree(a, n)
%SOLVE_LOW_DEGREE  Exact roots of a linear (n=1) or quadratic (n=2).
%
%   a is in leading-first order: a(1)·x^n + a(2)·x^(n-1) + … + a(n+1).

    res = complex(zeros(n, 1));

    if n == 1
        res(1) = -a(2) / a(1);

    else  % n == 2 quadratic: a(1)·x² + a(2)·x + a(3) = 0
        if a(2) == 0
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
            disc = 1 - 4 * a(1) * a(3) / a(2)^2;
            if disc < 0
                re_part = -a(2) / (2 * a(1));
                im_part =  a(2) * sqrt(-disc) / (2 * a(1));
                res(1) = complex(re_part,  im_part);
                res(2) = complex(re_part, -im_part);
            else
                % Two real roots — use Vieta's product to avoid cancellation
                res(1) = (-1 - sqrt(disc)) * a(2) / (2 * a(1));
                res(2) = a(3) / (a(1) * real(res(1)));
            end
        end
    end
end
