 using IntervalArithmetic





setprecision(BigFloat, 256)





function pow_log_interval(a, exponent)
    if inf(a) > 0
        return (a^exponent) * log(a)
    end
    if inf(exponent) <= 0
        return entireinterval(BigFloat)
    end
    if sup(a) == 0
        return interval(0)
    end
    ahi = sup(a)
    elo = inf(exponent)
    if !(ahi > 0 && elo > 0 && elo * log(ahi) + 1 <= 0)
        return entireinterval(BigFloat)
    end
    lower = (interval(ahi)^interval(elo)) * log(interval(ahi))
    return interval(inf(lower), 0)
end

C_interval(p, x) = (interval(1) + x^p)^(-interval(1) / p)

function C_scalar(p, x)
    pp = big(p)
    xx = big(x)
    return (big(1) + xx^pp)^(-big(1) / pp)
end

function F_constraint_scalar(p, s, t)
    pp = big(p)
    tt = big(t)
    ss = big(s)
    ct = C_scalar(pp, tt)
    cs = C_scalar(pp, ss)
    A = ct - cs
    B = tt * ct + ss * cs
    return A^pp + B^pp - 1
end

function F_constraint_interval(p, s, t)
    ct = C_interval(p, t)
    cs = C_interval(p, s)
    A = ct - cs
    B = t * ct + s * cs
    return A^p + B^p - interval(1)
end

function L_C_p(x, p)
    xp = x^p
    xp_logx = pow_log_interval(x, p)
    return log(interval(1) + xp) / (p^interval(2)) - xp_logx / (p * (interval(1) + xp))
end

function F_partials_interval(p, s, t)
    p_minus_1 = p - interval(1)
    one_over_p = interval(1) / p

    sp = s^p
    tp = t^p
    C1 = (interval(1) + sp)^(-one_over_p)
    C2 = (interval(1) + tp)^(-one_over_p)
    Ksigma = (interval(1) + sp)^(-interval(1) - one_over_p)
    Ktau = (interval(1) + tp)^(-interval(1) - one_over_p)

    A = C2 - C1
    B = t * C2 + s * C1

    A_pm1 = A^p_minus_1
    B_pm1 = B^p_minus_1
    s_pm1 = s^p_minus_1
    t_pm1 = t^p_minus_1

    H1 = B_pm1 + s_pm1 * A_pm1
    H2 = B_pm1 - t_pm1 * A_pm1

    Lsigma = L_C_p(s, p)
    Ltau = L_C_p(t, p)
    C1_p = C1 * Lsigma
    C2_p = C2 * Ltau
    A_p = C2_p - C1_p
    B_p = t * C2_p + s * C1_p

    Fp = A^p * log(A) + p * A_pm1 * A_p +
         B^p * log(B) + p * B_pm1 * B_p
    Fs = p * Ksigma * H1
    Ft = p * Ktau * H2
    return Fp, Fs, Ft
end

"""
    merge_intervals(intervals)

Merges a list of sorted, contiguous, or overlapping intervals into a minimal
list of disjoint intervals.
"""
function merge_intervals(intervals)
    if isempty(intervals)
        return []
    end

    sorted_intervals = sort(intervals, by = x -> inf(x))
    merged = [sorted_intervals[1]]

    for i in 2:length(sorted_intervals)
        current_interval = sorted_intervals[i]
        last_merged = merged[end]

        if sup(last_merged) >= inf(current_interval)
            new_hi = max(sup(last_merged), sup(current_interval))
            merged[end] = interval(inf(last_merged), new_hi)
        else
            push!(merged, current_interval)
        end
    end

    return merged
end

"""
    subdivide_tau_interval(tau_rect, num_subdivisions)

Splits a tau interval into a list of N subintervals.
"""
function subdivide_tau_interval(tau_rect::Interval{<:Real}, num_subdivisions::Int)
    sub_intervals = typeof(tau_rect)[]

    t_lo = inf(tau_rect)
    t_hi = sup(tau_rect)
    t_step = (t_hi - t_lo) / num_subdivisions

    tau_values_asc = collect(t_lo:t_step:t_hi)
    if last(tau_values_asc) < t_hi
        push!(tau_values_asc, t_hi)
    end

    for k in 1:(length(tau_values_asc) - 1)
        push!(sub_intervals, interval(tau_values_asc[k], tau_values_asc[k+1]))
    end

    return sub_intervals
end

sigma_p_interval(p) = (interval(2)^p - interval(1))^(interval(1) / p)





function sigma_scalar_from_tau_p2(p, t; lo=sqrt(big(3)), hi)
    root_zero_tol = big(1) / big(10)^90
    flo = F_constraint_scalar(p, lo, t)
    fhi = F_constraint_scalar(p, hi, t)
    abs(flo) < root_zero_tol && return lo
    abs(fhi) < root_zero_tol && return hi
    sign(flo) == sign(fhi) && error("sigma root not bracketed: p=$p t=$t flo=$flo fhi=$fhi lo=$lo hi=$hi")
    for _ in 1:180
        mid = (lo + hi) / 2
        fm = F_constraint_scalar(p, mid, t)
        if sign(flo) == sign(fm)
            lo = mid
            flo = fm
        else
            hi = mid
        end
    end
    return (lo + hi) / 2
end

function sigma_corner_bracket_p2(p, t; lo=sqrt(big(3)))
    tau_zero_tol = big(1) / big(10)^300
    root_radius_abs = big(1) / big(10)^80
    root_radius_rel = big(1) / big(10)^60

    hi_int = sigma_p_interval(interval(p))
    hi_scalar = sup(hi_int)
    abs(t) < tau_zero_tol && return hi_int

    root = sigma_scalar_from_tau_p2(p, t; lo=lo, hi=hi_scalar)
    radius = max(root_radius_abs, abs(root) * root_radius_rel)
    for _ in 1:80
        a = max(lo, root - radius)
        b = min(hi_scalar, root + radius)
        Fa = F_constraint_interval(interval(p), interval(a), interval(t))
        Fb = F_constraint_interval(interval(p), interval(b), interval(t))
        if sup(Fa) <= 0 && inf(Fb) >= 0
            return interval(a, b)
        end
        radius *= 2
    end
    error("could not certify sigma corner bracket: p=$p t=$t root=$root")
end

function sigma_box_from_ptau(pbox, tbox; lo=sqrt(big(3)))
    brackets = []
    for pp in (inf(pbox), sup(pbox)), tt in (inf(tbox), sup(tbox))
        push!(brackets, sigma_corner_bracket_p2(pp, tt; lo=lo))
    end
    merged_brackets = merge_intervals(brackets)
    return interval(inf(merged_brackets[1]), sup(merged_brackets[end]))
end

function sigma_box_certificate(pbox, tbox, sbox)
    Fp, Fs, Ft = F_partials_interval(pbox, sbox, tbox)
    return inf(Fs) > 0 && sup(Fp) < 0 && inf(Ft) > 0, Fp, Fs, Ft
end





function calculate_h_desingularized_interval(p, s, t, x)
    one = interval(1)
    two = interval(2)
    delta = p - two
    one_over_p = one / p
    p_minus_1 = p - one

    s_p = s^p
    t_p = (t^two) * x

    val_1_sp = one + s_p
    val_1_tp = one + t_p

    C1 = val_1_sp^(-one_over_p)
    C2 = val_1_tp^(-one_over_p)

    K_sigma = val_1_sp^(-one - one_over_p)
    K_tau = val_1_tp^(-one - one_over_p)

    A = C2 - C1
    B = t * C2 + s * C1

    s_pm1 = s^p_minus_1
    s_pm2 = s^delta
    t_pm1 = t * x
    t_pm2 = x

    A_pm1 = A^p_minus_1
    A_pm2 = A^delta
    B_pm1 = B^p_minus_1
    B_pm2 = B^delta

    D1 = one - s * t_pm1
    D2 = one - t * s_pm1

    H1 = B_pm1 + s_pm1 * A_pm1
    H2 = B_pm1 - t_pm1 * A_pm1

    S = K_sigma / K_tau
    R = H1 / H2
    T = -S * R

    A_s = s_pm1 * K_sigma
    A_t = -t_pm1 * K_tau
    B_s = K_sigma
    B_t = K_tau

    C1_s = -s_pm1 * K_sigma
    C2_t = -t_pm1 * K_tau

    D1_s = -t_pm1
    D2_s = -t * p_minus_1 * s_pm2
    D1_t = -s * p_minus_1 * t_pm2
    D2_t = -s_pm1

    H1_s = p_minus_1 * (B_pm2 * B_s + s_pm2 * A_pm1 + s_pm1 * A_pm2 * A_s)
    H2_s = p_minus_1 * (B_pm2 * B_s - t_pm1 * A_pm2 * A_s)
    H1_t = p_minus_1 * (B_pm2 * B_t + s_pm1 * A_pm2 * A_t)
    H2_t = p_minus_1 * (B_pm2 * B_t - t_pm2 * A_pm1 - t_pm1 * A_pm2 * A_t)

    g_s = C1_s * D1 * H1 + C1 * D1_s * H1 + C1 * D1 * H1_s -
          C2 * D2_s * H2 - C2 * D2 * H2_s

    g_t = C1 * D1_t * H1 + C1 * D1 * H1_t -
          C2_t * D2 * H2 - C2 * D2_t * H2 - C2 * D2 * H2_t

    return g_s + g_t * T
end

function calculate_h_partials_interval(p, s, t, x)
    with_partials(v, dp, ds, dx) = (v, dp, ds, dx)
    constant(v) = with_partials(v, interval(0), interval(0), interval(0))
    p_variable(v) = with_partials(v, interval(1), interval(0), interval(0))
    sigma_variable(v) = with_partials(v, interval(0), interval(1), interval(0))
    x_variable(v) = with_partials(v, interval(0), interval(0), interval(1))

    negate(a) = with_partials(-a[1], -a[2], -a[3], -a[4])
    add(a, b) = with_partials(a[1] + b[1], a[2] + b[2], a[3] + b[3], a[4] + b[4])
    subtract(a, b) = with_partials(a[1] - b[1], a[2] - b[2], a[3] - b[3], a[4] - b[4])

    function multiply(a, b)
        return with_partials(
            a[1] * b[1],
            a[2] * b[1] + a[1] * b[2],
            a[3] * b[1] + a[1] * b[3],
            a[4] * b[1] + a[1] * b[4],
        )
    end

    function reciprocal(a)
        return with_partials(
            inv(a[1]),
            -a[2] / (a[1]^interval(2)),
            -a[3] / (a[1]^interval(2)),
            -a[4] / (a[1]^interval(2)),
        )
    end

    divide(a, b) = multiply(a, reciprocal(b))

    function power_constant(a, exponent)
        y = a[1]^exponent
        if sup(abs(a[2])) == 0 && sup(abs(a[3])) == 0 && sup(abs(a[4])) == 0
            return constant(y)
        end
        factor = exponent * (a[1]^(exponent - interval(1)))
        return with_partials(y, factor * a[2], factor * a[3], factor * a[4])
    end

    function power(a, b)
        y = a[1]^b[1]
        loga = log(a[1])
        return with_partials(
            y,
            y * (b[2] * loga + b[1] * a[2] / a[1]),
            y * (b[3] * loga + b[1] * a[3] / a[1]),
            y * (b[4] * loga + b[1] * a[4] / a[1]),
        )
    end

    one = constant(interval(1))
    two = constant(interval(2))
    pp = p_variable(p)
    ss = sigma_variable(s)
    tt = constant(t)
    xx = x_variable(x)

    delta = subtract(pp, two)
    one_over_p = divide(one, pp)
    p_minus_1 = subtract(pp, one)

    s_p = power(ss, pp)
    t_p = multiply(power_constant(tt, interval(2)), xx)

    val_1_sp = add(one, s_p)
    val_1_tp = add(one, t_p)

    C1 = power(val_1_sp, negate(one_over_p))
    C2 = power(val_1_tp, negate(one_over_p))

    K_sigma = power(val_1_sp, subtract(negate(one), one_over_p))
    K_tau = power(val_1_tp, subtract(negate(one), one_over_p))

    A = subtract(C2, C1)
    B = add(multiply(tt, C2), multiply(ss, C1))

    s_pm1 = power(ss, p_minus_1)
    s_pm2 = power(ss, delta)
    t_pm1 = multiply(tt, xx)
    t_pm2 = xx

    A_pm1 = power(A, p_minus_1)
    A_pm2 = power(A, delta)
    B_pm1 = power(B, p_minus_1)
    B_pm2 = power(B, delta)

    D1 = subtract(one, multiply(ss, t_pm1))
    D2 = subtract(one, multiply(tt, s_pm1))

    H1 = add(B_pm1, multiply(s_pm1, A_pm1))
    H2 = subtract(B_pm1, multiply(t_pm1, A_pm1))

    S = divide(K_sigma, K_tau)
    R = divide(H1, H2)
    T = negate(multiply(S, R))

    A_s = multiply(s_pm1, K_sigma)
    A_t = negate(multiply(t_pm1, K_tau))
    B_s = K_sigma
    B_t = K_tau

    C1_s = negate(multiply(s_pm1, K_sigma))
    C2_t = negate(multiply(t_pm1, K_tau))

    D1_s = negate(t_pm1)
    D2_s = negate(multiply(multiply(tt, p_minus_1), s_pm2))
    D1_t = negate(multiply(multiply(ss, p_minus_1), t_pm2))
    D2_t = negate(s_pm1)

    H1_s = multiply(p_minus_1,
        add(add(multiply(B_pm2, B_s), multiply(s_pm2, A_pm1)),
             multiply(multiply(s_pm1, A_pm2), A_s)))
    H2_s = multiply(p_minus_1,
        subtract(multiply(B_pm2, B_s), multiply(multiply(t_pm1, A_pm2), A_s)))
    H1_t = multiply(p_minus_1,
        add(multiply(B_pm2, B_t), multiply(multiply(s_pm1, A_pm2), A_t)))
    H2_t = multiply(p_minus_1,
        subtract(subtract(multiply(B_pm2, B_t), multiply(t_pm2, A_pm1)),
             multiply(multiply(t_pm1, A_pm2), A_t)))

    g_s = subtract(
        add(add(multiply(multiply(C1_s, D1), H1), multiply(multiply(C1, D1_s), H1)),
             multiply(multiply(C1, D1), H1_s)),
        add(multiply(multiply(C2, D2_s), H2), multiply(multiply(C2, D2), H2_s)),
    )

    g_t = subtract(
        add(multiply(multiply(C1, D1_t), H1), multiply(multiply(C1, D1), H1_t)),
        add(add(multiply(multiply(C2_t, D2), H2), multiply(multiply(C2, D2_t), H2)),
             multiply(multiply(C2, D2), H2_t)),
    )

    h = add(g_s, multiply(g_t, T))
    return h[1], h[2], h[3], h[4]
end





function fixedtau_x_box(delta, tau)
    sup(delta) == 0 && return big(interval(1//1))
    inf(tau) <= 0 && return big(interval(0//1, 1//1))
    zlo = sup(delta) * log(inf(tau))
    zhi = inf(delta) * log(sup(tau))
    return exp(interval(zlo, zhi))
end

function fixedtau_x_delta_box(tau, x_lower=0)
    upper = x_lower <= 0 ? BigFloat(0) : x_lower * log(sup(tau))
    return interval(BigFloat(-Inf), upper)
end

function calculate_Hdelta_interval(delta, tau, x, sigma_margin)
    p = big(interval(2//1)) + delta
    s = sigma_box_from_ptau(p, tau; lo=sqrt(big(3)) - sigma_margin)
    branch_cert, Fp, Fs, Ft = sigma_box_certificate(p, tau, s)

    sigma_delta = -Fp / Fs
    x_delta = fixedtau_x_delta_box(tau, inf(x))
    _, hp, hs, hx = calculate_h_partials_interval(p, s, tau, x)
    H_delta = hp + hs * sigma_delta + hx * x_delta
    return H_delta, branch_cert, s, x_delta, sigma_delta, hp, hs, hx, Fp, Fs, Ft
end





function check_tiny_tau_collar(delta_max, tau_max, sigma_margin)
    delta_box = delta_max
    p_box = big(interval(2//1)) + delta_box
    tau_box = interval(big(0), sup(tau_max))
    sigma_box = sigma_box_from_ptau(p_box, tau_box; lo=sqrt(big(3)) - sigma_margin)

    x_low = big(interval(0//1, 1//2))
    h_low = calculate_h_desingularized_interval(p_box, sigma_box, tau_box, x_low)

    x_high = big(interval(1//2, 1//1))
    H_delta, branch_cert, sigma_high, x_delta, sigma_delta,
        hp, hs, hx, Fp, Fs, Ft =
        calculate_Hdelta_interval(delta_box, tau_box, x_high, sigma_margin)

    ok = sup(h_low) < 0 && sup(H_delta) < 0 && branch_cert &&
         inf(Fs) > 0 && inf(Ft) > 0 && sup(Fp) < 0

    return ok, h_low, H_delta, branch_cert, sigma_box, sigma_high,
        x_delta, sigma_delta, hp, hs, hx, Fp, Fs, Ft
end

function check_tau_boxes(;
        delta_max,
        tau_start,
        tau_end,
        tau_subdivisions,
        sigma_margin,
        verbose_every=10)

    failed_tau_intervals = Interval{BigFloat}[]
    checked_boxes = 0
    worst_upper_bound = -Inf

    tau_rect = interval(inf(tau_start), sup(tau_end))

    for tau_box in subdivide_tau_interval(tau_rect, tau_subdivisions)
        x_box = fixedtau_x_box(delta_max, tau_box)
        H_delta, branch_cert, sigma_box, x_delta, sigma_delta,
            hp, hs, hx, Fp, Fs, Ft =
            calculate_Hdelta_interval(delta_max, tau_box, x_box, sigma_margin)

        checked_boxes += 1
        worst_upper_bound = max(worst_upper_bound, sup(H_delta))
        println("  Tau=$tau_box X=$x_box S=$sigma_box => H_delta=$H_delta")

        if !(branch_cert && inf(Fs) > 0 && inf(Ft) > 0 && sup(Fp) < 0 && sup(H_delta) < 0)
            println("      !!! FAIL: tau box did not certify.")
            println("          branch_cert=$branch_cert Fp=$Fp Fs=$Fs Ft=$Ft")
            println("          x_delta=$x_delta sigma_delta=$sigma_delta")
            println("          hp=$hp hs=$hs hx=$hx")
            push!(failed_tau_intervals, tau_box)
        elseif verbose_every > 0 && checked_boxes % verbose_every == 0
            println("      certified boxes so far: $checked_boxes, worst upper bound: $worst_upper_bound")
        end
    end

    return isempty(failed_tau_intervals), failed_tau_intervals, checked_boxes, worst_upper_bound
end

function check_tau_coverage(delta_max, sigma_min, tau_end)
    p_box = big(interval(2//1)) + delta_max
    sigma_endpoint_hi = sup(sigma_p_interval(p_box))
    sigma_range = interval(inf(sigma_min), sigma_endpoint_hi)
    tau_range = interval(big(0), sup(tau_end))

    F_at_top = F_constraint_interval(p_box, sigma_min, tau_end)
    _, Fs, Ft = F_partials_interval(p_box, sigma_range, tau_range)
    ok = inf(F_at_top) > 0 && inf(Fs) > 0 && inf(Ft) > 0
    return ok, F_at_top, Fs, Ft, sigma_range, tau_range
end





# Certificate parameters.  These are exact rational interval endpoints.
delta_box = big(interval(0//1, 1//5000))  # 0 <= p-2 <= 2e-4
sigma_min = big(interval(43//25))         # sigma >= 1.72
tau_start = big(interval(1//(10^12)))     # positive boxes start
tau_end = big(interval(1//250))           # tau <= 0.004
tiny_tau_end = big(interval(1//(10^12)))  # tiny endpoint collar
tau_subdivisions = 40
sigma_margin = big(1) / big(50)

println("Calculating the p=2+ direct-tau interval certificate.")
println("  Julia version: ", VERSION)
println("  BigFloat precision: ", precision(BigFloat), " bits")
println("  p range:       [2, 2 + 1/5000]")
println("  sigma range:   [$sigma_min, sigma_p]")
println("  tau coverage:  [0, $tau_end]")
println("  tau boxes:     [$tau_start, $tau_end] split into $tau_subdivisions uniform pieces")
println("  Inequalities:  branch regularity, h<0 / H_delta<0, tau coverage")
println("="^40)

start_time = time()
proof_success = true

println("-"^150)
println("PART A: Tiny tau collar, 0 <= tau <= $tiny_tau_end")
println("-"^150)
tiny_ok, h_low, H_delta_high, tiny_branch_cert, sigma_low, sigma_high,
    x_delta, sigma_delta, hp_tiny, hs_tiny, hx_tiny, Fp_tiny, Fs_tiny, Ft_tiny =
    check_tiny_tau_collar(delta_box, tiny_tau_end, sigma_margin)
println("  sigma enclosure for x<=1/2: $sigma_low")
println("  h(x<=1/2)=$h_low")
println("  sigma enclosure for x>=1/2: $sigma_high")
println("  H_delta(x>=1/2)=$H_delta_high")
println("  branch_cert=$tiny_branch_cert Fp=$Fp_tiny Fs=$Fs_tiny Ft=$Ft_tiny")
println("  x_delta=$x_delta sigma_delta=$sigma_delta")
println("  hp=$hp_tiny hs=$hs_tiny hx=$hx_tiny")
println("  tiny tau collar ok: $tiny_ok")
proof_success &= tiny_ok

println("-"^150)
println("PART B: Positive tau boxes")
println("-"^150)
println("  Tau-interval        X-Enclosure       Sigma-Enclosure       H_delta Bound")
println("-"^150)
boxes_ok, failed_tau_intervals, checked_boxes, worst_H_delta =
    check_tau_boxes(
        delta_max=delta_box,
        tau_start=tau_start,
        tau_end=tau_end,
        tau_subdivisions=tau_subdivisions,
        sigma_margin=sigma_margin,
    )
proof_success &= boxes_ok

println("-"^150)
println("PART C: Tau coverage of the whole sigma strip")
println("-"^150)
coverage_ok, F_at_top, Fs_cover, Ft_cover, sigma_cover, tau_cover =
    check_tau_coverage(delta_box, sigma_min, tau_end)
println("  sigma range=$sigma_cover")
println("  tau range=$tau_cover")
println("  F(P, sigma_min, tau_end)=$F_at_top")
println("  F_sigma=$Fs_cover")
println("  F_tau=$Ft_cover")
println("  coverage ok: $coverage_ok")
proof_success &= coverage_ok

end_time = time()
println("="^40)
println("Final Report on Inconclusive Tau-Intervals:")
if proof_success && isempty(failed_tau_intervals)
    println("SUCCESS: the p=2+ direct-tau neighbourhood was verified.")
else
    println("FAILURE: at least one certificate block was inconclusive.")
    if !isempty(failed_tau_intervals)
        println("The following tau-ranges failed the H_delta check:")
        merged_failed = merge_intervals(failed_tau_intervals)
        for tau_int in merged_failed
            println("  ", tau_int)
        end
    end
end
println("Checked positive tau boxes: ", checked_boxes)
println("Worst H_delta upper bound: ", worst_H_delta)
println("Total time: ", round(end_time - start_time, digits=4), " seconds.")

exit(proof_success ? 0 : 1)
