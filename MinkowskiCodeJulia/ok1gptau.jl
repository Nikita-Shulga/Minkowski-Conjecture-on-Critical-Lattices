using IntervalArithmetic



"""
    find_tp_range(p_range, initial_tp_guess, precision, verbose)

Iteratively finds an interval for t_p from
    2*(1-t_p)^p = 1 + t_p^p.
"""
function find_tp_range(
    p_range::Interval{<:Real},
    initial_tp_guess::Interval{<:Real},
    precision::Real,
    verbose::Bool = false
)
    tp_interval = initial_tp_guess
    one_over_p = interval(1) / p_range
    two_pow_one_over_p = interval(2)^one_over_p

    for i in 1:100
        if verbose
            println("  t_p Iteration $i: ", tp_interval)
        end

        if isempty_interval(tp_interval)
            println("  Warning: t_p interval is empty.")
            return emptyinterval(Float64)
        end

        # In the p=1 neighbourhood t_p is close to 1/3.
        tp_interval = intersect_interval(tp_interval, interval(1//10, 36//100); dec=:auto)

        tp_pow_p = tp_interval^p_range
        one_plus_tp_pow_p = interval(1) + tp_pow_p
        term_pow = one_plus_tp_pow_p^one_over_p
        full_term = term_pow / two_pow_one_over_p
        new_tp_interval = interval(1) - full_term

        if diam(new_tp_interval) < precision
            if verbose
                println("  Converged t_p to $new_tp_interval")
            end
            return intersect_interval(tp_interval, new_tp_interval; dec=:auto)
        end

        tp_interval = intersect_interval(tp_interval, new_tp_interval; dec=:auto)
    end

    if verbose
        println("  t_p max iterations reached. Final interval: $tp_interval")
    end
    return tp_interval
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

        if sup(last_merged) >= (inf(current_interval) - 1e-12)
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



"""
    calculate_branch_gptau_interval(p_int, s_int, t_int)

Calculates the mixed derivative
    d^2/dp dtau g(p, sigma(p,tau), tau)
using only explicit analytic interval formulae.


All derivatives below are raw partial derivatives in the independent variables
p, sigma, tau. They are then combined using

    sigma_p   = -F_p/F_sigma,
    sigma_tau = -F_tau/F_sigma,
    sigma_pt  = -(F_pt + F_ps*sigma_tau + F_st*sigma_p
                  + F_ss*sigma_p*sigma_tau)/F_sigma,

and

    branch_g_pt = g_pt + g_ps*sigma_tau + g_st*sigma_p
                  + g_ss*sigma_p*sigma_tau + g_s*sigma_pt.
"""
function calculate_branch_gptau_interval(
    p_int::Interval{<:Real},
    s_int::Interval{<:Real},
    t_int::Interval{<:Real}
)
    # --------------------------------------------------------------------------
    # 1. Basic powers, logarithms and C/K factors
    # --------------------------------------------------------------------------
    p_minus_1 = p_int - interval(1)
    p_minus_2 = p_int - interval(2)
    p_minus_3 = p_int - interval(3)

    ln_s = log(s_int)
    ln_t = log(t_int)

    s_pow_p = s_int^p_int
    t_pow_p = t_int^p_int

    one_plus_sp = interval(1) + s_pow_p
    one_plus_tp = interval(1) + t_pow_p

    C1 = one_plus_sp^(-interval(1) / p_int)
    C2 = one_plus_tp^(-interval(1) / p_int)

    Ksigma = one_plus_sp^(-interval(1) - interval(1) / p_int)
    Ktau   = one_plus_tp^(-interval(1) - interval(1) / p_int)

    Lsigma = log(one_plus_sp) / (p_int^interval(2)) - (s_pow_p * ln_s) / (p_int * one_plus_sp)
    Ltau   = log(one_plus_tp) / (p_int^interval(2)) - (t_pow_p * ln_t) / (p_int * one_plus_tp)

    Lsigma_s = -s_int^p_minus_1 * ln_s / (one_plus_sp^interval(2))
    Ltau_t   = -t_int^p_minus_1 * ln_t / (one_plus_tp^interval(2))

    # C1 = C_p(sigma), C2 = C_p(tau)
    C1_p  = C1 * Lsigma
    C1_s  = -s_int^p_minus_1 * Ksigma
    C1_t  = interval(0)
    C1_ps = C1_s * Lsigma + C1 * Lsigma_s
    C1_pt = interval(0)
    C1_st = interval(0)
    C1_ss = -p_minus_1 * s_int^p_minus_2 * Ksigma +
            (p_int + interval(1)) * s_int^(interval(2)*p_int - interval(2)) *
            one_plus_sp^(-interval(2) - interval(1) / p_int)

    C2_p  = C2 * Ltau
    C2_s  = interval(0)
    C2_t  = -t_int^p_minus_1 * Ktau
    C2_ps = interval(0)
    C2_pt = C2_t * Ltau + C2 * Ltau_t
    C2_st = interval(0)
    C2_ss = interval(0)

    # --------------------------------------------------------------------------
    # 2. A and B and their raw partial derivatives
    # --------------------------------------------------------------------------
    A = C2 - C1
    A_p  = C2_p - C1_p
    A_s  = C2_s - C1_s
    A_t  = C2_t - C1_t
    A_ps = C2_ps - C1_ps
    A_pt = C2_pt - C1_pt
    A_st = C2_st - C1_st
    A_ss = C2_ss - C1_ss

    B = t_int * C2 + s_int * C1
    B_p  = t_int * C2_p + s_int * C1_p
    B_s  = C1 + s_int * C1_s
    B_t  = C2 + t_int * C2_t
    B_ps = C1_p + s_int * C1_ps
    B_pt = C2_p + t_int * C2_pt
    B_st = interval(0)
    B_ss = interval(2) * C1_s + s_int * C1_ss

    # --------------------------------------------------------------------------
    # 3. Powers A^p and B^p for the constraint F=A^p+B^p-1
    # --------------------------------------------------------------------------
    Apow = A^p_int
    Apow_p  = Apow * log(A) + p_int * A^p_minus_1 * A_p
    Apow_s  = p_int * A^p_minus_1 * A_s
    Apow_t  = p_int * A^p_minus_1 * A_t
    Apow_ps = A^p_minus_1 * A_s * (interval(1) + p_int * log(A)) +
              p_int * p_minus_1 * A^p_minus_2 * A_p * A_s +
              p_int * A^p_minus_1 * A_ps
    Apow_pt = A^p_minus_1 * A_t * (interval(1) + p_int * log(A)) +
              p_int * p_minus_1 * A^p_minus_2 * A_p * A_t +
              p_int * A^p_minus_1 * A_pt
    Apow_st = p_int * (p_minus_1 * A^p_minus_2 * A_s * A_t + A^p_minus_1 * A_st)
    Apow_ss = p_int * (p_minus_1 * A^p_minus_2 * A_s^interval(2) + A^p_minus_1 * A_ss)

    Bpowp = B^p_int
    Bpowp_p  = Bpowp * log(B) + p_int * B^p_minus_1 * B_p
    Bpowp_s  = p_int * B^p_minus_1 * B_s
    Bpowp_t  = p_int * B^p_minus_1 * B_t
    Bpowp_ps = B^p_minus_1 * B_s * (interval(1) + p_int * log(B)) +
               p_int * p_minus_1 * B^p_minus_2 * B_p * B_s +
               p_int * B^p_minus_1 * B_ps
    Bpowp_pt = B^p_minus_1 * B_t * (interval(1) + p_int * log(B)) +
               p_int * p_minus_1 * B^p_minus_2 * B_p * B_t +
               p_int * B^p_minus_1 * B_pt
    Bpowp_st = p_int * (p_minus_1 * B^p_minus_2 * B_s * B_t + B^p_minus_1 * B_st)
    Bpowp_ss = p_int * (p_minus_1 * B^p_minus_2 * B_s^interval(2) + B^p_minus_1 * B_ss)

    F_p  = Apow_p  + Bpowp_p
    F_s  = Apow_s  + Bpowp_s
    F_t  = Apow_t  + Bpowp_t
    F_ps = Apow_ps + Bpowp_ps
    F_pt = Apow_pt + Bpowp_pt
    F_st = Apow_st + Bpowp_st
    F_ss = Apow_ss + Bpowp_ss

    if in_interval(0, F_s) == true
        println("  Warning: F_sigma interval contains zero.")
        return emptyinterval(Float64)
    end

    sigma_p   = -F_p / F_s
    sigma_tau = -F_t / F_s
    sigma_pt  = -(F_pt + F_ps * sigma_tau + F_st * sigma_p + F_ss * sigma_p * sigma_tau) / F_s

    # --------------------------------------------------------------------------
    # 4. D1, D2 and their raw partial derivatives
    # --------------------------------------------------------------------------
    t_pow_pm1 = t_int^p_minus_1
    s_pow_pm1 = s_int^p_minus_1
    t_pow_pm2 = t_int^p_minus_2
    s_pow_pm2 = s_int^p_minus_2
    s_pow_pm3 = s_int^p_minus_3

    D1 = interval(1) - s_int * t_pow_pm1
    D1_p  = -s_int * t_pow_pm1 * ln_t
    D1_s  = -t_pow_pm1
    D1_t  = -s_int * p_minus_1 * t_pow_pm2
    D1_ps = -t_pow_pm1 * ln_t
    D1_pt = -s_int * t_pow_pm2 * (p_minus_1 * ln_t + interval(1))
    D1_st = -p_minus_1 * t_pow_pm2
    D1_ss = interval(0)

    D2 = interval(1) - t_int * s_pow_pm1
    D2_p  = -t_int * s_pow_pm1 * ln_s
    D2_s  = -t_int * p_minus_1 * s_pow_pm2
    D2_t  = -s_pow_pm1
    D2_ps = -t_int * s_pow_pm2 * (p_minus_1 * ln_s + interval(1))
    D2_pt = -s_pow_pm1 * ln_s
    D2_st = -p_minus_1 * s_pow_pm2
    D2_ss = -t_int * p_minus_1 * p_minus_2 * s_pow_pm3

    # --------------------------------------------------------------------------
    # 5. H1 and H2 via B^(p-1), (sigma*A)^(p-1), (tau*A)^(p-1)
    # --------------------------------------------------------------------------
    SA = s_int * A
    SA_p  = s_int * A_p
    SA_s  = A + s_int * A_s
    SA_t  = s_int * A_t
    SA_ps = A_p + s_int * A_ps
    SA_pt = s_int * A_pt
    SA_st = A_t + s_int * A_st
    SA_ss = interval(2) * A_s + s_int * A_ss

    TA = t_int * A
    TA_p  = t_int * A_p
    TA_s  = t_int * A_s
    TA_t  = A + t_int * A_t
    TA_ps = t_int * A_ps
    TA_pt = A_p + t_int * A_pt
    TA_st = A_s + t_int * A_st
    TA_ss = t_int * A_ss

    # PB = B^(p-1)
    PB = B^p_minus_1
    PB_p  = PB * log(B) + p_minus_1 * B^p_minus_2 * B_p
    PB_s  = p_minus_1 * B^p_minus_2 * B_s
    PB_t  = p_minus_1 * B^p_minus_2 * B_t
    PB_ps = B^p_minus_2 * B_s * (interval(1) + p_minus_1 * log(B)) +
            p_minus_1 * p_minus_2 * B^p_minus_3 * B_p * B_s +
            p_minus_1 * B^p_minus_2 * B_ps
    PB_pt = B^p_minus_2 * B_t * (interval(1) + p_minus_1 * log(B)) +
            p_minus_1 * p_minus_2 * B^p_minus_3 * B_p * B_t +
            p_minus_1 * B^p_minus_2 * B_pt
    PB_st = p_minus_1 * (p_minus_2 * B^p_minus_3 * B_s * B_t + B^p_minus_2 * B_st)
    PB_ss = p_minus_1 * (p_minus_2 * B^p_minus_3 * B_s^interval(2) + B^p_minus_2 * B_ss)

    # PS = (sigma*A)^(p-1)
    PS = SA^p_minus_1
    PS_p  = PS * log(SA) + p_minus_1 * SA^p_minus_2 * SA_p
    PS_s  = p_minus_1 * SA^p_minus_2 * SA_s
    PS_t  = p_minus_1 * SA^p_minus_2 * SA_t
    PS_ps = SA^p_minus_2 * SA_s * (interval(1) + p_minus_1 * log(SA)) +
            p_minus_1 * p_minus_2 * SA^p_minus_3 * SA_p * SA_s +
            p_minus_1 * SA^p_minus_2 * SA_ps
    PS_pt = SA^p_minus_2 * SA_t * (interval(1) + p_minus_1 * log(SA)) +
            p_minus_1 * p_minus_2 * SA^p_minus_3 * SA_p * SA_t +
            p_minus_1 * SA^p_minus_2 * SA_pt
    PS_st = p_minus_1 * (p_minus_2 * SA^p_minus_3 * SA_s * SA_t + SA^p_minus_2 * SA_st)
    PS_ss = p_minus_1 * (p_minus_2 * SA^p_minus_3 * SA_s^interval(2) + SA^p_minus_2 * SA_ss)

    # PT = (tau*A)^(p-1)
    PT = TA^p_minus_1
    PT_p  = PT * log(TA) + p_minus_1 * TA^p_minus_2 * TA_p
    PT_s  = p_minus_1 * TA^p_minus_2 * TA_s
    PT_t  = p_minus_1 * TA^p_minus_2 * TA_t
    PT_ps = TA^p_minus_2 * TA_s * (interval(1) + p_minus_1 * log(TA)) +
            p_minus_1 * p_minus_2 * TA^p_minus_3 * TA_p * TA_s +
            p_minus_1 * TA^p_minus_2 * TA_ps
    PT_pt = TA^p_minus_2 * TA_t * (interval(1) + p_minus_1 * log(TA)) +
            p_minus_1 * p_minus_2 * TA^p_minus_3 * TA_p * TA_t +
            p_minus_1 * TA^p_minus_2 * TA_pt
    PT_st = p_minus_1 * (p_minus_2 * TA^p_minus_3 * TA_s * TA_t + TA^p_minus_2 * TA_st)
    PT_ss = p_minus_1 * (p_minus_2 * TA^p_minus_3 * TA_s^interval(2) + TA^p_minus_2 * TA_ss)

    H1 = PB + PS
    H1_p  = PB_p  + PS_p
    H1_s  = PB_s  + PS_s
    H1_t  = PB_t  + PS_t
    H1_ps = PB_ps + PS_ps
    H1_pt = PB_pt + PS_pt
    H1_st = PB_st + PS_st
    H1_ss = PB_ss + PS_ss

    H2 = PB - PT
    H2_p  = PB_p  - PT_p
    H2_s  = PB_s  - PT_s
    H2_t  = PB_t  - PT_t
    H2_ps = PB_ps - PT_ps
    H2_pt = PB_pt - PT_pt
    H2_st = PB_st - PT_st
    H2_ss = PB_ss - PT_ss

    # --------------------------------------------------------------------------
    # 6. Raw partial derivatives of g = C1*D1*H1 - C2*D2*H2
    # --------------------------------------------------------------------------
    # U = C1*D1*H1
    U_s = C1_s * D1 * H1 + C1 * D1_s * H1 + C1 * D1 * H1_s

    U_ps = C1_ps * D1 * H1 + C1_p * D1_s * H1 + C1_p * D1 * H1_s +
           C1_s * D1_p * H1 + C1 * D1_ps * H1 + C1 * D1_p * H1_s +
           C1_s * D1 * H1_p + C1 * D1_s * H1_p + C1 * D1 * H1_ps

    U_pt = C1_pt * D1 * H1 + C1_p * D1_t * H1 + C1_p * D1 * H1_t +
           C1_t * D1_p * H1 + C1 * D1_pt * H1 + C1 * D1_p * H1_t +
           C1_t * D1 * H1_p + C1 * D1_t * H1_p + C1 * D1 * H1_pt

    U_st = C1_st * D1 * H1 + C1_s * D1_t * H1 + C1_s * D1 * H1_t +
           C1_t * D1_s * H1 + C1 * D1_st * H1 + C1 * D1_s * H1_t +
           C1_t * D1 * H1_s + C1 * D1_t * H1_s + C1 * D1 * H1_st

    U_ss = C1_ss * D1 * H1 + C1 * D1_ss * H1 + C1 * D1 * H1_ss +
           interval(2) * (C1_s * D1_s * H1 + C1_s * D1 * H1_s + C1 * D1_s * H1_s)

    # V = C2*D2*H2
    V_s = C2_s * D2 * H2 + C2 * D2_s * H2 + C2 * D2 * H2_s

    V_ps = C2_ps * D2 * H2 + C2_p * D2_s * H2 + C2_p * D2 * H2_s +
           C2_s * D2_p * H2 + C2 * D2_ps * H2 + C2 * D2_p * H2_s +
           C2_s * D2 * H2_p + C2 * D2_s * H2_p + C2 * D2 * H2_ps

    V_pt = C2_pt * D2 * H2 + C2_p * D2_t * H2 + C2_p * D2 * H2_t +
           C2_t * D2_p * H2 + C2 * D2_pt * H2 + C2 * D2_p * H2_t +
           C2_t * D2 * H2_p + C2 * D2_t * H2_p + C2 * D2 * H2_pt

    V_st = C2_st * D2 * H2 + C2_s * D2_t * H2 + C2_s * D2 * H2_t +
           C2_t * D2_s * H2 + C2 * D2_st * H2 + C2 * D2_s * H2_t +
           C2_t * D2 * H2_s + C2 * D2_t * H2_s + C2 * D2 * H2_st

    V_ss = C2_ss * D2 * H2 + C2 * D2_ss * H2 + C2 * D2 * H2_ss +
           interval(2) * (C2_s * D2_s * H2 + C2_s * D2 * H2_s + C2 * D2_s * H2_s)

    g_s  = U_s  - V_s
    g_ps = U_ps - V_ps
    g_pt = U_pt - V_pt
    g_st = U_st - V_st
    g_ss = U_ss - V_ss

    # --------------------------------------------------------------------------
    # 7. Final branch mixed derivative
    # --------------------------------------------------------------------------
    branch_g_pt = g_pt + g_ps * sigma_tau + g_st * sigma_p +
                  g_ss * sigma_p * sigma_tau + g_s * sigma_pt

    return branch_g_pt
end





const p_min_global = 1 // 1
const p_max_global = 101 // 100
const p_step = 1 // 2000

const tau_start = 1 // 5
const tau_end = 1 // 3
const tau_step = 1 // 2000

const sigma_fixed_start = 1 // 1
const sigma_fixed_end = 101377243 // 100000000

const tp_guess = interval(0, 36//100)
const target_precision = 1 // 10^9
const deriv_threshold = 1 // 10^9
const tau_subdivisions = 10

println("Calculating branch g_{p tau} over the p=1 direct-tau grid.")
println("  p range:     [$p_min_global, $p_max_global] with step $p_step")
println("  tau range:   [$tau_start, $tau_end] with step $tau_step")
println("  sigma range: [$sigma_fixed_start, $sigma_fixed_end] fixed on every box")
println("  Inequality to check: branch g_{p tau} > $deriv_threshold")
println("="^40)
println("-"^150)
println("  P-Rect             S-Enclosure       tau-interval      t_p-interval      branch g_{p tau} Bound")
println("-"^150)

start_time = time()

failed_p_intervals = Interval{Float64}[]
checked_leaf_boxes = 0
subdivision_count = 0
worst_bound = Inf

p_values_asc = collect(p_min_global:p_step:p_max_global)
if last(p_values_asc) < p_max_global
    push!(p_values_asc, p_max_global)
end

tau_values_asc = collect(tau_start:tau_step:tau_end)
if last(tau_values_asc) < tau_end
    push!(tau_values_asc, tau_end)
end

for i in 1:(length(p_values_asc) - 1)
    p_lo = p_values_asc[i]
    p_hi = p_values_asc[i+1]
    p_interval = interval(p_lo, p_hi)

    tp_result = find_tp_range(p_interval, tp_guess, target_precision, false)
    if isempty_interval(tp_result)
        println("!!! FAILED to find t_p for P-Interval: ", p_interval)
        push!(failed_p_intervals, p_interval)
        continue
    end

    s_interval = interval(sigma_fixed_start, sigma_fixed_end)
    local p_interval_failed = false

    for j in 1:(length(tau_values_asc) - 1)
        t_lo = tau_values_asc[j]
        t_hi = tau_values_asc[j+1]
        tau_interval = interval(t_lo, t_hi)

        deriv_res = calculate_branch_gptau_interval(p_interval, s_interval, tau_interval)

        if isempty_interval(deriv_res)
            println("!!! FAIL: Empty derivative interval for P=$p_interval, T=$tau_interval")
            p_interval_failed = true
            break
        end

        println("  P=$p_interval S=$s_interval T=$tau_interval => tp=$tp_result branch g_{p tau}=$deriv_res")

        if inf(deriv_res) > deriv_threshold
            global checked_leaf_boxes, worst_bound
            checked_leaf_boxes += 1
            worst_bound = min(worst_bound, inf(deriv_res))
            continue
        else
            println("--- Check failed on full tau interval. Subdividing... ---")
            global subdivision_count
            subdivision_count += 1

            tau_subs = subdivide_tau_interval(tau_interval, tau_subdivisions)
            all_subs_passed = true

            for (k, sub_tau) in enumerate(tau_subs)
                sub_res = calculate_branch_gptau_interval(p_interval, s_interval, sub_tau)

                if isempty_interval(sub_res)
                    println("      !!! SUB-FAIL: Empty derivative interval at P=$p_interval, T[$k]=$sub_tau")
                    all_subs_passed = false
                    break
                end

                println("      (Sub) P=$p_interval S=$s_interval T[$k]=$sub_tau => branch g_{p tau}=$sub_res")

                if inf(sub_res) > deriv_threshold
                    global checked_leaf_boxes, worst_bound
                    checked_leaf_boxes += 1
                    worst_bound = min(worst_bound, inf(sub_res))
                else
                    println("      !!! SUB-FAIL: P=$p_interval, T[$k]=$sub_tau")
                    println("          branch g_{p tau} = $sub_res (Wanted lo > $deriv_threshold)")
                    all_subs_passed = false
                    break
                end
            end

            if all_subs_passed
                continue
            else
                println("!!! FAIL: P=$p_interval, T=$tau_interval")
                println("    branch g_{p tau} = $deriv_res (Wanted lo > $deriv_threshold)")
                p_interval_failed = true
                break
            end
        end
    end

    if p_interval_failed
        push!(failed_p_intervals, p_interval)
    end
end

end_time = time()
println("="^40)
println("Final Report on Inconclusive P-Intervals:")

if isempty(failed_p_intervals)
    println("SUCCESS: All p-intervals were verified successfully.")
else
    merged_failed = merge_intervals(failed_p_intervals)
    println("The following p-ranges failed the check:")
    for p_int in merged_failed
        println("  ", p_int)
    end
end

println("Checked certified leaf boxes: ", checked_leaf_boxes)
println("Number of tau subdivisions used: ", subdivision_count)
println("Worst certified branch g_{p tau} lower bound: ", worst_bound)
println("Total time: ", round(end_time - start_time, digits=4), " seconds.")
