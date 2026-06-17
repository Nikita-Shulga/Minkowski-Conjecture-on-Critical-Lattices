using IntervalArithmetic



"""
    get_L_factorized(val_int, p_int)
    Calculates L = 1/p * [ ln(1+x^p)/p - x^p*ln(x)/(1+x^p) ].
"""
function get_L_factorized(val_int::Interval{<:Real}, p_int::Interval{<:Real})
    val_pow_p = val_int^p_int
    ln_one_plus = log(interval(1) + val_pow_p)
    
    term1 = ln_one_plus / p_int
    term_x_ln_x = t_pow_p_log_t_int(val_int, p_int)
    term2 = term_x_ln_x / (interval(1) + val_pow_p)
    
    return (term1 - term2) / p_int
end

"""
    t_pow_p_log_t_int(t_int::Interval{<:Real}, p_int::Interval{<:Real})

Calculates t^p * ln(t).
1. Neighborhood of 0: If t is within [0, e^(-1/p)], the function is strictly decreasing.
   We map [0, b] -> [b^p ln(b), 0].
2. Safe Region: If t is away from 0, use standard arithmetic.
"""
function t_pow_p_log_t_int(t_int::Interval{<:Real}, p_int::Interval{<:Real})
    if inf(t_int) > 0
        return (t_int^p_int) * log(t_int)
    end

    # Handle Neighborhood of 0
    safe_monotonic_bound = 1 // 10
    
    if sup(t_int) < safe_monotonic_bound
        # Domain: [0, t_hi] -> Range: [t_hi^p * ln(t_hi), 0]
        t_hi_iv = interval(sup(t_int))
        val_at_end = (t_hi_iv^p_int) * log(t_hi_iv)
        return interval(inf(val_at_end), 0)
    else
        return (t_int^p_int) * log(t_int)
    end
end

"""
    get_T2(val_int, p_int)
    Calculates explicit second partial T2(x) using the SAFE expanded form.
"""
function get_T2(val_int::Interval{<:Real}, p_int::Interval{<:Real})
    val_pow_p = val_int^p_int
    ln_one_plus = log(interval(1) + val_pow_p)
    one_plus_sq = (interval(1) + val_pow_p)
    
    # 1. First Term: (2 * x^p * ln x) / p^2 (1+..)
    term_x_ln_x = t_pow_p_log_t_int(val_int, p_int)
    term1 = (interval(2) * term_x_ln_x) / (p_int^interval(2) * one_plus_sq)
    
    # 2. Second Term: - (2 * ln(1+x^p)) / p^3
    term2 = (interval(2) * ln_one_plus) / (p_int^interval(3))
    
    # 3. Third Term: - x^p * (ln x)^interval(2) / ( p * (1+x^p)^interval(2) )
    p_half = p_int / interval(2)
    safe_sqrt_term = interval(0)
    
    # Check monotonicity for x^(p/2) ln x. Critical point exp(-2) ~ 0.135
    if sup(val_int) < 13 // 100
         safe_sqrt_term = t_pow_p_log_t_int(val_int, p_half)
    else
         # Fallback
         safe_sqrt_term = (val_int^p_half) * log(val_int)
    end

    term3 = (safe_sqrt_term^interval(2)) / (p_int * one_plus_sq^2)
    
    return term1 - term2 - term3
end

"""
    get_Psi(X, X_p, X_pp, p_int)
    Calculates Psi(X, X_p, X_pp) = d^2/dp^2 (X^p).
"""
function get_Psi(X, X_p, X_pp, p_int)
    X_pow_p = X^p_int
    ln_X = log(X)
    
    # Term 1: Squared bracket
    inner_sq = ln_X + (p_int * X_p / X)
    term_sq = inner_sq^interval(2)
    
    # Term 2: Linear bracket
    inner_lin = interval(2) * X_p + p_int * X_pp - (p_int * (X_p^interval(2)) / X)
    term_lin = inner_lin / X
    
    return X_pow_p * (term_sq + term_lin)
end



"""
    find_tp_range(p_range, initial_tp_guess, precision, verbose)

Iteratively finds the interval for t_p given an interval for p.
"""
function find_tp_range(
    p_range::Interval{<:Real}, # Accept any Real-based Interval
    initial_tp_guess::Interval{<:Real},
    precision::Real,           # Accept any Real (Float or Rational)
    verbose::Bool = false
)
    tp_interval = initial_tp_guess
    one_over_p = interval(1) / p_range
    two_pow_one_over_p = interval(2)^one_over_p

    for i in 1:100 # Max 100 iterations
        if verbose
            println("  t_p Iteration $i: ", tp_interval)
        end
        
        if isempty_interval(tp_interval)
            println("  Warning: t_p interval is empty.")
            return emptyinterval(Float64)
        end

        tp_interval = intersect_interval(tp_interval, interval(1//10, 36//100); dec=:auto)
       
        # Calculate the fixed point map: F(t) = 1 - 2^(-1/p) * (1 + t^p)^(1/p)
        tp_pow_p = tp_interval^p_range
        one_plus_tp_pow_p = interval(1) + tp_pow_p
        term_pow = one_plus_tp_pow_p^one_over_p
        full_term = term_pow / two_pow_one_over_p
        new_tp_interval = interval(1) - full_term
        
        # Check for convergence (width check)
        if diam(new_tp_interval) < precision
            if verbose
                println("  Converged t_p to $new_tp_interval")
            end
            return new_tp_interval
        end
        
        tp_interval = intersect_interval(tp_interval, new_tp_interval; dec=:auto)
    end
    
    if verbose
        println("  t_p max iterations reached. Final interval: $tp_interval")
    end
    return tp_interval
end

function get_tp_interval(p_int::Interval{<:Real})
    # Wrapped call to the new find_tp_range logic
    return find_tp_range(p_int, interval(0 // 1, 35 // 100), 1 // 10^13, false)
end

function calculate_d2Delta_p1_dp2(p_int::Interval{<:Real})
    t = get_tp_interval(p_int)

    
    one_minus_t = interval(1) - t
    
    ln_t = log(t)
    ln_one_minus_t = log(one_minus_t)
    
    # Derivatives of Theta 
    Theta_p = interval(2) * (one_minus_t^p_int) * ln_one_minus_t - (t^p_int) * ln_t
    Theta_t = -p_int * ( interval(2) * one_minus_t^(p_int - interval(1)) + t^(p_int - interval(1)) )
    
    
    dt_dp = -Theta_p / Theta_t

    
    Theta_pp = interval(2) * (one_minus_t^p_int) * (ln_one_minus_t^interval(2)) - (t^p_int) * (ln_t^interval(2))
    
    term_pt_1 = -interval(2) * (one_minus_t^(p_int - interval(1))) * (interval(1) + p_int * ln_one_minus_t)
    term_pt_2 = -(t^(p_int - interval(1))) * (interval(1) + p_int * ln_t)
    Theta_pt = term_pt_1 + term_pt_2
    
    term_tt_inner = interval(2) * (one_minus_t^(p_int - interval(2))) - (t^(p_int - interval(2)))
    Theta_tt = p_int * (p_int - interval(1)) * term_tt_inner
    
    num_d2t = Theta_pp + interval(2) * Theta_pt * dt_dp + Theta_tt * (dt_dp^interval(2))
    d2t_dp2 = -num_d2t / Theta_t
    
    # Base Delta(p,1) 
    val_4 = interval(4)
    term_4_pow = val_4^(-interval(1)/p_int) 
    Delta_p1 = term_4_pow * (interval(1) + t) / (interval(1) - t)
    
    ln4 = log(interval(4))
    
    term_A_inner = (ln4 / (p_int^interval(2))) + (interval(2) / (interval(1) - t^interval(2))) * dt_dp
    term_A = term_A_inner^interval(2)
    term_B = -(interval(2) * ln4) / (p_int^interval(3))
    
    num_C = (interval(1) - t^interval(2)) * d2t_dp2 + interval(2) * t * (dt_dp^interval(2))
    den_C = (interval(1) - t^interval(2))^interval(2)
    term_C = interval(2) * num_C / den_C
    
    final_ref = Delta_p1 * (term_A + term_B + term_C)
    
    return final_ref
end






function calculate_l1_pp_interval(
    p_int::Interval{<:Real},
    s_int::Interval{<:Real},
    t_int::Interval{<:Real}
)
    val_1_sp = interval(1) + s_int ^p_int
    val_1_tp = interval(1) + t_int^p_int
    
    C1 = val_1_sp^(-interval(1) / p_int)
    C2 = val_1_tp^(-interval(1) / p_int)
    
    A = C2 - C1
    B = t_int * C2 + s_int  * C1
    
    K_sigma = val_1_sp^(-interval(1) - interval(1)/p_int)
    Delta = (t_int + s_int ) * C1 * C2 

    # Gradient Factors
    L_sigma = get_L_factorized(s_int , p_int)
    L_tau   = get_L_factorized(t_int, p_int)

    M = (s_int ^(p_int - interval(1))) / val_1_sp - interval(1) / (s_int  + t_int)
     
    A_p = C2 * L_tau - C1 * L_sigma
    B_p = t_int * C2 * L_tau + s_int  * C1 * L_sigma
    
    # Implicit Constraints 
    H1 = (B^(p_int - interval(1))) + (s_int ^(p_int - interval(1))) * (A^(p_int - interval(1)))
    Denom_D = p_int * K_sigma * H1
    
    term_Np_1 = (A^p_int) * log(A)
    term_Np_2 = (B^p_int) * log(B)
    term_Np_3 = p_int * ( (A^(p_int - interval(1))) * A_p + (B^(p_int - interval(1))) * B_p )
    N_p = term_Np_1 + term_Np_2 + term_Np_3
    
    if 0 in Denom_D
        error("RIGOR FAIL: Denominator of sigma_p contains 0.")
    else
        sigma_p = -N_p / Denom_D
    end

    
    Q = L_sigma + L_tau - M * sigma_p

    # Second Derivative Factors 
    T2_sigma = get_T2(s_int , p_int)
    T2_tau   = get_T2(t_int, p_int)
    L_pp = T2_sigma + T2_tau
    
    T_cal_sigma = (L_sigma^interval(2)) + T2_sigma
    T_cal_tau   = (L_tau^interval(2))   + T2_tau
    
    A_pp = C2 * T_cal_tau - C1 * T_cal_sigma
    B_pp = t_int * C2 * T_cal_tau + s_int  * C1 * T_cal_sigma
    
    M_p = (s_int ^(p_int - interval(1)) * log(s_int )) / (val_1_sp^interval(2))
    term_Ms_1 = (s_int ^(p_int - interval(2))) * ((p_int - interval(1)) - s_int ^p_int) / (val_1_sp^interval(2))
    term_Ms_2 = interval(1) / ((s_int  + t_int)^interval(2))
    M_sigma = term_Ms_1 + term_Ms_2
    
    Psi_A = get_Psi(A, A_p, A_pp, p_int)
    Psi_B = get_Psi(B, B_p, B_pp, p_int)

    
    F_pp = Psi_A + Psi_B
    
    ln_Ks_p = log(val_1_sp)/(p_int^interval(2)) - ((p_int + interval(1))/p_int) * (s_int ^p_int * log(s_int ) / val_1_sp)
    term_H1p_B = (B^(p_int - interval(1))) * (log(B) + (p_int - interval(1))*B_p/B)
    term_sA = s_int  * A
    term_H1p_sA = (term_sA^(p_int - interval(1))) * (log(term_sA) + (p_int - interval(1))*A_p/A)
    H1_p = term_H1p_B + term_H1p_sA
    F_psigma = Denom_D * ( interval(1)/p_int + ln_Ks_p + H1_p/H1 )
    
    ln_Ks_sigma = -((p_int + interval(1)) * s_int ^(p_int - interval(1))) / val_1_sp
    term_H1s_1 = (B^(p_int - interval(2))) * K_sigma
    term_H1s_2 = (s_int ^(interval(2)*p_int - interval(2))) * (A^(p_int - interval(2))) * K_sigma
    term_H1s_3 = (s_int ^(p_int - interval(2))) * (A^(p_int - interval(1)))
    H1_sigma = (p_int - interval(1)) * (term_H1s_1 + term_H1s_2 + term_H1s_3)
    F_sigma_sigma = Denom_D * ( ln_Ks_sigma + H1_sigma / H1 ) 
    
    num_spp = F_pp + sigma_p * ( interval(2) * F_psigma + sigma_p * F_sigma_sigma )
    sigma_pp = -num_spp / Denom_D

    
    term_bracket = (Q^interval(2)) + L_pp - sigma_p * (interval(2) * M_p + sigma_p * M_sigma) - M * sigma_pp 
                   
    d2Delta_dp2 = Delta * term_bracket
    ref_term_val = calculate_d2Delta_p1_dp2(p_int)
    
    return (d2Delta_dp2 - ref_term_val), d2Delta_dp2, ref_term_val
end

"""
    merge_intervals(intervals::Vector{Interval{Float64}})

Merges a list of sorted, contiguous, or overlapping intervals into a minimal list of disjoint intervals.
"""
function merge_intervals(intervals::Vector{Interval{Float64}})
    if isempty(intervals)
        return Interval{Float64}[]
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





# Domain Constants
const p_min_global = 1 // 1
const p_max_global = 10001 // 10000
const p_step = 7 // 10000

const sigma_min_global = 1 // 1
const sigma_max_global = 101377243 // 100000000
const sigma_step = 2 // 10000


const tau_start = 0 // 1
const tau_end =  2 // 1000
const tau_step = 2 // 10000

const k_threshold = - 1 // 10^9

println("Calculating (l^(0))''_{pp} with DETAILED OUTPUT and P-FAIL LOGIC.")
println("  p range:     [$p_min_global, $p_max_global] with step $p_step")
println("  sigma range: [$sigma_min_global, $sigma_max_global] with step $sigma_step")
println("  tau range:   [$tau_start, $tau_end] with step $tau_step")
println("  Inequality to check: l1_pp < $k_threshold")
println("="^40)
println("-"^150)
println("  P-Rect             S-Rect             τ-sub-interval     d2Delta          RefTerm           l1_pp Bound")
println("-"^150)

start_time = time()

num_p_steps = Int(ceil((p_max_global - p_min_global) / p_step))
num_s_steps = Int(ceil((sigma_max_global - sigma_min_global) / sigma_step))
num_tau_steps = Int(round((tau_end - tau_start) / tau_step))

failed_p_intervals = Interval{Float64}[]

# Outer Loop: P
for i_p in 0:(num_p_steps - 1)
    p_lo = interval(p_min_global + i_p * p_step)
    p_hi = interval(min(p_min_global + (i_p + 1) * p_step, p_max_global))
    if p_hi <= p_lo; continue; end
    
    p_interval = interval(p_lo, p_hi)
    local p_interval_failed = false # Flag for this P-interval
    
    # Middle Loop: Sigma
    for i_s in 0:(num_s_steps - 1)
        s_lo = sigma_min_global + i_s * sigma_step
        s_hi = min(sigma_min_global + (i_s + 1) * sigma_step, sigma_max_global)
        if s_hi <= s_lo; continue; end
        s_interval = interval(s_lo, s_hi)
        
        # Inner Loop: Tau
        for i_t in 0:(num_tau_steps - 1)
            t_lo = tau_start + i_t * tau_step
            t_hi = tau_start + (i_t + 1) * tau_step
            tau_interval = interval(t_lo, t_hi)
             
            # CALCULATION
            (res_val, d2_val, ref_val) = calculate_l1_pp_interval(p_interval, s_interval, tau_interval)
            
            # PRINT EVERY STEP (Detailed mode) [cite: 32]
            println("  P=$p_interval S=$s_interval T=$tau_interval => d2=$d2_val ref=$ref_val l1_pp=$res_val")
            
            # CHECK FAILURE [cite: 33]
            if sup(res_val) >= k_threshold
                # If ANY sub-interval fails, the whole P-block is suspect.
                println("      !!! FAIL: T=$tau_interval (Hi >= $k_threshold)")
                
                # Mark failure
                p_interval_failed = true
                
                # Break out of Tau loop
                break 
            end
        end
        
        # If failure occurred in Tau loop, break out of Sigma loop too
        if p_interval_failed
            break
        end
    end
     
    # If P interval failed, store it
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

println("Total time: ", round(end_time - start_time, digits=4), " seconds.")