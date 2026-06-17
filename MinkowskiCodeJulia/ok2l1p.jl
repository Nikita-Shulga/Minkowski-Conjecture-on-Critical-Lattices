using IntervalArithmetic

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
        
        #
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

"""
    find_new_tau_range(p_int, s_int, t_up_val, precision, verbose)

Iteratively finds the bounds for the new τ variable.
"""
function find_new_tau_range(
    p_int::Interval{<:Real},
    s_int::Interval{<:Real},
    t_up_val::Real,
    precision::Real,
    verbose::Bool = false
)
    # ---------------------------------------------------------
    # PHASE 1: CANDIDATE GENERATION (Float64 heuristics)
    # ---------------------------------------------------------
    
    # Use scalar bounds
    op_scalar = interval(sup(p_int))
    up_scalar = interval(inf(p_int))
    os_scalar = interval(sup(s_int))
    us_scalar = interval(inf(s_int))
    
    utau = interval(0)
    otau = interval(t_up_val)   

    if verbose
        println("    Starting new τ iteration with p=$p_int, s=$s_int")
    end

    # --- Loop 1: Find Lower Bound Guess (utau) ---
    for i in 1:100
        utau_i = utau
        
        # Just calculate normally
        term1_l = (interval(1) + utau_i^op_scalar)^(interval(1) / op_scalar)
        term_a_l = (interval(1) + utau_i^op_scalar)^(-interval(1) / op_scalar)
        term_b_l = (interval(1) + os_scalar^up_scalar)^(-interval(1) / up_scalar) 
        base_l = max(interval(0), term_a_l - term_b_l) 
        term_c_l = base_l^up_scalar
        base2_l = max(interval(0), interval(1) - term_c_l) 
        term2_l = base2_l^(interval(1) / up_scalar)
        term3_l = os_scalar * (interval(1) + os_scalar^op_scalar)^(-interval(1) / op_scalar) 
        
        utau_new = term1_l * (term2_l - term3_l)

        if isnan(utau_new) || inf(utau_new) < 0
   		 utau_new = interval(0)
	end
        if sup(abs(utau_new - utau_i)) < precision
    		utau = utau_new
    		break
	end
        utau = utau_new
    end

    # --- Loop 2: Find Upper Bound Guess (otau) ---
    # REMOVED setrounding(Float64, RoundUp)
    for i in 1:100
        otau_i = otau
        
        # Just calculate normally
        term1_u = (interval(1) + otau_i^up_scalar)^(interval(1) / up_scalar)
        term_a_u = (interval(1) + otau_i^up_scalar)^(-interval(1) / up_scalar)
        term_b_u = (interval(1) + us_scalar^op_scalar)^(-interval(1) / op_scalar) 
        base_u = max(interval(0), term_a_u - term_b_u) 
        term_c_u = base_u^op_scalar
        base2_u = max(interval(0), interval(1) - term_c_u) 
        term2_u = base2_u^(interval(1) / op_scalar)
        term3_u = us_scalar * (interval(1) + us_scalar^up_scalar)^(-interval(1) / up_scalar) 
        
        otau_new = term1_u * (term2_u - term3_u)

        if isnan(otau_new); otau_new = otau_i; end
        if sup(abs(otau_new - otau_i)) < precision
    		otau = otau_new
    		break
	end
        otau = otau_new
    end

    # --- Float64 Safety Check ---
    safe_upper_limit = min(interval(t_up_val), interval(1)/interval(3))
    if sup(otau) > inf(safe_upper_limit)
    	otau = safe_upper_limit
    end

    if inf(otau) < sup(utau)
        println("    !!! RIGOR FAIL: Float64 bounds crossed (utau > otau).")
        println("      utau: $utau")
        println("      otau: $otau")
        return emptyinterval(Float64)
    end

    # ---------------------------------------------------------
    # PHASE 2: RIGOROUS VERIFICATION (Interval Arithmetic)
    # ---------------------------------------------------------
    
    # 1. Epsilon Inflation (Padding)
    # This step handles the tiny errors introduced by removing setrounding
    inflation_factor = interval(1 + 1e-5)
    epsilon_pad = interval(1e-7)

    utau_padded = max(interval(0), utau / inflation_factor - epsilon_pad)
    otau_padded = otau * inflation_factor + epsilon_pad
    
    candidate_I = interval(utau_padded, otau_padded)

    # 2. Convert parameters to Intervals
    op_I = interval(op_scalar)
    up_I = interval(up_scalar)
    os_I = interval(os_scalar)
    us_I = interval(us_scalar)

    # 3. Rigorous Evaluation
    utau_check_I = try
        term1_l = (interval(1) + candidate_I^op_I)^(interval(1) / op_I)
        term_a_l = (interval(1) + candidate_I^op_I)^(interval(-1) / op_I)
        term_b_l = (interval(1) + os_I^up_I)^(-interval(1) / up_I)
        base_l = max(interval(0), term_a_l - term_b_l) 
        term_c_l = base_l^up_I
        base2_l = max(interval(0), interval(1) - term_c_l)
        term2_l = base2_l^(interval(1) / up_I)
        term3_l = os_I * (interval(1) + os_I^op_I)^(interval(-1) / op_I)
        raw_res = term1_l * (term2_l - term3_l)
        max(interval(0), raw_res)
    catch
        return emptyinterval(Float64)
    end

    otau_check_I = try
        term1_u = (interval(1) + candidate_I^up_I)^(interval(1) / up_I)
        term_a_u = (interval(1) + candidate_I^up_I)^(-interval(1) / up_I)
        term_b_u = (interval(1) + us_I^op_I)^(-interval(1) / op_I)
        base_u = max(interval(0), term_a_u - term_b_u)
        term_c_u = base_u^op_I
        base2_u = max(interval(0), interval(1) - term_c_u)
        term2_u = base2_u^(interval(1) / op_I)
        term3_u = us_I * (interval(1) + us_I^up_I)^(-interval(1) / up_I)
        raw_res_u = term1_u * (term2_u - term3_u)
        max(interval(0), raw_res_u)
    catch
        return emptyinterval(Float64)
    end

    # 4. Construct Verified Result with SAFETY CLAMP
    raw_verified_I = interval(inf(utau_check_I), sup(otau_check_I))
    
    # The Safety Clamp:
    clamped_verified_I = intersect_interval(raw_verified_I, interval(0, safe_upper_limit); dec=:auto)

    # 5. Proof Check (Subset)
    if issubset_interval(clamped_verified_I, candidate_I)
        # Return the tightest valid interval
        return intersect_interval(clamped_verified_I, candidate_I; dec=:auto)
    else
        println("    !!! RIGOR FAIL: Verified interval is NOT inside candidate.")
        println("      Float64 Candidate (Pre-Check): $(interval(utau, otau))")
        println("      Padded Candidate  (Input to I): $candidate_I")
        println("      Rigorous Check    (Output of I): $raw_verified_I")
        println("      Clamped Check     (Safe Limit):  $clamped_verified_I")
        return emptyinterval(Float64)
    end
end




"""
    calculate_l1_derivative_interval(p_int, s_int, t_int, tp_int)

Calculates the derivative (l^(1))'_p.
Rigorous version:
1. Uses t_pow_p_log_t_safe for the t=0 singularity.
2. Direct calculation for A and B (assuming strictly separated tau and sigma).
"""
function calculate_l1_derivative_interval(
    p_int::Interval{<:Real},
    s_int::Interval{<:Real},
    t_int::Interval{<:Real},
    tp_int::Interval{<:Real}
)
    # --- 1. Common Helper Terms ---
    one_over_p = interval(1) / p_int
    p_minus_1 = p_int - interval(1)
    
    ln_s = log(s_int)
    
    val_1_sp = interval(1) + s_int^p_int
    val_1_tp = interval(1) + t_int^p_int 
    
    # C1, C2
    C1 = val_1_sp^(-one_over_p)
    C2 = val_1_tp^(-one_over_p)
    
    # A, B
    # Analytically A > 0 because tau < sigma.
    # We rely on the solver ensuring tau and sigma are distinct enough.
    A = C2 - C1
    B = t_int * C2 + s_int * C1
    
    # K_tau
    K_tau = val_1_tp^(-interval(1) - one_over_p)
    
    t_pow_p = t_int^p_int
    s_pow_p = s_int^p_int
    
    # H2 = B^(p-1) - t^(p-1) * A^(p-1)
    # Note: If A is accidentally non-positive due to float noise, A^(p-1) will error here.
    H2 = (B^p_minus_1) - (t_int^p_minus_1) * (A^p_minus_1)
    
    # --- 2. Calculate L_sigma, L_tau ---
    
    L_sigma = log(val_1_sp) / (p_int^2) - (s_pow_p * ln_s) / (p_int * val_1_sp)
    
    # L_tau: RIGOROUS CALCULATION (Keep this one!)
    # Handles the t=0 singularity via monotonicity
    term_singularity = t_pow_p_log_t_safe(t_int, p_int)
    L_tau = log(val_1_tp) / (p_int^2) - term_singularity / (p_int * val_1_tp)
    
    # --- 3. Calculate Partial derivatives of A and B wrt p ---
    dA_dp = C2 * L_tau - C1 * L_sigma
    dB_dp = t_int * C2 * L_tau + s_int * C1 * L_sigma
    
    # --- 4. Calculate N_p ---
    # N_p = A^p ln A + B^p ln B + ...
    term_N1 = (A^p_int) * log(A)
    term_N2 = (B^p_int) * log(B)
    term_N3 = p_int * (A^p_minus_1) * dA_dp
    term_N4 = p_int * (B^p_minus_1) * dB_dp
    
    N_p = term_N1 + term_N2 + term_N3 + term_N4
    
    # --- 5. Calculate T_p ---
    
    if 0.0 in H2
        return interval(-1e10, 1e10) # Division by zero fail-safe
    end
    
    T_p = - N_p / (p_int * K_tau * H2)
    
    # --- 6. Calculate dDelta(p, sigma) / dp ---
    Delta_ps = (t_int + s_int) * C2 * C1
    
    explicit_tau_term = (interval(1) / (s_int + t_int)) - ( (t_int^p_minus_1) / val_1_tp )
    
    dDelta_ps_dp = Delta_ps * ( (L_sigma + L_tau) + explicit_tau_term * T_p )
    
    # --- 7. Calculate dDelta(p, 1) / dp ---
    term_4_pow = interval(4)^(-one_over_p)
    term_tp_frac = (interval(1) + tp_int) / (interval(1) - tp_int)
    Delta_p1 = term_4_pow * term_tp_frac
    
    # Calculate dt_p / dp
    one_minus_tp = interval(1) - tp_int
    
    num_dtp_1 = interval(2) * (one_minus_tp^p_int) * log(one_minus_tp)
    num_dtp_2 = (tp_int^p_int) * log(tp_int)
    num_dtp = num_dtp_1 - num_dtp_2
    
    den_dtp_inner = interval(2) * (one_minus_tp^p_minus_1) + (tp_int^p_minus_1)
    den_dtp = p_int * den_dtp_inner
    
    dtp_dp = num_dtp / den_dtp
    
    bracket_term_1 = log(interval(4)) / (p_int^2)
    bracket_term_2 = (interval(2) / (interval(1) - tp_int^2)) * dtp_dp
    
    dDelta_p1_dp = Delta_p1 * (bracket_term_1 + bracket_term_2)
    
    # --- 8. Final Result ---
    l1_prime_p = dDelta_ps_dp - dDelta_p1_dp
    
    return l1_prime_p
end



"""
    t_pow_p_log_t_safe(t_int::Interval{Float64}, p_int::Interval{Float64})

Calculates t^p * ln(t).
1. Neighborhood of 0: If t is within [0, e^(-1/p)], the function is strictly decreasing.
   We map [0, b] -> [b^p ln(b), 0].
2. Safe Region: If t is away from 0, use standard arithmetic.
"""
function t_pow_p_log_t_safe(t_int::Interval{<:Real}, p_int::Interval{<:Real})
    if inf(t_int) > 0.0
        return (t_int^p_int) * log(t_int)
    end

    # Handle Neighborhood of 0
    safe_monotonic_bound = 1 // 10
    
    if sup(t_int) < safe_monotonic_bound
        # Domain: [0, t_hi] -> Range: [t_hi^p * ln(t_hi), 0]
        t_hi_iv = interval(sup(t_int))
        val_at_end = (t_hi_iv^p_int) * log(t_hi_iv)
        return interval(inf(val_at_end), 0.0)
    else
        return (t_int^p_int) * log(t_int)
    end
end




"""
    merge_intervals(intervals::Vector{Interval{Float64}})

Merges a list of sorted, contiguous, or overlapping intervals into a minimal list of disjoint intervals.
"""
function merge_intervals(intervals::Vector{Interval{Float64}})
    # FIX: Return empty list immediately if input is empty to prevent BoundsError
    if isempty(intervals)
        return Interval{Float64}[]
    end
    
    # Sort the intervals by their lower bound to ensure they are in order
    sorted_intervals = sort(intervals, by = x -> inf(x))
    
    merged = [sorted_intervals[1]]
    
    for i in 2:length(sorted_intervals)
        current_interval = sorted_intervals[i]
        last_merged = merged[end]
        
        # Check for overlap or contiguity (sup >= inf)
        if sup(last_merged) >= (inf(current_interval) - 1e-12)
            # Merge
            new_hi = max(sup(last_merged), sup(current_interval))
            merged[end] = interval(inf(last_merged), new_hi)
        else
            # No overlap, add as a new interval
            push!(merged, current_interval)
        end
    end
    
    return merged
end




"""
    subdivide_tau_interval(tau_rect, num_subdivisions)

Helper function to split a tau interval into a list of N sub-intervals.
"""
function subdivide_tau_interval(tau_rect::Interval{Float64}, num_subdivisions::Int)
    sub_intervals = Interval{Float64}[]
 
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



# --- Main Execution ---

# 1. Update Constants for New Region
const p_start = 199954 // 100000
const p_end = 2 // 1
const p_step = 1 // 100000

const sigma_lower_bound = 172 // 100
const sigma_step = 1 // 1000000


const tp_guess = interval(0, 36//100) # Initial guess for t_p
const target_precision = 1 // 10^9 # Iteration precision
const derivative_threshold = - 1 // 10^9 # Inequality: (l^1)'_p < -1e-9
const tau_subdivisions = 40 #  Number of subdivisions for tau
const tau_subdivision_threshold = 1 // 10^4 # Absolute width threshold


println("Calculating derivative bounds (l^(1))'_p over new grid.")
println("  p range: [$p_start, $p_end] (descending)")
println("  σ range: [$sigma_lower_bound, σ_p] (descending)")
println("  Step size: $p_step")
println("  Inequality to check: (l^(1))'_p < $derivative_threshold")
println("="^40)

start_time = time()
failed_p_intervals = Interval{Float64}[]

# 2. Determine p range (descending)
p_lo_val = min(p_start, p_end)
p_hi_val = max(p_start, p_end)

# Using a manual loop construction to avoid potential memory issues with 'collect' on very small steps
# but keeping the logic consistent with previous versions.
p_values_asc = collect(p_lo_val:p_step:p_hi_val)
if last(p_values_asc) < p_hi_val
    push!(p_values_asc, p_hi_val)
end
p_values_desc = reverse(p_values_asc)

for i in 1:(length(p_values_desc) - 1)
    p_hi = interval(p_values_desc[i])
    p_lo = interval(p_values_desc[i+1])
    p_interval = interval(p_lo, p_hi)

    # 1. Find t_p
    tp_result = find_tp_range(p_interval, tp_guess, target_precision, false)
    if isempty_interval(tp_result)
        println("!!! FAILED to find t_p for P-Interval: ", p_interval)
        push!(failed_p_intervals, p_interval)
        continue 
    end
    t_up_val = sup(tp_result)
    
    # 2. Determine Sigma Range for this P
    # σ_p = (2^p - 1)^(1/p)
    sigma_p_interval = (interval(2)^p_interval - interval(1))^(interval(1)/p_interval)
    sigma_start = sup(sigma_p_interval)
    
    # Construct Sigma Grid
    # Using explicit loop for small step handling
    sigma_values_desc = Float64[]
    let 
        curr = sigma_start
        while curr > sigma_lower_bound + 1e-12
            push!(sigma_values_desc, curr)
            curr -= sigma_step
        end
        push!(sigma_values_desc, sigma_lower_bound)
    end
    
    println("\nP-Interval: ", p_interval)
    println("  σ_p bound for this p: ", sigma_p_interval)
    println("-"^100)
    println("  P-Rect             S-Rect             τ-sub-interval     (l^1)'_p Bound")
    println("-"^100)

    local p_interval_failed = false

    # 3. Sigma Loop
    for j in 1:(length(sigma_values_desc) - 1)
        s_hi = sigma_values_desc[j]
        s_lo = sigma_values_desc[j+1]
        s_interval = interval(s_lo, s_hi)
        
        # 4. Find Tau
        tau_interval_rect = find_new_tau_range(p_interval, s_interval, t_up_val, target_precision, false)
        
        if isempty_interval(tau_interval_rect)
            println("!!! FAIL: Empty Tau for P=$p_interval, S=$s_interval")
            p_interval_failed = true
            break 
        end
        
        # 5. Check Inequality
        deriv_res = calculate_l1_derivative_interval(p_interval, s_interval, tau_interval_rect, tp_result)
        
        println("  P=$(p_interval) S=$(s_interval) τ=$(tau_interval_rect) => l'=$(deriv_res)")
        
        # Check if upper bound is less than threshold (strictly negative check)
        if sup(deriv_res) < derivative_threshold
            continue
        else
            # 6. Subdivision Fallback
            should_subdivide = (sup(tau_interval_rect) - inf(tau_interval_rect)) > tau_subdivision_threshold
            
            if should_subdivide
                println("--- Check failed on full tau interval. Subdividing... ---")
                tau_subs = subdivide_tau_interval(tau_interval_rect, tau_subdivisions)
                all_subs_passed = true
                
                for (k, sub_tau) in enumerate(tau_subs)
                    deriv_sub = calculate_l1_derivative_interval(p_interval, s_interval, sub_tau, tp_result)
                    
                    println("      (Sub) P=$p_interval S=$s_interval τ[$k]=$sub_tau => l'=$deriv_sub")
                    
                    if !(sup(deriv_sub) < derivative_threshold)
                        println("      !!! SUB-FAIL: P=$p_interval, S=$s_interval, τ_sub=$sub_tau")
                        println("          l' = $deriv_sub (Wanted hi < $derivative_threshold)")
                        all_subs_passed = false
                        break
                    end
                end
                
                if all_subs_passed
                    continue 
                end
            end
            
            println("!!! FAIL: P=$p_interval, S=$s_interval, τ=$tau_interval_rect")
            println("    l' = $deriv_res (Wanted hi < $derivative_threshold)")
            p_interval_failed = true
            break
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
    println("All p-intervals verified (l^1)'_p < -1e-9 successfully.")
else
    merged_failed = merge_intervals(failed_p_intervals)
    println("The following p-ranges failed the check:")
    for p_int in merged_failed
        println("  ", p_int)
    end
end
println("Total time: ", round(end_time - start_time, digits=4), " seconds.")
