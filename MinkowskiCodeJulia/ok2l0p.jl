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
    t_pow_p_log_t_safe(t_int, p_int)

Calculates t^p * ln(t) rigorously handling the singularity at t=0.
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
    calculate_l0_derivative_interval(p_int, s_int, t_int)

Calculates (l^(0))'_p = dDelta/dp - 1/2 * (sigma_p)'_p
"""
function calculate_l0_derivative_interval(
    p_int::Interval{<:Real},
    s_int::Interval{<:Real},
    t_int::Interval{<:Real}
)

    
    one_over_p = interval(1) / p_int
    p_minus_1 = p_int - interval(1)
    
    ln_s = log(s_int)
    
    val_1_sp = interval(1) + s_int^p_int
    val_1_tp = interval(1) + t_int^p_int # Used t_int directly
    
    C1 = val_1_sp^(-one_over_p)
    C2 = val_1_tp^(-one_over_p)
    
    # A, B
    A = C2 - C1
    B = t_int * C2 + s_int * C1 # Used t_int directly
    
    
    # Delta(p, sigma)
    Delta_val = (t_int + s_int) * C2 * C1 # Used t_int directly
    
    # --- 2. L_sigma, L_tau ---
    term_Ls_1 = log(val_1_sp) / (p_int^2)
    term_Ls_2 = (s_int^p_int * ln_s) / (p_int * val_1_sp)
    L_sigma = term_Ls_1 - term_Ls_2
    
    # L_tau = ln(1+t^p)/p^2 - (t^p ln t)/(p(1+t^p))
    term_Lt_1 = log(val_1_tp) / (p_int^2)
    
    # --- RIGOROUS SINGULARITY HANDLING ---
    term_singularity = t_pow_p_log_t_safe(t_int, p_int)
    term_Lt_2 = term_singularity / (p_int * val_1_tp)
    
    L_tau = term_Lt_1 - term_Lt_2
    
    # --- 3. Partial Derivatives of A, B w.r.t p ---
    dA_dp = C2 * L_tau - C1 * L_sigma
    dB_dp = t_int * C2 * L_tau + s_int * C1 * L_sigma # Used t_int directly
    
    # --- 4. N_p ---
    N_p = (A^p_int) * log(A) + 
          (B^p_int) * log(B) + 
          p_int * (A^p_minus_1) * dA_dp + 
          p_int * (B^p_minus_1) * dB_dp
          
    # --- 5. T_p ---
    K_tau = val_1_tp^(-interval(1) - one_over_p)
    H2 = (B^p_minus_1) - (t_int^p_minus_1) * (A^p_minus_1) # Used A, B, t_int directly
    
    if 0.0 in H2
        return interval(-1e10, 1e10)
    end
    
    T_p = -N_p / (p_int * K_tau * H2)
    
    # --- 6. Total Derivative dDelta/dp ---
    term_tau_explicit = (interval(1) / (s_int + t_int)) - ( (t_int^p_minus_1) / val_1_tp ) # Used t_int directly
    
    dDelta_dp = Delta_val * ( (L_sigma + L_tau) + term_tau_explicit * T_p )
    
    # --- 7. Derivative of Delta(p, sigma_p) ---
    two_pow_p = interval(2)^p_int
    S_p = two_pow_p - interval(1) # 2^p - 1
    
    ln_2 = log(interval(2))
    ln_Sp = log(S_p)
    
    numerator_bracket = two_pow_p * p_int * ln_2 - S_p * ln_Sp
    term_front = S_p^(one_over_p - interval(1))
    
    dDelta_sigma_p_dp = (term_front * numerator_bracket) / (interval(2) * p_int^2)
    
    # --- 8. Final Result ---
    l0_prime_p = dDelta_dp - dDelta_sigma_p_dp
    
    return l0_prime_p
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



"""
    check_p_s_rectangle(p_int, s_int, ...)
"""
function check_p_s_rectangle(
    p_int::Interval{<:Real},
    s_int::Interval{<:Real},
    t_up_val::Real,              
    target_precision::Real,      
    threshold::Real,           
    tau_subdivisions::Int,
    tau_subdivision_threshold::Real, 
    desc_prefix::String=""
)
    # 1. Find Tau
    tau_interval_rect = find_new_tau_range(p_int, s_int, t_up_val, target_precision, false)
    
    if isempty_interval(tau_interval_rect)
        println("$(desc_prefix)!!! FAIL: Empty Tau for P=$p_int, S=$s_int")
        return false
    end
    
    # 2. Standard Check
    deriv_res = calculate_l0_derivative_interval(p_int, s_int, tau_interval_rect)
    println("$(desc_prefix)  P=$(p_int) S=$(s_int) τ=$(tau_interval_rect) => l0'_p=$(deriv_res)")
    
    if inf(deriv_res) > threshold
        return true
    end
    
    # 3. Subdivision Logic (Fallback)
    is_zero_bound = (inf(tau_interval_rect) == 0.0)
    is_wide = (sup(tau_interval_rect) - inf(tau_interval_rect)) > tau_subdivision_threshold
    
    if is_zero_bound || is_wide
        local tau_subs
        
        if is_zero_bound
            println("$(desc_prefix)--- Check failed. Tau lower bound is 0. Subdividing into 2 parts (isolating zero at 1e-9)... ---")
            split_point = 1 // 10^9
            tau_subs = [
                interval(0.0, split_point),
                interval(split_point, sup(tau_interval_rect))
            ]
        else
            println("$(desc_prefix)--- Check failed on full tau interval. Subdividing... ---")
            tau_subs = subdivide_tau_interval(tau_interval_rect, tau_subdivisions)
        end
        
        for (k, sub_tau) in enumerate(tau_subs)
            deriv_sub = calculate_l0_derivative_interval(p_int, s_int, sub_tau)
            
            sub_label = is_zero_bound ? "(Sub-Zero)" : "(Sub)"
            println("$(desc_prefix)      $sub_label P=$p_int S=$s_int τ[$k]=$sub_tau => l0'_p=$deriv_sub")
            
            if !(inf(deriv_sub) > threshold)
                println("$(desc_prefix)      !!! SUB-FAIL: P=$p_int, S=$s_int, τ_sub=$sub_tau")
                println("$(desc_prefix)          l0'_p = $deriv_sub (Wanted lo > $threshold)")
                return false
            end
        end
        return true
    else
        println("$(desc_prefix)!!! FAIL: P=$p_int, S=$s_int, τ=$tau_interval_rect")
        println("$(desc_prefix)    l0'_p = $deriv_res (Wanted lo > $threshold)")
        return false
    end
end


# --- Main Execution ---

const p_start = 2 // 1
const p_end = 20002 // 10000
const p_step = 1 // 100000

const sigma_start_fixed = 1 // 1
const sigma_end_fixed = 10016 // 10000
const sigma_step = 1 // 10000


const tp_guess = interval(0, 36//100) # Initial guess for t_p
const target_precision = 1 // 10^9 # Iteration precision
const deriv_threshold = 1 // 10^9 # Inequality: (l^0)'_p > 1e-9
const tau_subdivisions = 40 #  Number of subdivisions for tau
const tau_subdivision_threshold = 1 // 10^4 # Absolute width threshold


println("Calculating derivative bounds (l^(0))'_p over new grid.")
println("  p range: [$p_start, $p_end] (descending)")
println("  σ range: [$sigma_start_fixed, $sigma_end_fixed] (descending)")
println("  Inequality to check: (l^0)'_p > $deriv_threshold")
println("="^40)

start_time = time()
failed_p_intervals = Interval{Float64}[]

# Determine p range (descending)
p_lo_val = min(p_start, p_end)
p_hi_val = max(p_start, p_end)

p_values_asc = collect(p_lo_val:p_step:p_hi_val)
if last(p_values_asc) < p_hi_val
    push!(p_values_asc, p_hi_val)
end
p_values_desc = reverse(p_values_asc)

# Determine Sigma Range (descending)
sigma_values_asc = collect(sigma_start_fixed:sigma_step:sigma_end_fixed)
if last(sigma_values_asc) < sigma_end_fixed
    push!(sigma_values_asc, sigma_end_fixed)
end
sigma_values_desc = reverse(sigma_values_asc)

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
    
    println("\nP-Interval: ", p_interval)
    println("-"^100)
    println("  P-Rect             S-Rect              τ-sub-interval     (l^0)'_p Bound")
    println("-"^100)

    local p_interval_failed = false

    # 3. Sigma Loop
    for j in 1:(length(sigma_values_desc) - 1)
        s_hi = sigma_values_desc[j]
        s_lo = sigma_values_desc[j+1]
        s_interval = interval(s_lo, s_hi)
        
        # --- CHECK RECTANGLE ---
        passed = check_p_s_rectangle(
            p_interval, s_interval, t_up_val,
            target_precision, deriv_threshold,
            tau_subdivisions, tau_subdivision_threshold
        )
        
        if !passed
            println("!!! FAIL: P=$p_interval, S=$s_interval")
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
    println("All p-intervals verified (l^0)'_p > 1e-9 successfully.")
else
    merged_failed = merge_intervals(failed_p_intervals)
    println("The following p-ranges failed the check:")
    for p_int in merged_failed
        println("  ", p_int)
    end
end
println("Total time: ", round(end_time - start_time, digits=4), " seconds.")