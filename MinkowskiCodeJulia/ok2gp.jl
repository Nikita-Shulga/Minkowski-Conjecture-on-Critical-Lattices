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
    calculate_gp_interval(p_int, s_int, t_int)

Calculates the interval for the total derivative dg/dp.
"""
function calculate_gp_interval(
    p_int::Interval{<:Real},
    s_int::Interval{<:Real},
    t_int::Interval{<:Real},
)
    # Ensure t_int is positive for log() operations

    # --- 1. Basic Terms ---
    ln_s = log(s_int)
    ln_t = log(t_int)
    
    val_1_sp = interval(1) + s_int^p_int
    val_1_tp = interval(1) + t_int^p_int
    
    C1 = val_1_sp^(-interval(1)/p_int)
    C2 = val_1_tp^(-interval(1)/p_int)
    
    # --- 2. p-Partial Derivatives of Basic Terms ---
    # C1_p
    term_c1_a = log(val_1_sp) / (p_int^2)
    term_c1_b = (s_int^p_int * ln_s) / (p_int * val_1_sp)
    C1_p = C1 * (term_c1_a - term_c1_b)

    # C2_p
    term_c2_a = log(val_1_tp) / (p_int^2)
    term_c2_b = (t_int^p_int * ln_t) / (p_int * val_1_tp)
    C2_p = C2 * (term_c2_a - term_c2_b)
    
    # D1, D2 and their p-partials
    t_pow_pm1 = t_int^(p_int - interval(1))
    s_pow_pm1 = s_int^(p_int - interval(1))
    
    D1 = interval(1) - s_int * t_pow_pm1
    D1_p = -s_int * t_pow_pm1 * ln_t
    
    D2 = interval(1) - t_int * s_pow_pm1
    D2_p = -t_int * s_pow_pm1 * ln_s
    
    # A, B and their p-partials
    A = C2 - C1
    A_p = C2_p - C1_p
    
    B = t_int * C2 + s_int * C1
    B_p = t_int * C2_p + s_int * C1_p
    
    # H1, H2 and their p-partials
    # Term 1: B^(p-1)
    B_pow_pm1 = B^(p_int - interval(1)) 
    B_term_p = B_pow_pm1 * log(B) + (p_int - interval(1)) * (B^(p_int - interval(2))) * B_p
    
    # Term 2: (sA)^(p-1)
    sA = s_int * A
    sA_pow_pm1 = sA^(p_int - interval(1))
    sA_term_p = sA_pow_pm1 * log(sA) + (p_int - interval(1)) * (sA^(p_int - interval(2))) * (s_int * A_p)
    
    H1 = B_pow_pm1 + sA_pow_pm1
    H1_p = B_term_p + sA_term_p
    
    # Term 3: (tA)^(p-1)
    tA = t_int * A
    tA_pow_pm1 = tA^(p_int - interval(1))
    tA_term_p = tA_pow_pm1 * log(tA) + (p_int - interval(1)) * (tA^(p_int - interval(2))) * (t_int * A_p)
    
    H2 = B_pow_pm1 - tA_pow_pm1
    H2_p = B_term_p - tA_term_p
    
    # --- 3. Compute k vectors ---
    k1 = C1_p * D1 * H1 - C2_p * D2 * H2
    k2 = C1 * D1_p * H1 - C2 * D2_p * H2
    k3 = C1 * D1 * H1_p - C2 * D2 * H2_p
    
    partial_g_p = k1 + k2 + k3
    
    # --- 4. Constraint numerator F_p ---
      N_p = A^p_int * log(A) + B^p_int * log(B) + p_int * A^(p_int - interval(1)) * A_p + p_int * B^(p_int - interval(1)) * B_p
    
    # --- 5. g_tau Calculation ---
    K_tau = val_1_tp^(-interval(1) - interval(1)/p_int)
    
    # Derivatives wrt tau
    C1_tau = 0.0
    C2_tau = - (t_int^(p_int - interval(1))) * K_tau
    
    D1_tau = -(p_int - interval(1)) * s_int * (t_int^(p_int - interval(2)))
    D2_tau = -s_pow_pm1
    
    dA_dtau = - (t_int^(p_int - interval(1))) * K_tau
    dB_dtau = K_tau
    
    # Explicit partial names
    dH1_dtau = (p_int - interval(1)) * ( (B^(p_int - interval(2))) * dB_dtau + s_pow_pm1 * (A^(p_int - interval(2))) * dA_dtau )
    dH2_dtau = (p_int - interval(1)) * ( (B^(p_int - interval(2))) * dB_dtau - (t_int^(p_int - interval(2))) * (A^(p_int - interval(1))) - (t_int^(p_int - interval(1))) * (A^(p_int - interval(2))) * dA_dtau )
    
    g_tau = C1 * D1_tau * H1 + 
            C1 * D1 * dH1_dtau - 
            C2_tau * D2 * H2 - 
            C2 * D2_tau * H2 - 
            C2 * D2 * dH2_dtau
            
    # --- 6. Final Total Derivative ---
    if in_interval(0,H2)== true 
        println("  Warning: H_2 interval contains zero.")
        return emptyinterval(Float64)
    end
    
    tau_p = -N_p / (p_int * K_tau * H2)
    
    total_gp = partial_g_p + tau_p * g_tau
    
    return total_gp
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

const p_start = 199 // 100
const p_end = 20076 // 10000
const p_step = 1 // 10000

const sigma_fixed_start = 102 // 100
const sigma_fixed_end = 172 // 100
const sigma_step = 1 // 10000


const tp_guess = interval(0, 36//100) # Initial guess for t_p
const target_precision = 1 // 10^9 # Iteration precision
const gp_threshold = 1 // 10^9 # Inequality: gp > 1e-9
const tau_subdivisions = 10 #  Number of subdivisions for tau
const tau_subdivision_threshold = 1 // 10^4 # Absolute width threshold

println("Calculating bounds for g'_p over (p, σ) grid.")
println("  p range: [$p_start, $p_end] (descending)")
println("  σ range: [$sigma_fixed_start, $sigma_fixed_end] (fixed rectangle)")
println("  Inequality to check: g'_p > $gp_threshold")
println("="^40)

start_time = time()
failed_p_intervals = Interval{Float64}[]

p_lo_val = min(p_start, p_end)
p_hi_val = max(p_start, p_end)

p_values_asc = collect(p_lo_val:p_step:p_hi_val)
if last(p_values_asc) < p_hi_val
    push!(p_values_asc, p_hi_val)
end
p_values_desc = reverse(p_values_asc)

sigma_values_asc = collect(sigma_fixed_start:sigma_step:sigma_fixed_end)
if last(sigma_values_asc) < sigma_fixed_end
    push!(sigma_values_asc, sigma_fixed_end)
end
sigma_values_desc = reverse(sigma_values_asc)


for i in 1:(length(p_values_desc) - 1)
    p_hi = interval(p_values_desc[i])
    p_lo = interval(p_values_desc[i+1])
    p_interval = interval(p_lo, p_hi)

    # 1. Find t_p (needed for otau_0)
    tp_result = find_tp_range(p_interval, tp_guess, target_precision, false)
    if isempty_interval(tp_result)
        println("!!! FAILED to find t_p for P-Interval: ", p_interval)
        push!(failed_p_intervals, p_interval)
        continue 
    end
    
    t_up_val = sup(tp_result)
    
    println("\nP-Interval: ", p_interval)
    println("-"^100)
    println("  P-Rect             S-Rect             τ-sub-interval     g'_p Bound")
    println("-"^100)

    local p_interval_failed = false

    # 2. Sigma Loop
    for j in 1:(length(sigma_values_desc) - 1)
        s_hi = sigma_values_desc[j]
        s_lo = sigma_values_desc[j+1]
        s_interval = interval(s_lo, s_hi)
        
        # 3. Find Tau
        tau_interval_rect = find_new_tau_range(p_interval, s_interval, t_up_val, target_precision, false)
        
        if isempty_interval(tau_interval_rect)
            println("!!! FAIL: Empty Tau for P=$p_interval, S=$s_interval")
            p_interval_failed = true
            break 
        end
        
        # 4. Check g'_p inequality
        gp_res = calculate_gp_interval(p_interval, s_interval, tau_interval_rect)
        
        # Catch the singularity from step 1
        if isempty_interval(gp_res)
            println("!!! FAIL: Singularity detected in H2 for P=$p_interval, S=$s_interval")
            p_interval_failed = true
            break 
        end
        
        # PRINT RESULT BEFORE CHECKING
        println("  P=$(p_interval) S=$(s_interval) τ=$(tau_interval_rect) => g'_p=$(gp_res)")
        
        if inf(gp_res) > gp_threshold
            continue
        else
            # 5. Subdivision Fallback
            should_subdivide = true
             
            if should_subdivide
                println("--- Check failed on full tau interval. Subdividing... ---")
                tau_subs = subdivide_tau_interval(tau_interval_rect, tau_subdivisions)
                all_subs_passed = true
                
                for (k, sub_tau) in enumerate(tau_subs)
                    gp_sub = calculate_gp_interval(p_interval, s_interval, sub_tau)
                    
                    #  Catch singularity in subdivision
                    if isempty_interval(gp_sub)
                        println("      !!! SUB-FAIL: Singularity in H2 at P=$p_interval, S=$s_interval, τ_sub=$sub_tau")
                        all_subs_passed = false
                        break
                    end
                    
                    println("      (Sub) P=$p_interval S=$s_interval τ[$k]=$sub_tau => g'_p=$gp_sub")
                    
                    if !(inf(gp_sub) > gp_threshold)
                        println("      !!! SUB-FAIL: P=$p_interval, S=$s_interval, τ_sub=$sub_tau")
                        println("          g'_p = $gp_sub (Wanted lo > $gp_threshold)")
                        all_subs_passed = false
                        break
                    end
                end
                
                if all_subs_passed
                     continue 
                end
            end
            
            println("!!! FAIL: P=$p_interval, S=$s_interval, τ=$tau_interval_rect")
            println("    g'_p = $gp_res (Wanted lo > $gp_threshold)")
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
    println("All p-intervals verified g'_p > 1e-9 successfully.")
else
    merged_failed = merge_intervals(failed_p_intervals)
    println("The following p-ranges failed the check:")
    for p_int in merged_failed
        println("  ", p_int)
    end
end
println("Total time: ", round(end_time - start_time, digits=4), " seconds.")