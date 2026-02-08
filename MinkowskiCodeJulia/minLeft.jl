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
    calculate_g_interval(p_int, s_int, t_int)

Calculates the interval for the function g(p,σ) given intervals for p, σ, and τ.
"""
function calculate_g_interval(
    p_int::Interval{<:Real},
    s_int::Interval{<:Real},
    t_int::Interval{<:Real}
)
    # Helper variables
    one_over_p = interval(1) / p_int
    p_minus_1 = p_int - interval(1)
    
    # C_1 and C_2 (used in A, B, and g)
    t_term_pow = (interval(1) + t_int^p_int)
    s_term_pow = (interval(1) + s_int^p_int)
    C_2 = t_term_pow^(-one_over_p) 
    C_1 = s_term_pow^(-one_over_p) 
    
    # A = C_2 - C_1
    A = C_2 - C_1
    
    # B = τ*C_2 + σ*C_1
    B = t_int * C_2 + s_int * C_1
    
    # Powers
    A_pow = A^p_minus_1
    B_pow = B^p_minus_1
    t_pow = t_int^p_minus_1
    s_pow = s_int^p_minus_1

    # D_1 and D_2
    D_1 = interval(1) - s_int * t_pow
    D_2 = interval(1) - t_int * s_pow

    # H_1 and H_2
    H_1 = B_pow + s_pow * A_pow
    H_2 = B_pow - t_pow * A_pow

    # g = C_1 * D_1 * H_1 - C_2 * D_2 * H_2
    G_left = C_1 * D_1 * H_1
    G_right = C_2 * D_2 * H_2
    
    g_interval = G_left - G_right
    
    return g_interval
end

"""
    calculate_h_interval(p_int, s_int, t_int)

Calculates the interval for the function h(p,σ) = Σ h_i.
"""


function calculate_h_interval(
    p_int::Interval{<:Real},
    s_int::Interval{<:Real},
    t_int::Interval{<:Real}
)
    # --- 1. Helper Variables ---
    one_over_p = interval(1) / p_int
    p_minus_1 = p_int - interval(1)
    p_minus_2 = p_int - interval(2)

    # --- 2. Component Functions (C, A, B, D, H) ---
    t_term_pow = (interval(1) + t_int^p_int)
    s_term_pow = (interval(1) + s_int^p_int)
    C_2 = t_term_pow^(-one_over_p) 
    C_1 = s_term_pow^(-one_over_p) 
    
    A = C_2 - C_1
    B = t_int * C_2 + s_int * C_1
    
    t_pow_pm1 = t_int^p_minus_1
    s_pow_pm1 = s_int^p_minus_1
    t_pow_pm2 = t_int^p_minus_2
    s_pow_pm2 = s_int^p_minus_2

    D_1 = interval(1) - s_int * t_pow_pm1
    D_2 = interval(1) - t_int * s_pow_pm1
    
    A_pow_pm1 = A^p_minus_1
    A_pow_pm2 = A^p_minus_2
    B_pow_pm1 = B^p_minus_1
    B_pow_pm2 = B^p_minus_2

    H_1 = B_pow_pm1 + s_pow_pm1 * A_pow_pm1
    H_2 = B_pow_pm1 - t_pow_pm1 * A_pow_pm1

    # --- 3. Derivative Shorthands (K) ---
    K_sigma = s_term_pow^(-interval(1) - one_over_p)
    K_tau   = t_term_pow^(-interval(1) - one_over_p)

    # --- 4. Component Partial Derivatives (A, B) ---
    dA_dsigma = s_pow_pm1 * K_sigma
    dA_dtau   = -t_pow_pm1 * K_tau
    dB_dsigma = K_sigma
    dB_dtau   = K_tau

    # --- 5. H-Partial Derivatives ---
    dH1_dsigma = (p_minus_1) * ( B_pow_pm2 * dB_dsigma + s_pow_pm2 * A_pow_pm1 + s_pow_pm1 * A_pow_pm2 * dA_dsigma )
    dH2_dsigma = (p_minus_1) * ( B_pow_pm2 * dB_dsigma - t_pow_pm1 * A_pow_pm2 * dA_dsigma )
    dH1_dtau   = (p_minus_1) * ( B_pow_pm2 * dB_dtau + s_pow_pm1 * A_pow_pm2 * dA_dtau )
    dH2_dtau   = (p_minus_1) * ( B_pow_pm2 * dB_dtau - t_pow_pm2 * A_pow_pm1 - t_pow_pm1 * A_pow_pm2 * dA_dtau )

    # --- 6. Tau Derivative (T) ---
    # We must check for H_2 being zero
    if in_interval(0,H_2)== true 
        println("  Warning: H_2 interval contains zero. Skipping h calculation.")
        return emptyinterval(Float64)
    end
    T = - (K_sigma / K_tau) * (H_1 / H_2)

    # --- 7. h_i Components ---
    h_1 = -s_pow_pm1 * K_sigma * D_1 * H_1 - t_pow_pm1 * K_sigma * D_2 * H_1
    h_2 = -t_pow_pm1 * C_1 * H_1
    h_3 = (p_minus_1) * t_int * s_pow_pm2 * C_2 * H_2
    h_4 = -(p_minus_1) * s_int * t_pow_pm2 * C_1 * H_1 * T
    h_5 = s_pow_pm1 * C_2 * H_2 * T
    h_6 = C_1 * D_1 * dH1_dsigma - C_2 * D_2 * dH2_dsigma
    h_7 = ( C_1 * D_1 * dH1_dtau - C_2 * D_2 * dH2_dtau ) * T

    # --- 8. Final Sum ---
    h_interval = h_1 + h_2 + h_3 + h_4 + h_5 + h_6 + h_7
    
    return h_interval
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

# 1. Define constants
const p_start = 198 // 100
const p_end = 199 // 100
const p_step = 5 // 100000
const sigma_step = 1 // 100000


const tp_guess = interval(0, 36//100) # Initial guess for t_p
const target_precision = 1 // 10^9 # Iteration precision
const g_h_precision = 1 // 10^9 # Comparison precision
const tau_subdivisions = 10 #  Number of subdivisions for tau
const tau_subdivision_threshold = 1 // 10^4 # Absolute width threshold

println("Calculating bounds for τ over (p, σ) grid.")
println("  p range: [$p_start, $p_end] (step $p_step), descending")
println("  σ range: [1, σ_p] (step $sigma_step), descending")
println("  Target precision: $target_precision")
println("  G/H Check precision: $g_h_precision")
println("  Tau subdivisions: $tau_subdivisions")
println("  Tau subdivision threshold: $tau_subdivision_threshold")
println("="^40)

# Record start time
start_time = time()

failed_p_intervals = Interval{Float64}[]

# --- Outer loop over p ---

# Determine the true min and max of the p-range
p_lo_val = min(p_start, p_end)
p_hi_val = max(p_start, p_end)

# Build the ascending range from the true min to the true max
p_values_asc = collect(p_lo_val:p_step:p_hi_val)
if last(p_values_asc) < p_hi_val
    push!(p_values_asc, p_hi_val) # Ensure we cover the full range
end

# Reverse the p_values array to iterate in descending order
p_values_desc = reverse(p_values_asc)

for i in 1:(length(p_values_desc) - 1)
    p_hi = interval(p_values_desc[i])   # Higher value
    p_lo = interval(p_values_desc[i+1]) # Lower value
    p_interval = interval(p_lo, p_hi)
    
    
    # --- 2. Calculate bounds for t_p, σ_p, and endpoint Deltas ---
    
    # 2a. Find t_p
    tp_result = find_tp_range(p_interval, tp_guess, target_precision, false)
    if isempty_interval(tp_result)
        println("!!! FAILED to find t_p for P-Interval: ", p_interval)
        push!(failed_p_intervals, p_interval)
        continue # Skip this p-interval
    end
    
    # otau_0 = t_{\up} = t_p(p_lo)
    t_up_val = sup(tp_result)
    
    # 2b. Find σ_p = (2^p - 1)^(1/p)
    sigma_p_lo = (interval(2)^p_lo - interval(1))^(interval(1) / p_lo)
    sigma_p_hi = (interval(2)^p_hi - interval(1))^(interval(1) / p_hi)
    sigma_p_interval = interval(sigma_p_lo, sigma_p_hi)
    
    # 2c. Calculate Full Intervals for Endpoint Deltas
    # Δ(p,1) = 4^(-1/p) * (1+t_p)/(1-t_p)
    delta_p1_interval = (interval(4)^(-interval(1) / p_interval)) * ( (interval(1) + tp_result) / (interval(1) - tp_result) )
    
    # Δ(p,σ_p) = (1/2) * σ_p
    delta_psigmap_interval = interval(1//2) * sigma_p_interval

    # 2d. Find the minimum of the endpoint deltas
    min_delta_endpoint = min(delta_p1_interval, delta_psigmap_interval)

    # --- 3. Setup Sigma Loop (Descending) ---
    sigma_start = 1//1
    sigma_end = sup(sigma_p_interval) # Loop up to the upper bound
    
    sigma_values_asc = collect(sigma_start:sigma_step:sigma_end)
    if last(sigma_values_asc) < sigma_end
        push!(sigma_values_asc, sigma_end)
    end
    sigma_values_desc = reverse(sigma_values_asc)
    
    # --- 4. Print intermediate results ---
    println("\nP-Interval: ", p_interval)
    println("    Found tp bound: ", tp_result)
    println("    Using otau_0 = t_{up} = ", t_up_val) # Print the initial tau bound
    println("    Found sigma_p bound: ", sigma_p_interval)
    println("    Found Delta(p,1) bound: ", delta_p1_interval)
    println("    Found Delta(p,sigma_p) bound: ", delta_psigmap_interval)
    println("    Min Delta Endpoint: ", min_delta_endpoint)
    println("    Sigma range: ", interval(sigma_start, sigma_end), " with ", length(sigma_values_desc)-1, " descending steps")
    println("-"^40)
    println("  P-Rect             S-Rect             τ-sub-interval     Delta/G/H Bound")
    println("-"^40)

    local first_delta_fail_printed = false
    local first_g_fail_printed = false
    local first_h_fail_printed = false
    local is_checking_g = false # State flag: Delta has failed
    local is_checking_h = false # State flag: G has failed

    # --- NEW: Flag to track if this specific p-interval has failed rigorously ---
    local p_interval_failed = false

    # --- 5. Inner loop over sigma (descending) ---
    for j in 1:(length(sigma_values_desc) - 1)
        s_hi = sigma_values_desc[j]   # Higher value
        s_lo = sigma_values_desc[j+1] # Lower value
        s_interval = interval(s_lo, s_hi)
        
        # 5a. Find the new τ bound for this rectangle
        tau_interval_rect = find_new_tau_range(p_interval, s_interval, t_up_val, target_precision, false) # Set verbose=false
        
        # 5b. Handle empty tau interval
        if isempty_interval(tau_interval_rect)
            if !first_delta_fail_printed; println("!!! FIRST FAIL for Delta: P=$(p_interval), S=$(s_interval) (Resulting tau was ∅)"); first_delta_fail_printed = true; end
            if is_checking_g && !first_g_fail_printed; println("!!! FIRST FAIL for G: P=$(p_interval), S=$(s_interval) (Resulting tau was ∅)"); first_g_fail_printed = true; end
            if is_checking_h && !first_h_fail_printed; println("!!! FIRST FAIL for H: P=$(p_interval), S=$(s_interval) (Resulting tau was ∅)"); first_h_fail_printed = true; end
            
            # Critical failure
            p_interval_failed = true
            break # Exit inner sigma loop
        end
        
        # 5c. Check if subdivision is needed *if* a check fails
        should_subdivide_if_fails = (inf(tau_interval_rect) * 2 < sup(tau_interval_rect)) ||
                                    ((sup(tau_interval_rect) - inf(tau_interval_rect)) > tau_subdivision_threshold)
        
        local rectangle_failed_at_some_tau_sub_interval = false

        # --- 5d. STAGE 1: DELTA CHECK ---
        if !is_checking_g
            delta_ps_interval = begin
                one_over_p = interval(1) / p_interval
                term1 = tau_interval_rect + s_interval
                term2 = (interval(1) + tau_interval_rect^p_interval)^(-one_over_p)
                term3 = (interval(1) + s_interval^p_interval)^(-one_over_p)
                term1 * term2 * term3
            end
            
            println("  P=$(p_interval) S=$(s_interval) τ[1]=$(tau_interval_rect) => Δ=$(delta_ps_interval)")
            
            inequality_holds = inf(delta_ps_interval) > sup(min_delta_endpoint)
            
            if !inequality_holds # Delta check failed
                if !first_delta_fail_printed; println("!!! FIRST FAIL for Delta: P=$(p_interval), S=$(s_interval), τ[1]=$(tau_interval_rect) (Δ.lo <= min_Δ.hi)"); first_delta_fail_printed = true; end
                
                if should_subdivide_if_fails
                    println("--- Check failed on full tau interval. Subdividing for Delta... ---")
                    tau_intervals_to_check = subdivide_tau_interval(tau_interval_rect, tau_subdivisions)
                    all_sub_intervals_passed = true
                    
                    for (k, current_tau_interval) in enumerate(tau_intervals_to_check)
                        one_over_p = interval(1) / p_interval
                        term1 = current_tau_interval + s_interval
                        term2 = (interval(1) + current_tau_interval^p_interval)^(-one_over_p)
                        term3 = (interval(1) + s_interval^p_interval)^(-one_over_p)
                        delta_ps_interval_sub = term1 * term2 * term3
                        
                        println("      (Sub-Check Δ) P=$(p_interval) S=$(s_interval) τ[$k]=$(current_tau_interval) => Δ=$(delta_ps_interval_sub)")
                        
                        if !(inf(delta_ps_interval_sub) > sup(min_delta_endpoint))
                            println("!!! SUB-FAIL for Delta: P=$(p_interval), S=$(s_interval), τ[$k]=$(current_tau_interval) (Δ.lo <= min_Δ.hi)")
                            all_sub_intervals_passed = false
                            break # Exit sub-interval loop
                        end
                    end
                    
                    if all_sub_intervals_passed
                        # Subdivision worked, this rectangle is fine for Delta.
                        continue # Go to next sigma rectangle
                    else
                        # Subdivision failed, this is a real Delta failure.
                        is_checking_g = true
                    end
                else
                    # Failed on a narrow interval, this is a real Delta failure.
                    is_checking_g = true
                end
            else
                # Delta check passed on the full interval.
                continue # Go to next sigma rectangle
            end
        end
        
        # --- 5e. STAGE 2: G CHECK ---
        if is_checking_g && !is_checking_h
            g_interval = calculate_g_interval(p_interval, s_interval, tau_interval_rect)
            println("      (Checking G) P=$(p_interval) S=$(s_interval) τ[1]=$(tau_interval_rect) => g_bound=$(g_interval)")
            
            g_inequality_holds = sup(g_interval) < -g_h_precision

            if !g_inequality_holds # G check failed
                if !first_g_fail_printed; println("!!! FIRST FAIL for G: P=$(p_interval), S=$(s_interval), τ[1]=$(tau_interval_rect) (sup(g) >= -$g_h_precision)"); first_g_fail_printed = true; end
                
                if should_subdivide_if_fails
                    println("--- G-Check failed on full tau interval. Subdividing for G... ---")
                    tau_intervals_to_check = subdivide_tau_interval(tau_interval_rect, tau_subdivisions)
                    all_sub_intervals_passed = true
                    
                    for (k, current_tau_interval) in enumerate(tau_intervals_to_check)
                        g_interval_sub = calculate_g_interval(p_interval, s_interval, current_tau_interval)
                        println("      (Sub-Check G) P=$(p_interval) S=$(s_interval) τ[$k]=$(current_tau_interval) => g_bound=$(g_interval_sub)")
                        
                        if !(sup(g_interval_sub) < -g_h_precision)
                            println("!!! SUB-FAIL for G: P=$(p_interval), S=$(s_interval), τ[$k]=$(current_tau_interval) (sup(g) >= -$g_h_precision)")
                            all_sub_intervals_passed = false
                            break # Exit sub-interval loop
                        end
                    end
                    
                    if all_sub_intervals_passed
                        continue # Go to next sigma
                    else
                        is_checking_h = true
                    end
                else
                    is_checking_h = true
                end
            else
                continue # G check passed, go to next sigma
            end
        end

        # --- 5f. STAGE 3: H CHECK ---
	if is_checking_h
  		  h_interval = calculate_h_interval(p_interval, s_interval, tau_interval_rect)
  	  
   			 # NEW: If h is empty, we couldn't prove it. Mark as fail.
    	if isempty_interval(h_interval)
       		 println("      !!! RIGOR FAIL: h_bound is empty (singularities).")
        	 p_interval_failed = true
        	break 
   	end

   	println("      (Checking H) P=$(p_interval) S=$(s_interval) τ[1]=$(tau_interval_rect) => h_bound=$(h_interval)")
   

            h_inequality_holds = sup(h_interval) < -g_h_precision

            if !h_inequality_holds # H check failed
                if !first_h_fail_printed; println("!!! FIRST FAIL for H: P=$(p_interval), S=$(s_interval), τ[1]=$(tau_interval_rect) (sup(h) >= -$g_h_precision)"); first_h_fail_printed = true; end
                
                if should_subdivide_if_fails
                    println("--- H-Check failed on full tau interval. Subdividing for H... ---")
                    tau_intervals_to_check = subdivide_tau_interval(tau_interval_rect, tau_subdivisions)
                    all_sub_intervals_passed = true
                    
                    for (k, current_tau_interval) in enumerate(tau_intervals_to_check)
                        h_interval_sub = calculate_h_interval(p_interval, s_interval, current_tau_interval)
                        println("      (Sub-Check H) P=$(p_interval) S=$(s_interval) τ[$k]=$(current_tau_interval) => h_bound=$(h_interval_sub)")
                        
                        if !(sup(h_interval_sub) < -g_h_precision)
                            println("!!! SUB-FAIL for H: P=$(p_interval), S=$(s_interval), τ[$k]=$(current_tau_interval) (sup(h) >= -$g_h_precision)")
                            all_sub_intervals_passed = false
                            break # Exit sub-interval loop
                        end
                    end
                    
                    if all_sub_intervals_passed
                        continue # Go to next sigma
                    else
                        rectangle_failed_at_some_tau_sub_interval = true
                    end
                else
                    rectangle_failed_at_some_tau_sub_interval = true
                end
            else
                continue # H check passed, go to next sigma
            end

            if rectangle_failed_at_some_tau_sub_interval
                p_interval_failed = true
                break # This is a real, final failure. Exit sigma loop.
            end
        end

    end # --- End of inner sigma loop ---

    # --- Record the failure if checking H (and subdivision) wasn't enough ---
    if p_interval_failed
        push!(failed_p_intervals, p_interval)
    end

end # --- End of outer p loop ---

# Calculate and print total time
end_time = time()
println("="^40)

# --- Final Report ---
println("Final Report on Inconclusive P-Intervals:")
if isempty(failed_p_intervals)
    println("All p-intervals were processed successfully.")
else
    # Merge the intervals first
    merged_failed_intervals = merge_intervals(failed_p_intervals)
    
    println("The following p-ranges were inconclusive (failed the h-check):")
    for p_int in merged_failed_intervals
        println("  ", p_int)
    end
end

println("Total time spent: ", round(end_time - start_time, digits=4), " seconds.")