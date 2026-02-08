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




# ==============================================================================
# 2. FORMATTED CHECKING LOGIC
# ==============================================================================

"""
    calculate_and_format_check(phase, p_int, s_int, t_int, min_delta, g_h_tol, is_sub, sub_idx=0)

Performs the calculation, prints the result with P-Phase prefix and format.
"""
function calculate_and_format_check(phase::Int, p_int, s_int, t_int, min_delta, g_h_tol, is_sub::Bool, sub_idx::Int=1)
    
    passed = false
    calc_val = emptyinterval(Float64)
    label = ""
    bound_name = ""
    fail_condition_msg = ""
    
    if phase == 1
        label = "H"
        bound_name = "h_bound"
        calc_val = calculate_h_interval(p_int, s_int, t_int)
        passed = (sup(calc_val) < -g_h_tol)
        fail_condition_msg = "(sup(h) >= -$g_h_tol)"
    elseif phase == 2
        label = "G"
        bound_name = "g_bound"
        calc_val = calculate_g_interval(p_int, s_int, t_int)
        passed = (inf(calc_val) > g_h_tol)
        fail_condition_msg = "(inf(g) <= $g_h_tol)"
    elseif phase == 3
        label = "Delta"
        bound_name = "Δ"
        one_over_p = interval(1) / p_int
        term1 = t_int + s_int
        term2 = (interval(1) + t_int^p_int)^(-one_over_p)
        term3 = (interval(1) + s_int^p_int)^(-one_over_p)
        calc_val = term1 * term2 * term3
        passed = (inf(calc_val) > sup(min_delta))
        fail_condition_msg = "(Δ.lo <= min_Δ.hi)"
    elseif phase == 4
        label = "G"
        bound_name = "g_bound"
        calc_val = calculate_g_interval(p_int, s_int, t_int)
        passed = (sup(calc_val) < -g_h_tol)
        fail_condition_msg = "(sup(g) >= -$g_h_tol)"
    elseif phase == 5
        label = "H"
        bound_name = "h_bound"
        calc_val = calculate_h_interval(p_int, s_int, t_int)
        passed = (sup(calc_val) < -g_h_tol)
        fail_condition_msg = "(sup(h) >= -$g_h_tol)"
    end
    
    # --- FORMATTING ---
    p_tag = "P$phase"
    indent = is_sub ? "      " : " "
    check_type = is_sub ? "(Sub-Check $label)" : "(Checking $label)"
    tau_str = is_sub ? "τ[$sub_idx]=$t_int" : "τ[1]=$t_int"
    
    # The main result line
    println("$indent$p_tag $check_type P=$p_int S=$s_int $tau_str => $bound_name=$calc_val")
    
    # Failure message (if failed)
    if !passed
        fail_prefix = is_sub ? "!!! SUB-FAIL" : "!!! FIRST FAIL"
        println("$fail_prefix for $label: P=$p_int, S=$s_int, $tau_str $fail_condition_msg")
    end
    
    return passed, calc_val
end

# ==============================================================================
# 3. MAIN EXECUTION
# ==============================================================================

const p_start = 258 // 100
const p_end   = 2642 // 1000
const p_step  = 1 // 10000
const sigma_step = 7 // 100000

const tp_guess = interval(0, 36//100) # Initial guess for t_p
const target_precision = 1 // 10^9 # Iteration precision
const g_h_precision = 1 // 10^9 # Comparison precision
const tau_subdivisions = 50 #  Number of subdivisions for tau
const tau_subdivision_threshold = 1 // 10^4 # Absolute width threshold

println("Checking Bell Curve Region ($p_start, $p_end)")
println("Sequence: h < -ε -> g > ε -> Delta -> g < -ε -> h < -ε")
println("Requirement: At least one rectangle MUST pass for EACH phase (1-5) per P-Interval.")

start_time = time()
failed_p_intervals = Interval{Float64}[]

# P Loop
p_lo_val = min(p_start, p_end); p_hi_val = max(p_start, p_end)
p_values_asc = collect(p_lo_val:p_step:p_hi_val)
if last(p_values_asc) < p_hi_val; push!(p_values_asc, p_hi_val); end
p_values_desc = reverse(p_values_asc)

for i in 1:(length(p_values_desc) - 1)
    p_hi = interval(p_values_desc[i]) 
    p_lo = interval(p_values_desc[i+1])
    p_interval = interval(p_lo, p_hi)

    tp_result = find_tp_range(p_interval, tp_guess, target_precision)
    if isempty_interval(tp_result)
        println("!!! FAILED to find t_p for P=$(p_interval)")
        push!(failed_p_intervals, p_interval); continue
    end
    
    t_up_val = sup(tp_result)
    sigma_p_lo = (interval(2)^p_lo - interval(1))^(interval(1) / p_lo)
    sigma_p_hi = (interval(2)^p_hi - interval(1))^(interval(1) / p_hi)
    sigma_p_interval = interval(sigma_p_lo, sigma_p_hi)
    
    delta_p1_interval = (interval(4)^(-interval(1) / p_interval)) * ( (interval(1) + tp_result) / (interval(1) - tp_result) )
    delta_psigmap_interval = interval(1//2) * sigma_p_interval
    min_delta_endpoint = min(delta_p1_interval, delta_psigmap_interval)
	
    sigma_start = 1.0
    sigma_end = sup(sigma_p_interval)
    s_vals = collect(sigma_start:sigma_step:sigma_end)
    if last(s_vals) < sigma_end; push!(s_vals, sigma_end); end
    s_vals_desc = reverse(s_vals)

    # --- HEADER PRINTING ---
    println("\nP-Interval: ", p_interval)
    println("    Found tp bound: ", tp_result)
    println("    Using otau_0 = t_{up} = ", t_up_val)
    println("    Found sigma_p bound: ", sigma_p_interval)
    println("    Min Delta Endpoint: ", min_delta_endpoint)
    println("    Sigma range: ", interval(1.0, sigma_end), " with ", length(s_vals_desc)-1, " descending steps")
    println("-"^40)
    println("  P-Rect             S-Rect             τ-sub-interval     Delta/G/H Bound")
    println("-"^40)

    local p_interval_failed = false
    local current_phase = 1 
    
    # --- NEW: Track which phases have been seen (successfully passed) ---
    local seen_phases = Set{Int}() 

    for j in 1:(length(s_vals_desc) - 1)
        s_hi = s_vals_desc[j]; s_lo = s_vals_desc[j+1]
        s_interval = interval(s_lo, s_hi)

        tau_rect = find_new_tau_range(p_interval, s_interval, t_up_val, target_precision)
        if isempty_interval(tau_rect)
            println("!!! FAIL: Empty Tau for P=$p_interval, S=$s_interval")
            p_interval_failed = true; break
        end

        # --- ATTEMPT CURRENT PHASE ON FULL RECTANGLE ---
        passed, _ = calculate_and_format_check(current_phase, p_interval, s_interval, tau_rect, min_delta_endpoint, g_h_precision, false)

        if passed
            push!(seen_phases, current_phase) # Track success
            continue 
        end

        # --- FAILED CURRENT PHASE. TRY NEXT PHASE? ---
        can_advance = false
        if current_phase < 5
             passed_next, _ = calculate_and_format_check(current_phase + 1, p_interval, s_interval, tau_rect, min_delta_endpoint, g_h_precision, false)
             if passed_next
                 # Track next phase success
                 push!(seen_phases, current_phase + 1)
                 current_phase += 1
                 can_advance = true
             end
        end

        if can_advance
            continue
        end

        # --- BOTH FAILED. SUBDIVIDE ---
        should_subdivide = (inf(tau_rect) * 2 < sup(tau_rect)) || ((sup(tau_rect) - inf(tau_rect)) > tau_subdivision_threshold)
        
        if !should_subdivide
             p_interval_failed = true; break
        end

        phase_name_map = Dict(1=>"H", 2=>"G", 3=>"Delta", 4=>"G", 5=>"H")
        curr_name = phase_name_map[current_phase]
        println("--- $curr_name-Check failed on full tau interval. Subdividing for $curr_name... ---")
        
        sub_taus = subdivide_tau_interval(tau_rect, tau_subdivisions)
        all_subs_passed = true

        for (k, stau) in enumerate(sub_taus)
            # 1. Try current phase on sub-interval
            spass, _ = calculate_and_format_check(current_phase, p_interval, s_interval, stau, min_delta_endpoint, g_h_precision, true, k)
            
            if spass
                push!(seen_phases, current_phase) # Track success
                continue 
            end

            # 2. If that fails, try next phase (Transition boundary)
            if current_phase < 5
                spass_next, _ = calculate_and_format_check(current_phase + 1, p_interval, s_interval, stau, min_delta_endpoint, g_h_precision, true, k)
                if spass_next
                    push!(seen_phases, current_phase + 1) # Track success
                    current_phase += 1
                    continue
                end
            end

            all_subs_passed = false
            break
        end

        if !all_subs_passed
            p_interval_failed = true
            break
        end

    end # End Sigma Loop

    # --- CHECK PHASE COMPLETENESS ---
    # If the calculation itself didn't fail, check if we missed any phases (1 through 5)
    if !p_interval_failed
        required_phases = Set([1, 2, 3, 4, 5])
        missing_phases = setdiff(required_phases, seen_phases)
        
        if !isempty(missing_phases)
            println("!!! FAIL for Phase Completeness: P=$p_interval. Missing Phases: $missing_phases")
            p_interval_failed = true
        end
    end

    if p_interval_failed
        push!(failed_p_intervals, p_interval)
    end

end # End P Loop

println("="^80)
println("Final Report:")
if isempty(failed_p_intervals)
    println("All p-intervals processed successfully (Verification of Phases 1-5 complete).")
else
    merged = merge_intervals(failed_p_intervals)
    println("Inconclusive p-ranges:")
    for m in merged; println("  $m"); end
end
println("Total time: ", round(time() - start_time, digits=4), "s")