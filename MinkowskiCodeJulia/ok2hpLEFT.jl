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
    t_pow_a_log_t_safe(t_int, a_int)

Calculates t^a * ln(t) rigorously handling the singularity at t=0.
- If a > 0: Limits to 0 at t=0.
- If a <= 0: Behaves according to standard interval arithmetic (unbounded at 0).
"""
function t_pow_a_log_t_safe(t_int::Interval{<:Real}, a_int::Interval{<:Real})
    # 1. Purely in the safe region (away from 0)
    if inf(t_int) > 0
        return (t_int^a_int) * log(t_int)
    end

    # 2. Check Exponent for Limit existence
    # The limit of t^a * ln(t) as t->0 is 0 ONLY if a > 0.
    if inf(a_int) <= 0
        # If exponent can be 0 or negative, term is unbounded. 
        # Return standard evaluation to capture the Infinity.
        return (t_int^a_int) * log(t_int)
    end

    # 3. Handle Neighborhood of 0 (Monotonic Decreasing Region)
    # Turning point is exp(-1/a). Conservatively check if t is small.
    safe_monotonic_bound = 1 // 10
    
    if sup(t_int) < safe_monotonic_bound
        # RIGOROUS ANALYTIC EXTENSION:
        # Domain: [0, t_hi] -> Range: [t_hi^a * ln(t_hi), 0]
        
        t_hi_iv = interval(sup(t_int))
        val_at_end = (t_hi_iv^a_int) * log(t_hi_iv)
        
        # Return the interval [min_val, 0.0]
        return interval(inf(val_at_end), 0)
    else
        return (t_int^a_int) * log(t_int)
    end
end


"""
    calculate_hp_derivative_interval(p_int, s_int, t_int)

Calculates h'_p = d(sum h_i)/dp using total derivatives w.r.t p.
"""
function calculate_hp_derivative_interval(
    p_int::Interval{<:Real},
    s_int::Interval{<:Real},
    t_int::Interval{<:Real}
)
    # We use t_int directly.
    
    # --- 0. Basic Constants & Terms ---
    one_over_p = interval(1) / p_int
    p_minus_1 = p_int - interval(1)
    p_minus_2 = p_int - interval(2)
    
    ln_s = log(s_int)
    # ln_t = log(t_int)  <-- Implicitly handled in specific terms via function call
    
    val_1_sp = interval(1) + s_int^p_int
    val_1_tp = interval(1) + t_int^p_int
    
    C1 = val_1_sp^(-one_over_p)
    C2 = val_1_tp^(-one_over_p)
    
    K_sigma = val_1_sp^(-interval(1) - one_over_p)
    K_tau = val_1_tp^(-interval(1) - one_over_p)
    
    # A, B
    A = C2 - C1
    B = t_int * C2 + s_int * C1 # Used t_int
    
    # Intervals A and B are used directly.
    
    t_pow_pm1 = t_int^p_minus_1
    t_pow_pm2 = t_int^p_minus_2
    s_pow_pm1 = s_int^p_minus_1
    s_pow_pm2 = s_int^p_minus_2
    
    A_pow_pm1 = A^p_minus_1
    A_pow_pm2 = A^p_minus_2
    B_pow_pm1 = B^p_minus_1
    B_pow_pm2 = B^p_minus_2
    
    D1 = interval(1) - s_int * t_pow_pm1
    D2 = interval(1) - t_int * s_pow_pm1
    
    H1 = B_pow_pm1 + s_pow_pm1 * A_pow_pm1
    H2 = B_pow_pm1 - t_pow_pm1 * A_pow_pm1
    
    # Check H2 safety for division
    if in_interval(0,H2)== true 
        println("  Warning: H_2 interval contains zero.")
        return emptyinterval(Float64)
    end
    
    # --- 1. Basic p-Derivatives (Partial) ---
    
    # (C1)_p
    term_C1_bracket = (log(val_1_sp) / (p_int^interval(2))) - ( (s_int^p_int * ln_s) / (p_int * val_1_sp) )
    C1_p = C1 * term_C1_bracket
    
    # (C2)_p_part
    # Safe handling of t^p * ln(t)
    term_tp_ln_t = t_pow_a_log_t_safe(t_int, p_int)
    term_C2_bracket = (log(val_1_tp) / (p_int^interval(2))) - ( term_tp_ln_t / (p_int * val_1_tp) )
    C2_p_part = C2 * term_C2_bracket
    
    # (K_sigma)_p
    term_Ks_bracket = (log(val_1_sp) / (p_int^interval(2))) - ( (interval(1) + one_over_p) * (s_int^p_int * ln_s) / val_1_sp )
    Ks_p = K_sigma * term_Ks_bracket
    
    # (K_tau)_p_part
    # Safe handling of t^p * ln(t)
    term_Kt_bracket = (log(val_1_tp) / (p_int^interval(2))) - ( (interval(1) + one_over_p) * term_tp_ln_t / val_1_tp )
    Kt_p_part = K_tau * term_Kt_bracket
    
    # --- 2. Calculate Tau_p ---
    
    # A_p_part, B_p_part
    Ap_part = C2_p_part - C1_p
    Bp_part = t_int * C2_p_part + s_int * C1_p # Used t_int
    
    # Phi_p
    # Using A directly (assumes A > 0, which is true for disjoint sigma/tau)
    term_Phi_p_A = (A^p_int) * ( log(A) + p_int * (Ap_part / A) )
    term_Phi_p_B = (B^p_int) * ( log(B) + p_int * (Bp_part / B) )
    Phi_p = term_Phi_p_A + term_Phi_p_B
    
    # Phi_tau
    Phi_tau = p_int * K_tau * ( B_pow_pm1 - A_pow_pm1 * t_pow_pm1 )
    
    # tau_p
    if 0.0 in Phi_tau
        return emptyinterval(Float64)
    end
    tau_p = - Phi_p / Phi_tau
    
    # --- 3. Total p-Derivatives ---
    
    # (C2)_p = (C2)_p_part - t^(p-1) K_tau tau_p
    C2_p = C2_p_part - t_pow_pm1 * K_tau * tau_p
    
    # (K_tau)_p = (K_tau)_p_part - (p+1) * (t^(p-1)/(1+t^p)) * K_tau * tau_p
    Kt_p = Kt_p_part - (p_int + interval(1)) * ( (t_pow_pm1 / val_1_tp) * K_tau * tau_p )
    
    # A_p = Ap_part - t^(p-1) K_tau tau_p
    A_p = Ap_part - t_pow_pm1 * K_tau * tau_p
    
    # B_p = Bp_part + K_tau tau_p
    B_p = Bp_part + K_tau * tau_p
    
    # --- 4. Power Derivatives ---
    
    # (sigma^(p-1))_p
    s_pow_pm1_p = s_pow_pm1 * ln_s
    # (sigma^(p-2))_p
    s_pow_pm2_p = s_pow_pm2 * ln_s
    
    # (tau^(p-1))_p = t^(p-1) (ln t + (p-1) tau_p / t)
    # Decomposed to: t^(p-1) ln t + (p-1) * tau_p * t^(p-2)
    # Safe handling of t^(p-1) * ln(t)
    t_pow_pm1_log_t = t_pow_a_log_t_safe(t_int, p_minus_1)
    t_pow_pm1_p = t_pow_pm1_log_t + p_minus_1 * tau_p * t_pow_pm2
    
    # (tau^(p-2))_p = t^(p-2) (ln t + (p-2) tau_p / t)
    # Decomposed to: t^(p-2) ln t + (p-2) * tau_p * t^(p-3)
    # Safe handling of t^(p-2) * ln(t) (Handles p-2 < 0 correctly)
    t_pow_pm2_log_t = t_pow_a_log_t_safe(t_int, p_minus_2)
    t_pow_pm3 = t_int^(p_int - interval(3))
    t_pow_pm2_p = t_pow_pm2_log_t + p_minus_2 * tau_p * t_pow_pm3
    
    # (B^(p-1))_p = B^(p-1) ln B + (p-1) B^(p-2) B_p
    B_pow_pm1_p = B_pow_pm1 * log(B) + p_minus_1 * B_pow_pm2 * B_p
    
    # (A^(p-1))_p = A^(p-1) ln A + (p-1) A^(p-2) A_p
    A_pow_pm1_p = A_pow_pm1 * log(A) + p_minus_1 * A_pow_pm2 * A_p
    
    # (B^(p-2))_p = B^(p-2) ln B + (p-2) B^(p-3) B_p
    B_pow_pm3 = B^(p_int - interval(3))
    B_pow_pm2_p = B_pow_pm2 * log(B) + p_minus_2 * B_pow_pm3 * B_p
    
    # (A^(p-2))_p = A^(p-2) ln A + (p-2) A^(p-3) A_p
    A_pow_pm3 = A^(p_int - interval(3))
    A_pow_pm2_p = A_pow_pm2 * log(A) + p_minus_2 * A_pow_pm3 * A_p
    
    # --- 5. H Derivatives ---
    
    # (H1)_p
    H1_p = B_pow_pm1_p + s_pow_pm1_p * A_pow_pm1 + s_pow_pm1 * A_pow_pm1_p
    
    # (H2)_p
    H2_p = B_pow_pm1_p - t_pow_pm1_p * A_pow_pm1 - t_pow_pm1 * A_pow_pm1_p
    
    # --- 6. T_p ---
    # S = K_sigma / K_tau
    S = K_sigma / K_tau
    S_p = (Ks_p * K_tau - K_sigma * Kt_p) / (K_tau^interval(2))
    
    # R = H1 / H2
    R = H1 / H2
    R_p = (H1_p * H2 - H1 * H2_p) / (H2^interval(2))
    
    # T = -S * R => T_p = - (S_p * R + S * R_p)
    T = -S * R
    T_p = - (S_p * R + S * R_p)
    
    # --- 7. Partial Derivatives Helper (S1, S2, S3, S4) and p-derivatives ---
    
    dB_dsigma = K_sigma
    dB_dsigma_p = Ks_p  # (K_sigma)_p
    
    dA_dsigma = s_pow_pm1 * K_sigma
    dA_dsigma_p = s_pow_pm1_p * K_sigma + s_pow_pm1 * Ks_p
    
    dB_dtau = K_tau
    dB_dtau_p = Kt_p    # (K_tau)_p
    
    dA_dtau = -t_pow_pm1 * K_tau
    dA_dtau_p = - (t_pow_pm1_p * K_tau + t_pow_pm1 * Kt_p)
    
    # S1
    S1 = B_pow_pm2 * dB_dsigma + s_pow_pm2 * A_pow_pm1 + s_pow_pm1 * A_pow_pm2 * dA_dsigma
    # (S1)_p
    term_S1_p_1 = B_pow_pm2_p * dB_dsigma + B_pow_pm2 * dB_dsigma_p
    term_S1_p_2 = s_pow_pm2_p * A_pow_pm1 + s_pow_pm2 * A_pow_pm1_p
    term_S1_p_3 = s_pow_pm1_p * A_pow_pm2 * dA_dsigma + s_pow_pm1 * A_pow_pm2_p * dA_dsigma + s_pow_pm1 * A_pow_pm2 * dA_dsigma_p
    S1_p = term_S1_p_1 + term_S1_p_2 + term_S1_p_3
    
    # S2
    S2 = B_pow_pm2 * dB_dsigma - t_pow_pm1 * A_pow_pm2 * dA_dsigma
    # (S2)_p
    term_S2_p_1 = B_pow_pm2_p * dB_dsigma + B_pow_pm2 * dB_dsigma_p
    term_S2_p_2 = t_pow_pm1_p * A_pow_pm2 * dA_dsigma + t_pow_pm1 * A_pow_pm2_p * dA_dsigma + t_pow_pm1 * A_pow_pm2 * dA_dsigma_p
    S2_p = term_S2_p_1 - term_S2_p_2
    
    # S3
    S3 = B_pow_pm2 * dB_dtau + s_pow_pm1 * A_pow_pm2 * dA_dtau
    # (S3)_p
    term_S3_p_1 = B_pow_pm2_p * dB_dtau + B_pow_pm2 * dB_dtau_p
    term_S3_p_2 = s_pow_pm1_p * A_pow_pm2 * dA_dtau + s_pow_pm1 * A_pow_pm2_p * dA_dtau + s_pow_pm1 * A_pow_pm2 * dA_dtau_p
    S3_p = term_S3_p_1 + term_S3_p_2
    
    # S4
    S4 = B_pow_pm2 * dB_dtau - t_pow_pm2 * A_pow_pm1 - t_pow_pm1 * A_pow_pm2 * dA_dtau
    # (S4)_p
    term_S4_p_1 = B_pow_pm2_p * dB_dtau + B_pow_pm2 * dB_dtau_p
    term_S4_p_2 = t_pow_pm2_p * A_pow_pm1 + t_pow_pm2 * A_pow_pm1_p
    term_S4_p_3 = t_pow_pm1_p * A_pow_pm2 * dA_dtau + t_pow_pm1 * A_pow_pm2_p * dA_dtau + t_pow_pm1 * A_pow_pm2 * dA_dtau_p
    S4_p = term_S4_p_1 - term_S4_p_2 - term_S4_p_3
    
    # H partials p-derivatives
    H1_s = p_minus_1 * S1
    H1_s_p = S1 + p_minus_1 * S1_p
    
    H2_s = p_minus_1 * S2
    H2_s_p = S2 + p_minus_1 * S2_p
    
    H1_t = p_minus_1 * S3
    H1_t_p = S3 + p_minus_1 * S3_p
    
    H2_t = p_minus_1 * S4
    H2_t_p = S4 + p_minus_1 * S4_p
    
    # --- 8. Calculate (h_i)_p ---
    
    # h1 = -sigma^(p-1) K_sigma D1 H1 - tau^(p-1) K_sigma D2 H1
    # D derivatives:
    # (D1)_p = -sigma t^(p-1)_p
    D1_p = -s_int * t_pow_pm1_p
    
    # (D2)_p = -(tau)_p sigma^(p-1) - tau (sigma^(p-1))_p
    D2_p = -( tau_p * s_pow_pm1 + t_int * s_pow_pm1_p ) # Used t_int
    
    # U1 = s^(p-1) K_sigma D1 H1
    U1 = s_pow_pm1 * K_sigma * D1 * H1
    U1_p = s_pow_pm1_p * K_sigma * D1 * H1 +
           s_pow_pm1 * Ks_p * D1 * H1 +
           s_pow_pm1 * K_sigma * D1_p * H1 +
           s_pow_pm1 * K_sigma * D1 * H1_p
           
    # U2 = t^(p-1) K_sigma D2 H1
    U2 = t_pow_pm1 * K_sigma * D2 * H1
    U2_p = t_pow_pm1_p * K_sigma * D2 * H1 +
           t_pow_pm1 * Ks_p * D2 * H1 +
           t_pow_pm1 * K_sigma * D2_p * H1 +
           t_pow_pm1 * K_sigma * D2 * H1_p
           
    h1_p = - (U1_p + U2_p)
    
    # h2 = - t^(p-1) C1 H1
    h2_p = - ( t_pow_pm1_p * C1 * H1 +
               t_pow_pm1 * C1_p * H1 +
               t_pow_pm1 * C1 * H1_p )
               
    # h3 = (p-1) tau sigma^(p-2) C2 H2
    F = p_minus_1; G = t_int; H = s_pow_pm2; I = C2; J = H2
    G_p = tau_p
    
    term_h3_1 = interval(1) * G * H * I * J
    term_h3_2 = F * G_p * H * I * J
    term_h3_3 = F * G * s_pow_pm2_p * I * J
    term_h3_4 = F * G * H * C2_p * J
    term_h3_5 = F * G * H * I * H2_p
    
    h3_p = term_h3_1 + term_h3_2 + term_h3_3 + term_h3_4 + term_h3_5
    
    # h4 = -(p-1) sigma tau^(p-2) C1 H1 T
    F = p_minus_1; G = s_int; H = t_pow_pm2; I = C1; J = H1; Q = T
    
    term_h4_1 = interval(1) * G * H * I * J * Q
    # G_p is 0
    term_h4_3 = F * G * t_pow_pm2_p * I * J * Q
    term_h4_4 = F * G * H * C1_p * J * Q
    term_h4_5 = F * G * H * I * H1_p * Q
    term_h4_6 = F * G * H * I * J * T_p
    
    h4_p = - (term_h4_1 + term_h4_3 + term_h4_4 + term_h4_5 + term_h4_6)
    
    # h5 = sigma^(p-1) C2 H2 T
    R_val = s_pow_pm1; S_val = C2; U_val = H2; V_val = T
    
    term_h5_1 = s_pow_pm1_p * S_val * U_val * V_val
    term_h5_2 = R_val * C2_p * U_val * V_val
    term_h5_3 = R_val * S_val * H2_p * V_val
    term_h5_4 = R_val * S_val * U_val * T_p
    
    h5_p = term_h5_1 + term_h5_2 + term_h5_3 + term_h5_4
    
    # h6 = C1 D1 H1_sigma - C2 D2 H2_sigma
    E1_p = C1_p * D1 * H1_s +
           C1 * D1_p * H1_s +
           C1 * D1 * H1_s_p
           
    E2_p = C2_p * D2 * H2_s +
           C2 * D2_p * H2_s +
           C2 * D2 * H2_s_p
           
    h6_p = E1_p - E2_p
    
    # h7 = (C1 D1 H1_tau - C2 D2 H2_tau) * T
    t1_p = C1_p * D1 * H1_t + C1 * D1_p * H1_t + C1 * D1 * H1_t_p
    t2_p = C2_p * D2 * H2_t + C2 * D2_p * H2_t + C2 * D2 * H2_t_p
    
    E_val = C1 * D1 * H1_t - C2 * D2 * H2_t
    E_p = t1_p - t2_p
    
    h7_p = E_p * T + E_val * T_p
    
    # --- 9. Final Sum ---
    hp_total = h1_p + h2_p + h3_p + h4_p + h5_p + h6_p + h7_p
    
    return hp_total
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
            println("!!! FAIL: Empty Tau for P=$p_interval, S=$s_interval")
            return false
    end
    
    # 2. Standard Check
    deriv_res = calculate_hp_derivative_interval(p_int, s_int, tau_interval_rect)
    # Catch the singularity from step 1
    if isempty_interval(deriv_res)
            println("!!! FAIL: Singularity detected in H2 for P=$p_interval, S=$s_interval")
            return false
    end
        
    # PRINT RESULT BEFORE CHECKING
    println("$(desc_prefix)  P=$(p_int) S=$(s_int) τ=$(tau_interval_rect) => h'_p=$(deriv_res)")



    # CHECK: lower bound > threshold (strictly positive)
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
                interval(0, split_point),
                interval(split_point, sup(tau_interval_rect))
            ]
        else
            println("$(desc_prefix)--- Check failed on full tau interval. Subdividing... ---")
            tau_subs = subdivide_tau_interval(tau_interval_rect, tau_subdivisions)
        end
        
        for (k, sub_tau) in enumerate(tau_subs)
            deriv_sub = calculate_hp_derivative_interval(p_int, s_int, sub_tau)
            
            sub_label = is_zero_bound ? "(Sub-Zero)" : "(Sub)"
            println("$(desc_prefix)      $sub_label P=$p_int S=$s_int τ[$k]=$sub_tau => h'_p=$deriv_sub")
            
            if !(inf(deriv_sub) > threshold)
                println("$(desc_prefix)      !!! SUB-FAIL: P=$p_int, S=$s_int, τ_sub=$sub_tau")
                println("$(desc_prefix)          h'_p = $deriv_sub (Wanted lo > $threshold)")
                return false
            end
        end
        return true
    else
        println("$(desc_prefix)!!! FAIL: P=$p_int, S=$s_int, τ=$tau_interval_rect")
        println("$(desc_prefix)    h'_p = $deriv_res (Wanted lo > $threshold)")
        return false
    end
end


# --- Main Execution ---

const p_start = 199954 // 100000
const p_end = 2 // 1
const p_step = 1 // 100000

const sigma_start_fixed = 1 // 1
const sigma_end_fixed = 10016 // 10000
const sigma_step = 1 // 100000


const tp_guess = interval(0, 36//100) # Initial guess for t_p
const target_precision = 1 // 10^9 # Iteration precision
const deriv_threshold  =  1 // 10^9 # Inequality: h > 1e-9
const tau_subdivisions = 40 #  Number of subdivisions for tau
const tau_subdivision_threshold = 1 // 10^4 # Absolute width threshold


println("Calculating derivative bounds h'_p over new grid.")
println("  p range: [$p_start, $p_end] (descending)")
println("  σ range: [$sigma_start_fixed, $sigma_end_fixed] (descending)")
println("  Inequality to check: h'_p > $deriv_threshold")
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
    
    # 2. Determine Sigma Range for this P (Not used dynamically here, fixed range used)
    
    println("\nP-Interval: ", p_interval)
    println("-"^100)
    println("  P-Rect             S-Rect              τ-sub-interval     h'_p Bound")
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
    println("All p-intervals verified h'_p > 1e-9 successfully.")
else
    merged_failed = merge_intervals(failed_p_intervals)
    println("The following p-ranges failed the check:")
    for p_int in merged_failed
        println("  ", p_int)
    end
end
println("Total time: ", round(end_time - start_time, digits=4), " seconds.")
