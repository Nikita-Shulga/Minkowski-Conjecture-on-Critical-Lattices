using IntervalArithmetic

"""
    get_L(val_int, p_int)

Calculates the logarithmic derivative term L = d/dp [ -1/p * ln(1 + val^p) ].
Derived from C_p = C * L.
"""
function get_L(val_int::Interval{<:Real}, p_int::Interval{<:Real})
    # L = ln(1+x^p)/p^2 - (x^p * ln(x)) / (p * (1+x^p))
    
    val_pow_p = val_int^p_int
    ln_val = log(val_int)
    ln_one_plus_pow = log(interval(1) + val_pow_p)
    
    term1 = ln_one_plus_pow / (p_int^2)
    term2 = (val_pow_p * ln_val) / (p_int * (interval(1) + val_pow_p))
    
    return term1 - term2
end

"""
    calculate_k_interval(p_int, s_int, t_int)

Calculates the interval for k = d/dp g(p, τ).
Includes k1, k2, k3, k4 and the definition of N_p for sigma_p.
"""
function calculate_k_interval(
    p_int::Interval{<:Real},
    s_int::Interval{<:Real},
    t_int::Interval{<:Real}
)
    # Ensure strict positivity
    
    # --- 1. Basic Terms ---
    one_over_p = interval(1) / p_int
    p_minus_1 = p_int - interval(1)
    p_minus_2 = p_int - interval(2)
    
    val_1_sp = interval(1) + s_int^p_int
    val_1_tp = interval(1) + t_int^p_int
    
    # C1, C2
    C1 = val_1_sp^(-one_over_p)
    C2 = val_1_tp^(-one_over_p)
    
    # A, B
    A = C2 - C1
    
    B = t_int * C2 + s_int * C1
    
    # D1, D2
    t_pow_pm1 = t_int^p_minus_1
    s_pow_pm1 = s_int^p_minus_1
    
    D1 = interval(1) - s_int * t_pow_pm1
    D2 = interval(1) - t_int * s_pow_pm1
    
    # H1, H2
    A_pow_pm1 = A^p_minus_1
    B_pow_pm1 = B^p_minus_1
    
    term_s_A = s_pow_pm1 * A_pow_pm1
    term_t_A = t_pow_pm1 * A_pow_pm1
    
    H1 = B_pow_pm1 + term_s_A
    H2 = B_pow_pm1 - term_t_A

    # --- 2. New p-Derivative Terms ---
    
    # L_sigma, L_tau
    L_sigma = get_L(s_int, p_int)
    L_tau   = get_L(t_int, p_int)
    
    # C1_p, C2_p
    C1_p = C1 * L_sigma
    C2_p = C2 * L_tau
    
    # A_p, B_p
    A_p = C2 * L_tau - C1 * L_sigma
    B_p = t_int * C2 * L_tau + s_int * C1 * L_sigma
    
    # D1_p, D2_p
    ln_tau = log(t_int)
    ln_sigma = log(s_int)
    
    D1_p = -s_int * t_pow_pm1 * ln_tau
    D2_p = -t_int * s_pow_pm1 * ln_sigma
    
    # H1_p, H2_p preparation
    ln_B = log(B)
    ln_A = log(A)
    ln_s = log(s_int)
    ln_t = log(t_int)
    
    B_pow_pm2 = B^p_minus_2
    A_pow_pm2 = A^p_minus_2
    
    term_H_common_B = B_pow_pm1 * ln_B + p_minus_1 * B_pow_pm2 * B_p
    term_H_common_A_inner = A_pow_pm1 * ln_A + p_minus_1 * A_pow_pm2 * A_p
    
    # H1_p
    term_H1_s = s_pow_pm1 * A_pow_pm1 * ln_s + s_pow_pm1 * term_H_common_A_inner
    H1_p = term_H_common_B + term_H1_s
    
    # H2_p
    term_H2_t = t_pow_pm1 * A_pow_pm1 * ln_t + t_pow_pm1 * term_H_common_A_inner
    H2_p = term_H_common_B - term_H2_t

    # --- 3. Compute k1, k2, k3 ---
    
    # k1 = C1_p * D1 * H1 - C2_p * D2 * H2
    k1 = C1_p * D1 * H1 - C2_p * D2 * H2
    
    # k2 = C1 * D1_p * H1 - C2 * D2_p * H2
    k2 = C1 * D1_p * H1 - C2 * D2_p * H2
    
    # k3 = C1 * D1 * H1_p - C2 * D2 * H2_p
    k3 = C1 * D1 * H1_p - C2 * D2 * H2_p
    
    # --- 4. Compute k4 (Sigma Partial Terms) ---
    
    K_sigma = val_1_sp^(-interval(1) - one_over_p)
    K_tau   = val_1_tp^(-interval(1) - one_over_p)
    
    # Component Partial Derivatives w.r.t Sigma
    dC1_dsigma = -s_int^(p_int - interval(1)) * K_sigma
    dC2_dsigma = interval(0.0)
    
    dD1_dsigma = -t_pow_pm1
    dD2_dsigma = -(p_minus_1) * t_int * s_int^p_minus_2
    
    dA_dsigma = s_pow_pm1 * K_sigma
    dB_dsigma = K_sigma
    
    # dH1/dsigma
    term_dH1_1 = B_pow_pm2 * dB_dsigma
    term_dH1_2 = s_int^p_minus_2 * A_pow_pm1
    term_dH1_3 = s_pow_pm1 * A_pow_pm2 * dA_dsigma
    dH1_dsigma = p_minus_1 * (term_dH1_1 + term_dH1_2 + term_dH1_3)
    
    # dH2/dsigma
    term_dH2_1 = B_pow_pm2 * dB_dsigma
    term_dH2_2 = t_pow_pm1 * A_pow_pm2 * dA_dsigma
    dH2_dsigma = p_minus_1 * (term_dH2_1 - term_dH2_2)
    
    # Assemble partial g / partial sigma
    # dg/ds = (C1)_s D1 H1 + C1 (D1)_s H1 + C1 D1 (H1)_s - [ (C2)_s D2 H2 + C2 (D2)_s H2 + C2 D2 (H2)_s ]
    
    # Term 1 group (positive part of g)
    part_pos_1 = dC1_dsigma * D1 * H1
    part_pos_2 = C1 * dD1_dsigma * H1
    part_pos_3 = C1 * D1 * dH1_dsigma
    dg_pos = part_pos_1 + part_pos_2 + part_pos_3
    
    # Term 2 group (negative part of g)
    # Note: dC2_dsigma is 0, so first sub-term is 0
    part_neg_2 = C2 * dD2_dsigma * H2
    part_neg_3 = C2 * D2 * dH2_dsigma
    dg_neg = part_neg_2 + part_neg_3
    
    partial_g_sigma = dg_pos - dg_neg
    
    # --- Calculate Sigma_p ---
    
    # N_p = A^p ln A + B^p ln B + p A^{p-1} A_p + p B^{p-1} B_p
    
    # Reconstruct A^p and B^p from existing powers
    A_pow_p = A * A_pow_pm1 
    B_pow_p = B * B_pow_pm1
    
    # N_p Terms
    term_Np_1 = A_pow_p * ln_A
    term_Np_2 = B_pow_p * ln_B
    term_Np_3 = p_int * A_pow_pm1 * A_p
    term_Np_4 = p_int * B_pow_pm1 * B_p
    
    N_p = term_Np_1 + term_Np_2 + term_Np_3 + term_Np_4
    
    # sigma_p = - N_p / (p * K_sigma * H1)
    denominator_sigma_p = p_int * K_sigma * H1
    
    # RIGOR FIX: If denominator contains 0, we cannot bound sigma_p.
    if 0.0 in denominator_sigma_p
        return emptyinterval(Float64)
    end
    
    sigma_p = -N_p / denominator_sigma_p
    
    k4 = partial_g_sigma * sigma_p
    
    # --- 5. Total Sum ---
    k_total = k1 + k2 + k3 + k4
    
    return k_total
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




# --- Main Execution ---

# 1. Update Constants for New Region
const p_min = 1 // 1
const p_max = 101 // 100
const p_interval_fixed = interval(p_min, p_max)

const sigma_min = 1 // 1
const sigma_max = 101377243 // 100000000
const sigma_interval_fixed = interval(sigma_min, sigma_max)

const tau_start = 2 // 1000
const tau_end = 2 // 10
const tau_step = 1 // 100000
const k_threshold = - 1 // 10^9 # Inequality: k < -1e-9

println("Calculating bounds for k = d/dp g(p,τ).")
println("  p fixed range: $p_interval_fixed")
println("  σ fixed range: $sigma_interval_fixed")
println("  τ range: [$tau_start, $tau_end] with step $tau_step")
println("  Inequality to check: k < $k_threshold")
println("="^40)
println("-"^100)
println("  P-Rect             S-Rect             τ-sub-interval     k Bound")
println("-"^100)

start_time = time()
failed_tau_intervals = Interval{Float64}[]

# 2. Iterate over Tau
# Calculate number of steps
num_steps = Int(round((tau_end - tau_start) / tau_step))

for i in 0:(num_steps - 1)
    t_lo = tau_start + i * tau_step
    # Explicitly snap the last upper bound to the exact end constant
    if i == num_steps - 1
   	 t_hi = tau_end
    else
    	t_hi = tau_start + (i + 1) * tau_step
    end
    
    # Create the interval for this step
    tau_interval = interval(t_lo, t_hi)
    
    # Calculate k
    k_res = calculate_k_interval(p_interval_fixed, sigma_interval_fixed, tau_interval)
    
    # --- RIGOR FIX: Check for Singularity Failure ---
    if isempty_interval(k_res)
        println("      !!! FAIL: Singularity detected (sigma_p unbounded) at τ=$tau_interval")
        push!(failed_tau_intervals, tau_interval)
        continue
    end
    # ------------------------------------------------

    # PRINT RESULTS FOR EVERY STEP
    println("  P=$(p_interval_fixed) S=$(sigma_interval_fixed) τ=$(tau_interval) => k=$(k_res)")


    # Check Inequality k < threshold (strictly negative)
    if sup(k_res) < k_threshold
        # Passed
        continue
    else
        # Failed
        println("      !!! FAIL: τ=$tau_interval")
        println("          k = $k_res (Wanted hi < $k_threshold)")
        push!(failed_tau_intervals, tau_interval)
    end
end

end_time = time()
println("="^40)
println("Final Report on Inconclusive Tau-Intervals:")
if isempty(failed_tau_intervals)
    println("All tau-intervals verified k < -1e-9 successfully.")
else
    merged_failed = merge_intervals(failed_tau_intervals)
    println("The following tau-ranges failed the check:")
    for t_int in merged_failed
        println("  ", t_int)
    end
end
println("Total time: ", round(end_time - start_time, digits=4), " seconds.")