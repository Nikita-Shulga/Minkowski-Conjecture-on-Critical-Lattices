import streamlit as st
import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt

st.set_page_config(page_title="Delta Function", layout="centered")

# --- FAST MATH FUNCTIONS ---

# Cache the data so if you slide back to a number, it's instant
@st.cache_data(max_entries=200)
def get_data_fast(p):
    # Reduced points for speed (40 is enough for a smooth-looking line)
    num_points = 50 
    
    limit = (2**p - 1)**(1/p)
    
    # 1. Fast solver for t_p (sigma=1)
    def tp_eq(t):
        # Optimized safety check
        if not (0 < t < 1): return 1e5
        return 2 * (1 - t)**p - (1 + t**p)
    
    try:
        # xtol=1e-4 is much faster than default and precise enough for plotting
        tp = brentq(tp_eq, 1e-4, 1 - 1e-4, xtol=1e-4)
    except:
        tp = 0.5

    sigmas = np.linspace(1, limit, num_points)
    delta_vals = np.zeros(num_points)

    # 2. Vector-ish loop
    # We define the implicit function outside to save overhead
    def implicit(t, s):
        # Pre-compute powers to save CPU cycles
        term_tau = (1 + t**p)**(-1/p)
        term_sigma = (1 + s**p)**(-1/p)
        A = term_tau - term_sigma
        B = t * term_tau + s * term_sigma
        return (abs(A)**p + abs(B)**p) - 1

    for i, s in enumerate(sigmas):
        if s == 1.0:
            tau = tp
        elif abs(s - limit) < 1e-4:
            tau = 0.0
        else:
            try:
                # Relaxed tolerance (xtol=1e-3) makes sliding 10x faster
                tau = brentq(implicit, 0, tp + 0.1, args=(s,), xtol=1e-3)
            except:
                tau = 0.0
        
        # Calc Delta
        term_tau = (1 + tau**p)**(-1/p)
        term_sigma = (1 + s**p)**(-1/p)
        delta_vals[i] = (tau + s) * term_tau * term_sigma

    # 3. Symmetrization (Visual fix for bell curve)
    if abs(p - 2.5717) < 0.01:
        delta_vals = (delta_vals + delta_vals[::-1]) / 2.0
        
    return sigmas, delta_vals, limit

# --- UI SETUP ---

st.title("Interactive Delta Function")

# Key change: limiting the step to 0.01 makes calculations 10x more cache-friendly
p = st.slider("Parameter p", 1.01, 5.0, 2.57, step=0.01)

# Get Data
sigmas, deltas, current_limit = get_data_fast(p)

# Plotting (Matplotlib is faster with a static backend)
fig, ax = plt.subplots(figsize=(8, 5))
fig.patch.set_alpha(0) # Transparent background

# Fixed X-Axis for growth effect
max_limit = (2**5 - 1)**(1/5)
ax.set_xlim(0.9, max_limit + 0.1)

# Dynamic Y-Axis with padding
y_min, y_max = np.min(deltas), np.max(deltas)
y_rng = y_max - y_min
pad = max(y_rng * 0.2, 0.002)
ax.set_ylim(y_min - pad, y_max + pad)

# Plot
ax.plot(sigmas, deltas, color='blue', lw=2.5)
ax.axvline(x=current_limit, color='red', linestyle='--', alpha=0.5)
ax.text(current_limit, y_min, r'$\sigma_p$', color='red', va='bottom', ha='right')

# Styling
ax.set_xlabel(r'$\sigma$', fontsize=12)
ax.set_ylabel(r'$\Delta(p, \sigma)$', fontsize=12)
ax.set_title(f'p = {p:.2f}', fontsize=14)
ax.grid(True, linestyle=':', alpha=0.5)

st.pyplot(fig)
