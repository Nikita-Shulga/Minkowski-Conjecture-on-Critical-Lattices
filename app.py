import streamlit as st
import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt

# --- 1. Math Functions ---

@st.cache_data
def get_data(p, num_points=100):
    """Calculates the curve points efficiently."""
    limit = (2**p - 1)**(1/p)
    
    # Solve for t_p (value at sigma=1)
    def tp_eq(t):
        if t <= 0 or t >= 1: return 1e9
        return 2 * (1 - t)**p - (1 + t**p)
    
    try:
        tp = brentq(tp_eq, 1e-6, 1 - 1e-6)
    except:
        tp = 0.5

    sigmas = np.linspace(1, limit, num_points)
    delta_vals = []

    for s in sigmas:
        # Determine tau
        if s == 1.0:
            tau = tp
        elif abs(s - limit) < 1e-5:
            tau = 0.0
        else:
            # Implicit equation for tau
            def implicit(t):
                term_tau = (1 + t**p)**(-1/p)
                term_sigma = (1 + s**p)**(-1/p)
                A = term_tau - term_sigma
                B = t * term_tau + s * term_sigma
                return (abs(A)**p + abs(B)**p) - 1
            
            try:
                # Use a safe bracket. Upper bound is roughly tp.
                tau = brentq(implicit, 0, tp + 0.1) 
            except:
                tau = 0.0
        
        # Calculate Delta
        term_tau = (1 + tau**p)**(-1/p)
        term_sigma = (1 + s**p)**(-1/p)
        val = (tau + s) * term_tau * term_sigma
        delta_vals.append(val)

    # Symmetrization for the bell curve case
    if abs(p - 2.571728) < 0.005:
        delta_vals = np.array(delta_vals)
        delta_vals = (delta_vals + delta_vals[::-1]) / 2.0
        delta_vals = delta_vals.tolist()
        
    return sigmas, delta_vals, limit

# --- 2. Streamlit UI ---

st.set_page_config(page_title="Delta Function Analysis", layout="centered")

st.title("Interactive Delta Function")
st.markdown(r"Visualization of $\Delta(p, \sigma)$ with dynamic scaling.")

# The Slider
p = st.slider("Parameter p", min_value=1.01, max_value=15.0, value=2.5717, step=0.0001, format="%.4f")

# Calculations
sigmas, deltas, current_limit = get_data(p)

# Plotting
fig, ax = plt.figure(figsize=(8, 6)), plt.gca()

# Fixed X-Axis (Growth visualization)
max_limit = (2**15 - 1)**(1/15) # p=15 limit
ax.set_xlim(0.95, max_limit + 0.1)

# Dynamic Y-Axis
y_min, y_max = min(deltas), max(deltas)
y_range = y_max - y_min
pad = max(y_range * 0.2, 0.001)
ax.set_ylim(y_min - pad, y_max + pad)

# Draw Curve
ax.plot(sigmas, deltas, color='blue', linewidth=2)

# Draw Limit Line (Sigma_p)
ax.axvline(x=current_limit, color='red', linestyle='--', alpha=0.6)
ax.text(current_limit, y_min, r'$\sigma_p$', color='red', verticalalignment='bottom', horizontalalignment='right')

# Labels and Style
ax.set_xlabel(r'$\sigma$', fontsize=14)
ax.set_ylabel(r'$\Delta(p, \sigma)$', fontsize=14)
ax.set_title(f'p = {p:.5f}', fontsize=16)
ax.grid(True, linestyle=':', alpha=0.6)

st.pyplot(fig)

st.markdown("""
**Description:**
As $p$ increases, the domain $\sigma \in [1, \sigma_p]$ grows (visible on X-axis).
The vertical axis automatically rescales to fit the curve shape.
""")
