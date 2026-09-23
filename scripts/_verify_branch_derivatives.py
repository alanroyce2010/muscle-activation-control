"""One-off sympy check of the by-hand branch-derivative formulas used in
Phase 2 (docs/theory/phase2_linearization.md). Not part of the package --
run once to verify, keep for anyone who wants to re-derive/audit later.
"""

import sympy as sp

# --- Activation ODE branches ---
a, u, tau_act, tau_deact, astar = sp.symbols(
    "a u tau_act tau_deact astar", positive=True
)

adot_activating = (u - a) / (tau_act * (sp.Rational(1, 2) + sp.Rational(3, 2) * a))
adot_deactivating = (u - a) * (sp.Rational(1, 2) + sp.Rational(3, 2) * a) / tau_deact

da_da_activating = sp.diff(adot_activating, a).subs({u: astar, a: astar}).simplify()
da_du_activating = sp.diff(adot_activating, u).subs({u: astar, a: astar}).simplify()
da_da_deactivating = sp.diff(adot_deactivating, a).subs({u: astar, a: astar}).simplify()
da_du_deactivating = sp.diff(adot_deactivating, u).subs({u: astar, a: astar}).simplify()

print("=== activation ODE partials at equilibrium (u=a=astar) ===")
print("activating branch:  d(adot)/da =", da_da_activating, " d(adot)/du =", da_du_activating)
print("deactivating branch: d(adot)/da =", da_da_deactivating, " d(adot)/du =", da_du_deactivating)

hand_da_da_activating = -1 / (tau_act * (sp.Rational(1, 2) + sp.Rational(3, 2) * astar))
hand_da_du_activating = 1 / (tau_act * (sp.Rational(1, 2) + sp.Rational(3, 2) * astar))
hand_da_da_deactivating = -(sp.Rational(1, 2) + sp.Rational(3, 2) * astar) / tau_deact
hand_da_du_deactivating = (sp.Rational(1, 2) + sp.Rational(3, 2) * astar) / tau_deact

assert sp.simplify(da_da_activating - hand_da_da_activating) == 0
assert sp.simplify(da_du_activating - hand_da_du_activating) == 0
assert sp.simplify(da_da_deactivating - hand_da_da_deactivating) == 0
assert sp.simplify(da_du_deactivating - hand_da_du_deactivating) == 0
print("hand-derived formulas MATCH sympy for all 4 activation partials.\n")

# --- Force-velocity branches (Katz 1939 asymmetry via Thelen 2003 Eq. 6-7) ---
v, af, flen, k, afl = sp.symbols("v Af Flen K afl", positive=True)  # afl = a*f_l

f_shortening = afl * (v + k) / (k - v / af)
c = v * (2 + 2 / af) / (flen - 1)
f_lengthening = afl * (c * flen + k) / (k + c)

df_dv_short_at0 = sp.diff(f_shortening, v).subs(v, 0).simplify()
df_dv_long_at0 = sp.diff(f_lengthening, v).subs(v, 0).simplify()

print("=== force-velocity slope at v_norm=0 ===")
print("shortening-side slope:", df_dv_short_at0)
print("lengthening-side slope:", df_dv_long_at0)
ratio = sp.simplify(df_dv_long_at0 / df_dv_short_at0)
print("ratio (lengthening/shortening):", ratio)
assert ratio == 2
print("confirmed: lengthening-side slope is exactly 2x shortening-side slope"
      " (matches Thelen 2003's stated Katz 1939 citation).")

hand_short = afl * (1 + 1 / af) / k
hand_long = afl * (2 + 2 / af) / k
assert sp.simplify(df_dv_short_at0 - hand_short) == 0
assert sp.simplify(df_dv_long_at0 - hand_long) == 0
print("hand-derived formulas MATCH sympy for both force-velocity slopes.")
