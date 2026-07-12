"""
SCIM Fault Transient Torque — calculation package.

Modules
-------
config  : Load and normalize motor parameters from a JSONC input file.
circuit : Equivalent-circuit steady-state analysis and slip derivation.
quick   : Closed-form analytical short-circuit torque estimate.
dq      : State-space dq-model time-domain simulation.
sweep   : Sweep fault inception angle to find worst-case torque.
compare : Compute comparison metrics between quick and dq results.

Usage
-----
    from scim_calc.config import load_jsonc, normalize_params
    from scim_calc.quick import run_quick_calc
    from scim_calc.dq import run_dq_simulation

    cfg = load_jsonc("input.jsonc")
    p = normalize_params(cfg)
    q = run_quick_calc(p)
    d = run_dq_simulation(p)
"""
