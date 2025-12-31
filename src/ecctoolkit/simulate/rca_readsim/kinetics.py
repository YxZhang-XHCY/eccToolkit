"""
Kinetic RCA repeat simulation (Gillespie/CTMC).

This model converts reaction time and state transitions into total
polymerized length, then maps to repeat count by length.
"""

import numpy as np

from .config import SimConfig


STATE_FREE = 0
STATE_BOUND = 1
STATE_ELONG = 2
STATE_PAUSE = 3
STATE_OFF = 4
STATE_INACT = 5


def _length_scale(length_ratio: float, exponent: float, params, inverse: bool = False) -> float:
    if exponent == 0:
        return 1.0

    scale = 1.0
    if length_ratio > 1.0:
        mode = params.len_scale_mode
        if mode == "exp":
            scale = float(np.exp(exponent * (length_ratio - 1.0)))
        elif mode == "saturating":
            tau = max(params.len_scale_tau, 1e-6)
            scale = 1.0 + exponent * (1.0 - float(np.exp(-(length_ratio - 1.0) / tau)))
        else:
            scale = length_ratio ** exponent

    if inverse:
        if scale <= 0:
            return 0.0
        return 1.0 / scale
    return scale


def _compute_v_eff(
    base_v: float,
    length_ratio: float,
    params,
    apply_length: bool = True,
) -> float:
    substrate_factor = 1.0
    if params.dntp is not None and params.k_dntp > 0:
        substrate_factor = params.dntp / (params.k_dntp + params.dntp)

    inhib_factor = 1.0
    if params.ppi is not None and params.k_ppi > 0:
        inhib_factor = 1.0 / (1.0 + params.ppi / params.k_ppi)

    length_factor = 1.0
    if apply_length:
        length_factor = _length_scale(length_ratio, params.len_exp_v, params, inverse=True)
    return max(0.0, base_v * substrate_factor * inhib_factor * length_factor)


def simulate_repeat_count(
    length: int,
    config: SimConfig,
    rng: np.random.Generator,
) -> int:
    params = config.kinetics
    if length <= 0:
        return config.R_min
    if params.rebinding_target != "same_molecule":
        raise NotImplementedError(
            f"rebinding_target={params.rebinding_target} is not supported yet"
        )

    length_ratio = max(length / config.L_ref, 1e-6)

    penalty_mode = params.length_penalty_mode
    apply_on = penalty_mode == "full"
    apply_pause = penalty_mode in ("full", "pause_off")
    apply_off = penalty_mode in ("full", "pause_off", "off_only")
    apply_v = penalty_mode in ("full", "speed_only")

    k_on = params.k_on
    if apply_on:
        k_on *= _length_scale(length_ratio, params.len_exp_on, params, inverse=True)

    k_pause = params.k_pause
    if apply_pause:
        k_pause *= _length_scale(length_ratio, params.len_exp_pause, params)

    k_off = params.k_off
    if apply_off:
        k_off *= _length_scale(length_ratio, params.len_exp_off, params)
    if params.k_off_pause is not None:
        k_off_pause = params.k_off_pause
    else:
        k_off_pause = k_off * params.off_pause_multiplier

    v_eff = _compute_v_eff(
        params.v_nt_per_min,
        length_ratio if apply_v else 1.0,
        params,
        apply_length=apply_v,
    )

    t = 0.0
    p_nt = 0.0
    state = STATE_FREE
    events = 0

    while t < params.reaction_time_min and state != STATE_INACT:
        if events >= params.max_events:
            break

        if state == STATE_FREE or state == STATE_OFF:
            total_rate = k_on + params.k_inact_free
            if total_rate <= 0:
                break
            dt = rng.exponential(1.0 / total_rate)
            dt_effective = min(dt, params.reaction_time_min - t)
            t += dt_effective
            if t >= params.reaction_time_min:
                break
            r = rng.random() * total_rate
            state = STATE_BOUND if r < k_on else STATE_INACT
            events += 1
            continue

        if state == STATE_BOUND:
            total_rate = params.k_start + params.k_inact_bound
            if total_rate <= 0:
                break
            dt = rng.exponential(1.0 / total_rate)
            dt_effective = min(dt, params.reaction_time_min - t)
            t += dt_effective
            if t >= params.reaction_time_min:
                break
            r = rng.random() * total_rate
            state = STATE_ELONG if r < params.k_start else STATE_INACT
            events += 1
            continue

        if state == STATE_ELONG:
            total_rate = k_pause + k_off + params.k_inact_bound
            if total_rate <= 0:
                break
            dt = rng.exponential(1.0 / total_rate)
            dt_effective = min(dt, params.reaction_time_min - t)
            if dt_effective > 0 and v_eff > 0:
                p_nt += v_eff * dt_effective
            t += dt_effective
            if t >= params.reaction_time_min:
                break
            r = rng.random() * total_rate
            if r < k_pause:
                state = STATE_PAUSE
            elif r < k_pause + k_off:
                state = STATE_OFF
            else:
                state = STATE_INACT
            events += 1
            continue

        if state == STATE_PAUSE:
            total_rate = params.k_resume + k_off_pause + params.k_inact_bound
            if total_rate <= 0:
                break
            dt = rng.exponential(1.0 / total_rate)
            dt_effective = min(dt, params.reaction_time_min - t)
            t += dt_effective
            if t >= params.reaction_time_min:
                break
            r = rng.random() * total_rate
            if r < params.k_resume:
                state = STATE_ELONG
            elif r < params.k_resume + k_off_pause:
                state = STATE_OFF
            else:
                state = STATE_INACT
            events += 1
            continue

        break

    repeat_raw = p_nt / length
    repeat_floor = int(np.floor(repeat_raw))
    remainder = repeat_raw - repeat_floor
    if remainder > 0 and rng.random() < remainder:
        repeat_floor += 1
    repeat_est = max(config.R_min, repeat_floor)
    return min(config.R_max, repeat_est)
