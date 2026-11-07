# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama
"""Write the golden vectors for dynamics-core from an independent reference.

The oracle shares no floating-point code with the Rust crate. It starts from the
exact values of the double-precision inputs and works in Python's decimal module
at 400 significant digits, then rounds each result once to the nearest double:

- peaks come from iterating the recurrence A(n+1) = A(n)·r + h, not from the
  closed form the Rust code uses (Lean proves the two equal);
- pulses_to_fraction is ceil(-ln(1 - p)·τ/T), then checked against its
  definition, 1 - r^n >= p at n but not at n - 1;
- everything else is the defining formula, with decimal exp and ln.

Run:

    python vectors/generate.py
"""

import json
import math
from decimal import ROUND_CEILING, Decimal, getcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
getcontext().prec = 400
TINY = Decimal("1e-100")  # below this, exp(x) - 1 and ln(1 + x) use their series
EXACT = {"rel": 0, "abs": 0}
AGREE = "dynamics:IMPLEMENTATIONS-AGREE"


def write(
    name: str,
    description: str,
    cases: list[dict],
    rel: float = 1e-13,
    abs_: float = 0.0,
) -> None:
    doc = {
        "schema": "golden-vectors@1",
        "description": description,
        "tolerance": {"rel": rel, "abs": abs_},
        "cases": cases,
    }
    (HERE / f"{name}.json").write_text(
        json.dumps(doc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


# The reference implementation (exact inputs, 400 digits) -------------------------


def series(x: Decimal, term) -> Decimal:
    """Sum of term(x, k) for k = 1, 2, ... until the terms stop mattering."""
    total, k = Decimal(0), 1
    while True:
        t = term(x, k)
        if t == 0 or abs(t) < abs(total) * Decimal("1e-420"):
            return total
        total, k = total + t, k + 1


def expm1(x: Decimal) -> Decimal:
    if abs(x) < TINY:
        return series(x, lambda x, k: x**k / math.factorial(k))
    return x.exp() - 1


def ln1p(x: Decimal) -> Decimal:
    if abs(x) < TINY:
        return series(x, lambda x, k: (-1) ** (k + 1) * x**k / k)
    return (1 + x).ln()


def time_constants(period: float, tau: float) -> Decimal:
    return Decimal(period) / Decimal(tau)


def decay_factor(dt: float, tau: float) -> float:
    return float((-time_constants(dt, tau)).exp())


def relax(a: float, dt: float, tau: float) -> float:
    return float(Decimal(a) * (-time_constants(dt, tau)).exp())


def half_life(tau: float) -> float:
    return float(Decimal(tau) * Decimal(2).ln())


def peak_after_pulses(h: float, period: float, tau: float, n: int) -> float:
    r, a = (-time_constants(period, tau)).exp(), Decimal(0)
    for _ in range(n):  # the recurrence itself: relax for one period, then add h
        a = a * r + Decimal(h)
    return float(a)


def ratio(x: Decimal) -> Decimal:
    """1 / (1 - e^(-x))."""
    return 1 / -expm1(-x)


def summation_ratio(period: float, tau: float) -> float:
    return float(ratio(time_constants(period, tau)))


def summation_ratio_at_hz(hz: float, tau: float) -> float:
    return float(ratio(1 / (Decimal(hz) * Decimal(tau))))  # the exact period 1/hz


def steady_state_peak(h: float, period: float, tau: float) -> float:
    return float(Decimal(h) * ratio(time_constants(period, tau)))


def steady_state_trough(h: float, period: float, tau: float) -> float:
    x = time_constants(period, tau)
    return float(Decimal(h) * (-x).exp() * ratio(x))


def reaches(p: float, x: Decimal, n: int) -> bool:
    """The definition: the n-th peak is at least p of the steady peak, 1 - r^n >= p."""
    return -expm1(-n * x) >= Decimal(p)


def pulses_to_fraction(p: float, period: float, tau: float) -> int:
    x = time_constants(period, tau)
    n = max(1, int((-ln1p(-Decimal(p)) / x).to_integral_value(rounding=ROUND_CEILING)))
    assert reaches(p, x, n) and (n == 1 or not reaches(p, x, n - 1)), (
        p,
        period,
        tau,
        n,
    )
    return n


def spike_shape(s: float, rise: float, decay: float) -> float:
    s_ = Decimal(s)
    return float(-expm1(-s_ / Decimal(rise)) * (-s_ / Decimal(decay)).exp())


def spike_peak_time(rise: float, decay: float) -> float:
    return float(Decimal(rise) * ln1p(Decimal(decay) / Decimal(rise)))


ORACLE = {
    f.__name__: f
    for f in (
        decay_factor,
        relax,
        half_life,
        peak_after_pulses,
        steady_state_peak,
        steady_state_trough,
        summation_ratio,
        summation_ratio_at_hz,
        pulses_to_fraction,
        spike_shape,
        spike_peak_time,
    )
}


def naive_pulses(p: float, period: float, tau: float) -> int:
    """The obvious floating-point formula, ceil(-ln(1 - p)·τ/T), which is off by one on boundaries."""
    return max(1, math.ceil(-math.log1p(-p) / (period / tau)))


def case(
    claim: str, fn: str, *args, note: str | None = None, tolerance: dict | None = None
) -> dict:
    """A case whose expected value comes from the oracle."""
    c = {"claim": claim, "fn": fn, "args": list(args), "expected": ORACLE[fn](*args)}
    if tolerance is not None:
        c["tolerance"] = tolerance
    if note:
        c["note"] = note
    return c


def error(fn: str, *args, note: str) -> dict:
    return {
        "claim": "dynamics:INPUT-VALIDATION",
        "fn": fn,
        "args": list(args),
        "error": True,
        "note": note,
    }


# The cases --------------------------------------------------------------------

TAU = 75.0  # s, the CaMKII-like time constant of the project's egg model (dynamics:TAU-75S)
TEN_MIN, THIRTY_FIVE_MIN = (
    600.0,
    2100.0,
)  # s, the human spike intervals (Taylor et al. 1993)
HUMAN = "dynamics:HUMAN-SPIKES-BARELY-ADD-UP"


def relaxation() -> list[dict]:
    compose = "dynamics:RELAX-COMPOSE"
    return [
        case(compose, "decay_factor", 0.0, TAU, note="no time, no decay"),
        case(compose, "decay_factor", TAU, TAU, note="one time constant: e^-1"),
        case(
            HUMAN,
            "decay_factor",
            TEN_MIN,
            TAU,
            note="r = e^-8 for spikes 10 minutes apart",
        ),
        case(
            HUMAN,
            "decay_factor",
            THIRTY_FIVE_MIN,
            TAU,
            note="r = e^-28 for spikes 35 minutes apart",
        ),
        case(compose, "decay_factor", 1.0, 1e6),
        case(compose, "decay_factor", 1e-12, 1.0),
        case(
            compose, "decay_factor", 800.0, 1.0, note="e^-800 is below every double: 0"
        ),
        case(compose, "relax", 2.5, 30.0, TAU, note="30 s ..."),
        case(compose, "relax", 2.5, 45.0, TAU, note="... then 45 s more ..."),
        case(compose, "relax", 2.5, 75.0, TAU, note="... is the same as 75 s at once"),
        case(compose, "relax", 0.0, 5.0, 1.0, note="an empty integrator stays empty"),
        case(compose, "relax", 1e300, 1.0, 1.0),
        case(
            compose,
            "relax",
            1e300,
            750.0,
            1.0,
            note="e^-750 alone is 0 in f64, but the state left, 1.9e-26, is not",
        ),
        case(
            compose,
            "relax",
            1e10,
            720.0,
            1.0,
            note="e^-720 alone is subnormal, but the state left is a normal double",
        ),
        case(compose, "relax", 1e300, 1600.0, 1.0, note="below every double: 0"),
        case("dynamics:HALF-LIFE", "half_life", TAU, note="75·ln 2 = 51.99 s"),
        case("dynamics:HALF-LIFE", "half_life", 1.0),
        case(
            "dynamics:HALF-LIFE",
            "decay_factor",
            half_life(TAU),
            TAU,
            note="halves after one half-life",
        ),
    ]


def pulses() -> list[dict]:
    peak, steady, monotone = (
        "dynamics:PEAK-AFTER-N",
        "dynamics:STEADY-STATE",
        "dynamics:SUMMATION-MONOTONE",
    )
    return [
        case(
            peak, "peak_after_pulses", 1.0, TEN_MIN, TAU, 1, note="the first peak is h"
        ),
        case(
            HUMAN,
            "peak_after_pulses",
            1.0,
            TEN_MIN,
            TAU,
            2,
            note="1 + e^-8: the second spike barely adds",
        ),
        case(
            HUMAN,
            "peak_after_pulses",
            1.0,
            TEN_MIN,
            TAU,
            12,
            note="12 spikes, as in the egg model",
        ),
        *(case(peak, "peak_after_pulses", 1.0, 37.5, TAU, n) for n in (2, 5, 10, 50)),
        case(
            peak,
            "peak_after_pulses",
            0.3,
            1.0,
            TAU,
            100,
            note="pulses much faster than the time constant",
        ),
        case(
            peak,
            "peak_after_pulses",
            1.0,
            1e-9,
            1.0,
            1000,
            note="almost no decay: close to n·h",
        ),
        case(
            peak,
            "peak_after_pulses",
            2.0,
            1e-30,
            1.0,
            7,
            note="no decay at all in f64: exactly 14",
        ),
        case(steady, "steady_state_peak", 1.0, TEN_MIN, TAU, note="1 / (1 - e^-8)"),
        case(steady, "steady_state_peak", 1.0, 37.5, TAU),
        case(steady, "steady_state_peak", 2.0, 1.0, 100.0),
        case(
            steady,
            "steady_state_peak",
            1.0,
            1e-6,
            1.0,
            note="period a millionth of tau: about 1e6",
        ),
        case(
            steady,
            "steady_state_peak",
            1e-20,
            1e-320,
            1.0,
            note="the summation ratio alone overflows; h/(1 - r) does not",
        ),
        case(
            "dynamics:STEADY-TROUGH",
            "steady_state_trough",
            1.0,
            TEN_MIN,
            TAU,
            note="e^-8 / (1 - e^-8)",
        ),
        case("dynamics:STEADY-TROUGH", "steady_state_trough", 1.0, 37.5, TAU),
        case("dynamics:STEADY-TROUGH", "steady_state_trough", 2.0, 1.0, 100.0),
        case(
            "dynamics:STEADY-TROUGH",
            "steady_state_trough",
            1e300,
            750.0,
            1.0,
            note="r = e^-750 alone is 0 in f64",
        ),
        case(
            "dynamics:STEADY-TROUGH",
            "steady_state_trough",
            1e-20,
            1e-320,
            1.0,
            note="about h·τ/T, as T/τ is tiny",
        ),
        case(
            "dynamics:STEADY-TROUGH",
            "steady_state_trough",
            2.3e-316,
            1e-8,
            1.0,
            note="a subnormal h, divided once: a normal result with all its digits",
        ),
        case(HUMAN, "summation_ratio", TEN_MIN, TAU, note="below 1.0004"),
        case(HUMAN, "summation_ratio", THIRTY_FIVE_MIN, TAU, note="1 + 7e-13"),
        case(
            monotone,
            "summation_ratio",
            TAU,
            TAU,
            note="one pulse per time constant: 1.58",
        ),
        case(
            monotone,
            "summation_ratio",
            7.5,
            TAU,
            note="ten pulses per time constant: 10.5",
        ),
        case(monotone, "summation_ratio", 1e-8, 1.0, note="about 1e8"),
        case(monotone, "summation_ratio", 1000.0, 1.0, note="1 + e^-1000 is 1 in f64"),
        case(
            HUMAN,
            "summation_ratio_at_hz",
            1 / TEN_MIN,
            TAU,
            note="1/600 Hz, through rhythm-units",
        ),
        case(HUMAN, "summation_ratio_at_hz", 1 / THIRTY_FIVE_MIN, TAU),
        case("dynamics:SUMMATION-AT-HZ", "summation_ratio_at_hz", 1.0, 0.5),
        case(
            "dynamics:SUMMATION-AT-HZ",
            "summation_ratio_at_hz",
            10.0,
            TAU,
            note="10 Hz pulses add up 750 times",
        ),
    ]


def thresholds() -> list[dict]:
    exact = "dynamics:PULSES-TO-FRACTION-EXACT"

    def boundary(period: float, tau: float, n: int) -> list[dict]:
        """p at, just below and just above 1 - r^n, the fraction reached exactly by the n-th peak."""
        on = float(-expm1(-n * time_constants(period, tau)))
        out = []
        for p, where in (
            (math.nextafter(on, 0.0), "just below"),
            (on, "the double nearest"),
            (math.nextafter(on, 1.0), "just above"),
        ):
            c = case(
                exact,
                "pulses_to_fraction",
                p,
                period,
                tau,
                note=f"p {where} 1 - r^{n} for T/τ = {period}/{tau}",
            )
            if naive_pulses(p, period, tau) != c["expected"]:
                c["note"] += (
                    f"; the floating-point formula gives {naive_pulses(p, period, tau)}"
                )
            out.append(c)
        return out

    tie = 2.0**-300
    return [
        case(
            exact,
            "pulses_to_fraction",
            0.9,
            37.5,
            TAU,
            note="pulses every τ/2 pass 90 % at the fifth",
        ),
        case(
            exact,
            "pulses_to_fraction",
            0.5,
            TEN_MIN,
            TAU,
            note="the first human spike already gives 99.97 %",
        ),
        case(exact, "pulses_to_fraction", 0.99, 1.0, TAU),
        case(
            exact,
            "pulses_to_fraction",
            0.632,
            TAU,
            TAU,
            note="1 - e^-1 = 0.63212 > 0.632",
        ),
        case(exact, "pulses_to_fraction", 0.999999, 1.0, 1000.0),
        case(
            exact,
            "pulses_to_fraction",
            0.5,
            1e-10,
            1.0,
            note="about 6.9 billion pulses",
        ),
        case(exact, "pulses_to_fraction", 1e-300, 1.0, 1.0),
        case(
            exact,
            "pulses_to_fraction",
            1e-300,
            1e-310,
            1.0,
            note="a subnormal period: about 1e10 pulses",
        ),
        case(
            exact,
            "pulses_to_fraction",
            tie,
            tie,
            1.0,
            note="an exact tie, T/τ = p: one pulse gives 1 - e^-p < p, so two are needed",
        ),
        *boundary(37.5, TAU, 5),
        *boundary(TEN_MIN, TAU, 1),
        *boundary(1.0, 3.0, 7),
        *boundary(0.1, 1.0, 20),
        *boundary(0.001, 1.0, 3),
        *boundary(1e-6, 1.0, 693147),
    ]


def spikes() -> list[dict]:
    peak = "dynamics:SPIKE-PEAK"
    human_peak, mouse_peak = spike_peak_time(15.0, 60.0), spike_peak_time(5.0, 25.0)
    for rise, decay, t in (
        (15.0, 60.0, human_peak),
        (5.0, 25.0, mouse_peak),
    ):  # the peak is a maximum
        top = Decimal(spike_shape(t, rise, decay))
        assert all(
            Decimal(spike_shape(t * f, rise, decay)) < top for f in (0.999, 1.001)
        )
    return [
        case(
            peak,
            "spike_peak_time",
            15.0,
            60.0,
            note="the egg model's human spike: 15·ln 5 = 24.14 s",
        ),
        case(
            peak,
            "spike_peak_time",
            5.0,
            25.0,
            note="the egg model's mouse spike: 5·ln 6 = 8.96 s",
        ),
        case(peak, "spike_peak_time", 1.0, 1.0, note="ln 2"),
        case(peak, "spike_peak_time", 60.0, 15.0),
        case(
            peak,
            "spike_peak_time",
            1e300,
            1e-10,
            note="decay/rise is subnormal: the peak time is decay",
        ),
        case(peak, "spike_peak_time", 1e-300, 1e300, note="decay/rise overflows"),
        case(AGREE, "spike_shape", 0.0, 15.0, 60.0, note="a spike starts at 0"),
        case(AGREE, "spike_shape", 15.0, 15.0, 60.0),
        case(
            peak,
            "spike_shape",
            human_peak,
            15.0,
            60.0,
            note="the top: 0.8·5^-0.25 = 0.535",
        ),
        case(AGREE, "spike_shape", 120.0, 15.0, 60.0),
        case(AGREE, "spike_shape", 1e4, 15.0, 60.0, note="long after: e^-166.7"),
        case(peak, "spike_shape", mouse_peak, 5.0, 25.0),
    ]


def errors() -> list[dict]:
    return [
        error("decay_factor", -1.0, TAU, note="dt must be >= 0"),
        error("decay_factor", 1.0, 0.0, note="tau must be > 0"),
        error("decay_factor", 1.0, -TAU, note="tau must be > 0"),
        error("relax", -0.5, 1.0, 1.0, note="a state must be >= 0"),
        error("relax", 1.0, -1.0, 1.0, note="dt must be >= 0"),
        error("half_life", 0.0, note="tau must be > 0"),
        error("half_life", -TAU, note="tau must be > 0"),
        error("peak_after_pulses", 0.0, TEN_MIN, TAU, 3, note="h must be > 0"),
        error("peak_after_pulses", 1.0, 0.0, TAU, 3, note="the period must be > 0"),
        error("peak_after_pulses", 1.0, TEN_MIN, 0.0, 3, note="tau must be > 0"),
        error(
            "peak_after_pulses",
            1.0,
            TEN_MIN,
            TAU,
            0,
            note="counting starts at the first pulse",
        ),
        error(
            "peak_after_pulses",
            1.0,
            TEN_MIN,
            TAU,
            -2,
            note="counting starts at the first pulse",
        ),
        error(
            "peak_after_pulses", 1e300, 1e-9, 1.0, 1000000000, note="the peak overflows"
        ),
        error("steady_state_peak", -1.0, 1.0, 1.0, note="h must be > 0"),
        error("steady_state_peak", 1e300, 1e-9, 1.0, note="the steady peak overflows"),
        error("steady_state_trough", 1.0, 0.0, 1.0, note="the period must be > 0"),
        error("summation_ratio", 0.0, TAU, note="the period must be > 0"),
        error("summation_ratio", TEN_MIN, -TAU, note="tau must be > 0"),
        error("summation_ratio", 1e-320, 1.0, note="1 / (1 - r) overflows"),
        error("summation_ratio_at_hz", 0.0, TAU, note="hz must be > 0"),
        error("summation_ratio_at_hz", -1.0, TAU, note="hz must be > 0"),
        error("summation_ratio_at_hz", 1e-310, TAU, note="the period 1/hz overflows"),
        error("summation_ratio_at_hz", 1.0, 0.0, note="tau must be > 0"),
        error("pulses_to_fraction", 0.0, 1.0, 1.0, note="p must be > 0"),
        error(
            "pulses_to_fraction",
            1.0,
            1.0,
            1.0,
            note="the steady state itself is never reached",
        ),
        error("pulses_to_fraction", -0.5, 1.0, 1.0, note="p must be in (0, 1)"),
        error("pulses_to_fraction", 1.5, 1.0, 1.0, note="p must be in (0, 1)"),
        error("pulses_to_fraction", 0.5, 0.0, 1.0, note="the period must be > 0"),
        error("pulses_to_fraction", 0.5, 1.0, -1.0, note="tau must be > 0"),
        error(
            "pulses_to_fraction",
            0.999999,
            1e-15,
            1.0,
            note="1.4e16 pulses: more than 2^53",
        ),
        error("pulses_to_fraction", 0.5, 1e-300, 1.0, note="far more than 2^53 pulses"),
        error("spike_shape", -1.0, 15.0, 60.0, note="s must be >= 0"),
        error("spike_shape", 1.0, 0.0, 60.0, note="rise must be > 0"),
        error("spike_shape", 1.0, 15.0, -60.0, note="decay must be > 0"),
        error("spike_peak_time", 0.0, 60.0, note="rise must be > 0"),
        error("spike_peak_time", 15.0, 0.0, note="decay must be > 0"),
    ]


if __name__ == "__main__":
    # exp(-x) turns the rounding of x into a relative error of about x·1e-16, so 1e-13 leaves a wide margin
    # up to x ≈ 170, the largest argument here. Pulse counts are whole numbers and must match exactly.
    groups = {
        "relaxation": (
            "Exponential relaxation: decay_factor, relax and half_life.",
            relaxation(),
            1e-13,
        ),
        "pulses": (
            "Periodic pulses: peaks, the steady state and the summation ratio.",
            pulses(),
            1e-13,
        ),
        "thresholds": (
            "pulses_to_fraction, exact at the boundaries.",
            thresholds(),
            0.0,
        ),
        "spikes": ("The spike shape and its peak time.", spikes(), 1e-13),
        "errors": (
            "Invalid input is rejected with an error in every binding.",
            errors(),
            1e-13,
        ),
    }
    for name, (description, cases, rel) in groups.items():
        write(name, description, cases, rel=rel)
        print(f"wrote {len(cases):3} cases to vectors/{name}.json")
