# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama
"""The shared golden vectors, checked against the Python bindings."""

from pathlib import Path

import dynamics_core
import pytest
from evidence_kernel.vectors import check_dir

VECTORS = Path(__file__).resolve().parents[3] / "vectors"


def test_golden_vectors():
    assert check_dir(VECTORS, dynamics_core) > 0


@pytest.mark.parametrize(
    ("n", "error"),
    [
        (0, ValueError),
        (2.5, TypeError),
        (2**63, OverflowError),
        (-(2**63) - 1, OverflowError),
    ],
)
def test_bad_pulse_counts_raise(n, error):
    with pytest.raises(error):
        dynamics_core.peak_after_pulses(1.0, 1.0, 1.0, n)
