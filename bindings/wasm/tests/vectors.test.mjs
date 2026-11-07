// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama
// The shared golden vectors, checked against the WebAssembly build.

import assert from "node:assert/strict";
import { createRequire } from "node:module";
import { test } from "node:test";
import { checkDir } from "golden-vectors";

const impl = createRequire(import.meta.url)("../pkg/node/dynamics_core_wasm.js");

test("golden vectors", () => {
  assert.ok(checkDir(new URL("../../../vectors/", import.meta.url), impl) > 0);
});

test("pulse counts are whole numbers from 1 to 2^53, never wrapped or truncated", () => {
  assert.equal(impl.peak_after_pulses(1, 1, 1, 2), 1 + Math.exp(-1));
  const steady = impl.steady_state_peak(1, 1, 1);
  for (const n of [2 ** 32 + 2, 3e9, 2 ** 53]) assert.equal(impl.peak_after_pulses(1, 1, 1, n), steady);
  for (const n of [2.5, 0, -1, NaN, Infinity, 2 ** 53 + 2]) {
    assert.throws(() => impl.peak_after_pulses(1, 1, 1, n), /^Error: n must be /, `n = ${n}`);
  }
});
