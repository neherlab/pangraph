/* Copyright (c) 2018-2019 Evan Nemerson <evan@nemerson.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#if !defined(SIMDE__INSIDE_NEON_H)
#  error Do not include simde/arm/neon/float32x4.h directly; use simde/arm/neon.h.
#endif

#if !defined(SIMDE__NEON_FLOAT32X4_H)
#define SIMDE__NEON_FLOAT32X4_H

typedef union {
#if defined(SIMDE_VECTOR_SUBSCRIPT)
  simde_float32   f32 SIMDE_VECTOR(16) SIMDE_MAY_ALIAS;
#else
  simde_float32   f32[4];
#endif

#if defined(SIMDE_NEON_NATIVE)
  float32x4_t     n;
#endif

#if defined(SIMDE_NEON_MMX)
  __m64           mmx[2];
#endif
#if defined(SIMDE_NEON_SSE)
  __m128          sse;
#endif
} simde_float32x4_t;

#if defined(SIMDE_NEON_NATIVE)
HEDLEY_STATIC_ASSERT(sizeof(float32x4_t) == sizeof(simde_float32x4_t), "float32x4_t size doesn't match simde_float32x4_t size");
#endif
HEDLEY_STATIC_ASSERT(16 == sizeof(simde_float32x4_t), "simde_float32x4_t size incorrect");

SIMDE__FUNCTION_ATTRIBUTES
simde_float32x4_t
simde_vaddq_f32(simde_float32x4_t a, simde_float32x4_t b) {
  simde_float32x4_t r;
#if defined(SIMDE_NEON_NATIVE)
  r.n = vaddq_f32(a.n, b.n);
#else
  SIMDE__VECTORIZE
  for (size_t i = 0 ; i < (sizeof(r.f32) / sizeof(r.f32[0])) ; i++) {
    r.f32[i] = a.f32[i] + b.f32[i];
  }
#endif
  return r;
}

SIMDE__FUNCTION_ATTRIBUTES
simde_float32x4_t
simde_vld1q_f32 (simde_float32 const ptr[4]) {
  simde_float32x4_t r;
#if defined(SIMDE_NEON_NATIVE)
  r.n = vld1q_f32(ptr);
#else
  SIMDE__VECTORIZE
  for (size_t i = 0 ; i < (sizeof(r.f32) / sizeof(r.f32[0])) ; i++) {
    r.f32[i] = ptr[i];
  }
#endif
  return r;
}

SIMDE__FUNCTION_ATTRIBUTES
simde_float32x4_t
simde_x_vloadq_f32 (simde_float32 l0, simde_float32 l1, simde_float32 l2, simde_float32 l3) {
  simde_float32 v[] = { l0, l1, l2, l3 };
  return simde_vld1q_f32(v);
}

SIMDE__FUNCTION_ATTRIBUTES
simde_float32x4_t
simde_vdupq_n_f32 (simde_float32 value) {
  simde_float32x4_t r;
#if defined(SIMDE_NEON_NATIVE)
  r.n = vdupq_n_f32(value);
#else
  SIMDE__VECTORIZE
  for (size_t i = 0 ; i < (sizeof(r.f32) / sizeof(r.f32[0])) ; i++) {
    r.f32[i] = value;
  }
#endif
  return r;
}

SIMDE__FUNCTION_ATTRIBUTES
simde_float32x4_t
simde_vmulq_f32(simde_float32x4_t a, simde_float32x4_t b) {
  simde_float32x4_t r;
#if defined(SIMDE_NEON_NATIVE)
  r.n = vmulq_f32(a.n, b.n);
#else
  SIMDE__VECTORIZE
  for (size_t i = 0 ; i < (sizeof(r.f32) / sizeof(r.f32[0])) ; i++) {
    r.f32[i] = a.f32[i] * b.f32[i];
  }
#endif
  return r;
}

SIMDE__FUNCTION_ATTRIBUTES
simde_float32x4_t
simde_vsubq_f32(simde_float32x4_t a, simde_float32x4_t b) {
  simde_float32x4_t r;
#if defined(SIMDE_NEON_NATIVE)
  r.n = vsubq_f32(a.n, b.n);
#else
  SIMDE__VECTORIZE
  for (size_t i = 0 ; i < (sizeof(r.f32) / sizeof(r.f32[0])) ; i++) {
    r.f32[i] = a.f32[i] - b.f32[i];
  }
#endif
  return r;
}

#endif
