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
#  error Do not include simde/arm/neon/uint16x4.h directly; use simde/arm/neon.h.
#endif

#if !defined(SIMDE__NEON_UINT16X4_H)
#define SIMDE__NEON_UINT16X4_H

typedef union {
#if defined(SIMDE_VECTOR_SUBSCRIPT)
  uint16_t         u16 SIMDE_VECTOR(8) SIMDE_MAY_ALIAS;
#else
  uint16_t         u16[4];
#endif

#if defined(SIMDE_NEON_NATIVE)
  uint16x4_t       n;
#endif

#if defined(SIMDE_NEON_MMX)
  __m64           mmx;
#endif
} simde_uint16x4_t;

#if defined(SIMDE_NEON_NATIVE)
HEDLEY_STATIC_ASSERT(sizeof(uint16x4_t) == sizeof(simde_uint16x4_t), "uint16x4_t size doesn't match simde_uint16x4_t size");
#endif
HEDLEY_STATIC_ASSERT(8 == sizeof(simde_uint16x4_t), "simde_uint16x4_t size incorrect");

SIMDE__FUNCTION_ATTRIBUTES
simde_uint16x4_t
simde_vadd_u16(simde_uint16x4_t a, simde_uint16x4_t b) {
  simde_uint16x4_t r;
#if defined(SIMDE_NEON_NATIVE)
  r.n = vadd_u16(a.n, b.n);
#else
  SIMDE__VECTORIZE
  for (size_t i = 0 ; i < (sizeof(r.u16) / sizeof(r.u16[0])) ; i++) {
    r.u16[i] = a.u16[i] + b.u16[i];
  }
#endif
  return r;
}

SIMDE__FUNCTION_ATTRIBUTES
simde_uint16x4_t
simde_vld1_u16 (uint16_t const ptr[4]) {
  simde_uint16x4_t r;
#if defined(SIMDE_NEON_NATIVE)
  r.n = vld1_u16(ptr);
#else
  SIMDE__VECTORIZE
  for (size_t i = 0 ; i < (sizeof(r.u16) / sizeof(r.u16[0])) ; i++) {
    r.u16[i] = ptr[i];
  }
#endif
  return r;
}

SIMDE__FUNCTION_ATTRIBUTES
simde_uint16x4_t
simde_x_vload_u16 (uint16_t l0, uint16_t l1, uint16_t l2, uint16_t l3) {
  uint16_t v[] = { l0, l1, l2, l3 };
  return simde_vld1_u16(v);
}

SIMDE__FUNCTION_ATTRIBUTES
simde_uint16x4_t
simde_vdup_n_u16 (uint16_t value) {
  simde_uint16x4_t r;
#if defined(SIMDE_NEON_NATIVE)
  r.n = vdup_n_u16(value);
#else
  SIMDE__VECTORIZE
  for (size_t i = 0 ; i < (sizeof(r.u16) / sizeof(r.u16[0])) ; i++) {
    r.u16[i] = value;
  }
#endif
  return r;
}

SIMDE__FUNCTION_ATTRIBUTES
simde_uint16x4_t
simde_vmul_u16(simde_uint16x4_t a, simde_uint16x4_t b) {
  simde_uint16x4_t r;
#if defined(SIMDE_NEON_NATIVE)
  r.n = vmul_u16(a.n, b.n);
#else
  SIMDE__VECTORIZE
  for (size_t i = 0 ; i < (sizeof(r.u16) / sizeof(r.u16[0])) ; i++) {
    r.u16[i] = a.u16[i] * b.u16[i];
  }
#endif
  return r;
}

SIMDE__FUNCTION_ATTRIBUTES
simde_uint16x4_t
simde_vsub_u16(simde_uint16x4_t a, simde_uint16x4_t b) {
  simde_uint16x4_t r;
#if defined(SIMDE_NEON_NATIVE)
  r.n = vsub_u16(a.n, b.n);
#else
  SIMDE__VECTORIZE
  for (size_t i = 0 ; i < (sizeof(r.u16) / sizeof(r.u16[0])) ; i++) {
    r.u16[i] = a.u16[i] - b.u16[i];
  }
#endif
  return r;
}

#endif
