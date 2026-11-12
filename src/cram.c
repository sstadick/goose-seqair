#include <stdint.h>
#include <string.h>
// Shared native primitives pinned to goose-simd v0.3.0; see vendor/goose-simd.
#include "../vendor/goose-simd/include/goosimd_expand.h"
#include "../vendor/goose-simd/include/goosimd_nibbles.h"
#ifndef _WIN32
#include <dlfcn.h>
#include <pthread.h>
#endif

#pragma pack(push, 1)
typedef struct { uint8_t *data; int64_t len; } cram_slice;
typedef struct { const uint32_t *data; int64_t len; } cram_words;
#pragma pack(pop)

// Seqair 8401961 rans_nx16.rs: packed scalar four-state decoding and SIMD
// thirty-two-state decoding, with checked scalar tails. Keep the portable
// path for unsupported targets; ISA kernels are selected once per block.
#if defined(__clang__) && !defined(GOOSIMD_SCALAR) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#define NX16_AVX2 1

// Seqair's RENORM_SHUFFLE8, generated as constant initializers: the next
// u16 belongs to each marked lane in lane order; other bytes are zeroed.
#define NX16_BYTE(m, l, b) ((m) & (1u << (l)) ? 2 * __builtin_popcount((m) & ((1u << (l)) - 1)) + (b) : 128)
#define NX16_LANE(m, l) NX16_BYTE(m, l, 0), NX16_BYTE(m, l, 1), 128, 128
#define NX16_ROW(m) { NX16_LANE(m, 0), NX16_LANE(m, 1), NX16_LANE(m, 2), NX16_LANE(m, 3), \
                      NX16_LANE(m, 4), NX16_LANE(m, 5), NX16_LANE(m, 6), NX16_LANE(m, 7) }
#define NX16_ROWS(m) NX16_ROW(m), NX16_ROW(m+1), NX16_ROW(m+2), NX16_ROW(m+3), \
                     NX16_ROW(m+4), NX16_ROW(m+5), NX16_ROW(m+6), NX16_ROW(m+7), \
                     NX16_ROW(m+8), NX16_ROW(m+9), NX16_ROW(m+10), NX16_ROW(m+11), \
                     NX16_ROW(m+12), NX16_ROW(m+13), NX16_ROW(m+14), NX16_ROW(m+15)
static const uint8_t nx16_shuffle[256][32] = {
    NX16_ROWS(0), NX16_ROWS(16), NX16_ROWS(32), NX16_ROWS(48),
    NX16_ROWS(64), NX16_ROWS(80), NX16_ROWS(96), NX16_ROWS(112),
    NX16_ROWS(128), NX16_ROWS(144), NX16_ROWS(160), NX16_ROWS(176),
    NX16_ROWS(192), NX16_ROWS(208), NX16_ROWS(224), NX16_ROWS(240),
};
#undef NX16_ROWS
#undef NX16_ROW
#undef NX16_LANE
#undef NX16_BYTE

__attribute__((target("sse2")))
static void nx16_flush_tile(const uint8_t tile[32][32], uint8_t *output, int64_t width, int64_t begin) {
    // Four 16x16 quarters; four rounds of byte interleaving transpose each.
    for (unsigned a = 0; a < 2; ++a) {
        for (unsigned b = 0; b < 2; ++b) {
            __m128i rows[16], next[16];
            for (unsigned k = 0; k < 16; ++k)
                rows[k] = _mm_loadu_si128((const __m128i *)&tile[a * 16 + k][b * 16]);
            for (unsigned round = 0; round < 4; ++round) {
                for (unsigned k = 0; k < 8; ++k) {
                    next[2*k] = _mm_unpacklo_epi8(rows[k], rows[k+8]);
                    next[2*k+1] = _mm_unpackhi_epi8(rows[k], rows[k+8]);
                }
                for (unsigned k = 0; k < 16; ++k) rows[k] = next[k];
            }
            for (unsigned k = 0; k < 16; ++k)
                _mm_storeu_si128((__m128i *)(output + (b * 16 + k) * width + begin + a * 16), rows[k]);
        }
    }
}

__attribute__((target("avx2")))
static int64_t nx16_vector(cram_slice input, cram_slice output, const uint32_t *restrict table,
                           unsigned bits, int order, uint32_t *restrict states, uint32_t *restrict contexts, int64_t *at) {
    uint8_t tile[32][32];
    __m256i slots[4];
    const uint8_t *restrict source = input.data;
    uint8_t *restrict destination = output.data;
    const __m256i mask = _mm256_set1_epi32((1u << bits) - 1);
    const __m256i bytes = _mm256_set1_epi32(order ? 255 : 0), zero = _mm256_setzero_si256();
    const __m128i shift = _mm_cvtsi32_si128(bits);
    const __m256i narrow = _mm256_broadcastsi128_si256(_mm_setr_epi8(0,4,8,12,-128,-128,-128,-128,-128,-128,-128,-128,-128,-128,-128,-128));
    int64_t width = output.len / 32, i = 0, position = *at;
    for (; i < width && input.len - position >= 80; ++i) {
        uint8_t *out = order ? tile[i & 31] : destination + i * 32;
        for (unsigned group = 0; group < 4; ++group) {
            __m256i x = _mm256_loadu_si256((const __m256i *)(states + 8 * group));
            __m256i ctx = _mm256_loadu_si256((const __m256i *)(contexts + 8 * group));
            __m256i index = _mm256_or_si256(ctx, _mm256_and_si256(x, mask));
            slots[group] = _mm256_i32gather_epi32((const int *)table, index, 4);
            __m256i sym = _mm256_shuffle_epi8(slots[group], narrow);
            __m128i eight = _mm_unpacklo_epi32(_mm256_castsi256_si128(sym), _mm256_extracti128_si256(sym, 1));
            _mm_storel_epi64((__m128i *)(out + 8 * group), eight);
        }
        for (unsigned group = 0; group < 4; ++group) {
            __m256i s = slots[group], x = _mm256_loadu_si256((const __m256i *)(states + 8 * group));
            __m256i frequency = _mm256_add_epi32(_mm256_srli_epi32(s, 20), _mm256_set1_epi32(1));
            x = _mm256_add_epi32(_mm256_mullo_epi32(frequency, _mm256_srl_epi32(x, shift)),
                                _mm256_and_si256(_mm256_srli_epi32(s, 8), _mm256_set1_epi32(4095)));
            // Unsigned x < 32768; signed comparison misclassifies large states.
            __m256i low = _mm256_cmpeq_epi32(_mm256_srli_epi32(x, 15), zero);
            unsigned m = (unsigned)_mm256_movemask_ps(_mm256_castsi256_ps(low));
            __m256i raw = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)(source + position)));
            __m256i norm = _mm256_shuffle_epi8(raw, _mm256_loadu_si256((const __m256i *)nx16_shuffle[m]));
            x = _mm256_blendv_epi8(x, _mm256_or_si256(_mm256_slli_epi32(x, 16), norm), low);
            position += 2 * __builtin_popcount(m);
            _mm256_storeu_si256((__m256i *)(states + 8 * group), x);
            __m256i ctx = _mm256_sll_epi32(_mm256_and_si256(s, bytes), shift);
            _mm256_storeu_si256((__m256i *)(contexts + 8 * group), ctx);
        }
        if (order && (i & 31) == 31) nx16_flush_tile(tile, destination, width, i - 31);
    }
    if (order) {
        for (unsigned row = 0; row < (unsigned)(i & 31); ++row)
            for (unsigned lane = 0; lane < 32; ++lane)
                destination[lane * width + i - (i & 31) + row] = tile[row][lane];
    }
    *at = position;
    return i;
}

// Seqair renorm16: widen sixteen little-endian u16s, expand into marked
// lanes with vpexpandd, then add to the selected shifted states.
__attribute__((target("avx512f,avx512bw,avx512dq,avx512vl")))
static int64_t nx16_avx512(cram_slice input, cram_slice output, const uint32_t *restrict table,
                          unsigned bits, int order, uint32_t *restrict states, uint32_t *restrict contexts, int64_t *at) {
    uint8_t tile[32][32];
    __m512i slots[2];
    __m512i mask = _mm512_set1_epi32((1u << bits) - 1), bytes = _mm512_set1_epi32(order ? 255 : 0);
    __m128i shift = _mm_cvtsi32_si128(bits);
    int64_t width = output.len / 32, i = 0, position = *at;
    for (; i < width && input.len - position >= 80; ++i) {
        uint8_t *out = order ? tile[i & 31] : output.data + i * 32;
        for (unsigned group = 0; group < 2; ++group) {
            __m512i x = _mm512_loadu_si512(states + 16 * group);
            __m512i ctx = _mm512_loadu_si512(contexts + 16 * group);
            __m512i index = _mm512_or_si512(ctx, _mm512_and_si512(x, mask));
            slots[group] = _mm512_i32gather_epi32(index, table, 4);
            _mm_storeu_si128((__m128i *)(out + 16 * group), _mm512_cvtepi32_epi8(slots[group]));
        }
        for (unsigned group = 0; group < 2; ++group) {
            __m512i s = slots[group], x = _mm512_loadu_si512(states + 16 * group);
            __m512i frequency = _mm512_add_epi32(_mm512_srli_epi32(s, 20), _mm512_set1_epi32(1));
            x = _mm512_add_epi32(_mm512_mullo_epi32(frequency, _mm512_srl_epi32(x, shift)),
                                 _mm512_and_si512(_mm512_srli_epi32(s, 8), _mm512_set1_epi32(4095)));
            __mmask16 low = _mm512_cmplt_epu32_mask(x, _mm512_set1_epi32(32768));
            __m512i words = _mm512_cvtepu16_epi32(_mm256_loadu_si256((const __m256i *)(input.data + position)));
            __m512i fresh = goosimd_expand_u32x16(words, low);
            x = _mm512_add_epi32(_mm512_mask_slli_epi32(x, low, x, 16), fresh);
            position += 2 * __builtin_popcount((unsigned)low);
            _mm512_storeu_si512(states + 16 * group, x);
            _mm512_storeu_si512(contexts + 16 * group, _mm512_sll_epi32(_mm512_and_si512(s, bytes), shift));
        }
        if (order && (i & 31) == 31) nx16_flush_tile(tile, output.data, width, i - 31);
    }
    if (order) {
        for (unsigned row = 0; row < (unsigned)(i & 31); ++row)
            for (unsigned lane = 0; lane < 32; ++lane)
                output.data[lane * width + i - (i & 31) + row] = tile[row][lane];
    }
    *at = position; return i;
}
#elif defined(__aarch64__) && !defined(GOOSIMD_SCALAR) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#include <arm_neon.h>
#define NX16_NEON 1

static void nx16_flush_tile(const uint8_t tile[32][32], uint8_t *output, int64_t width, int64_t begin) {
    for (unsigned a = 0; a < 2; ++a) for (unsigned b = 0; b < 2; ++b) {
        uint8x16_t rows[16], next[16];
        for (unsigned k = 0; k < 16; ++k) rows[k] = vld1q_u8(&tile[a * 16 + k][b * 16]);
        for (unsigned round = 0; round < 4; ++round) {
            for (unsigned k = 0; k < 8; ++k) {
                next[2*k] = vzip1q_u8(rows[k], rows[k+8]);
                next[2*k+1] = vzip2q_u8(rows[k], rows[k+8]);
            }
            for (unsigned k = 0; k < 16; ++k) rows[k] = next[k];
        }
        for (unsigned k = 0; k < 16; ++k) vst1q_u8(output + (b * 16 + k) * width + begin + a * 16, rows[k]);
    }
}

static int64_t nx16_neon(cram_slice input, cram_slice output, const uint32_t *restrict table,
                         unsigned bits, int order, uint32_t *restrict states, uint32_t *restrict contexts, int64_t *at) {
    static const uint32_t weights[4] = {1,2,4,8};
    static const uint8_t used[16] = {0,2,2,4,2,4,4,6,2,4,4,6,4,6,6,8};
    uint8_t tile[32][32];
    uint32_t next[32];
    uint32x4_t mask = vdupq_n_u32((1u << bits) - 1), bytes = vdupq_n_u32(order ? 255 : 0);
    int32x4_t right = vdupq_n_s32(-(int)bits), left = vdupq_n_s32((int)bits);
    int64_t width = output.len / 32, i = 0, position = *at;
    for (unsigned j = 0; j < 32; ++j) next[j] = contexts[j] | (states[j] & ((1u << bits) - 1));
    for (; i < width && input.len - position >= 80; ++i) {
        uint8_t *out = order ? tile[i & 31] : output.data + i * 32;
        // Sixteen states per batch, as upstream: keep the four vectors live
        // rather than spilling thirty-two gathered slots on a 4-lane target.
        for (unsigned batch = 0; batch < 32; batch += 16) {
            uint32_t slots[16];
            for (unsigned j = 0; j < 16; ++j) slots[j] = table[next[batch+j]];
            uint32x4_t sv[4];
            for (unsigned g = 0; g < 4; ++g) sv[g] = vld1q_u32(slots + 4*g);
            uint16x8_t lo = vcombine_u16(vmovn_u32(sv[0]), vmovn_u32(sv[1]));
            uint16x8_t hi = vcombine_u16(vmovn_u32(sv[2]), vmovn_u32(sv[3]));
            vst1q_u8(out + batch, vcombine_u8(vmovn_u16(lo), vmovn_u16(hi)));
            for (unsigned g = 0; g < 4; ++g) {
                unsigned lane = batch + 4*g;
                uint32x4_t s = sv[g], x = vld1q_u32(states + lane);
                uint32x4_t frequency = vaddq_u32(vshrq_n_u32(s, 20), vdupq_n_u32(1));
                x = vaddq_u32(vmulq_u32(frequency, vshlq_u32(x, right)), vandq_u32(vshrq_n_u32(s, 8), vdupq_n_u32(4095)));
                uint32x4_t low = vcltq_u32(x, vdupq_n_u32(32768));
                unsigned m = vaddvq_u32(vandq_u32(low, vld1q_u32(weights)));
                uint32x4_t words = vmovl_u16(vld1_u16((const uint16_t *)(input.data + position)));
                uint32x4_t fresh = goosimd_expand_u32x4(words, m);
                x = vaddq_u32(vbslq_u32(low, vshlq_n_u32(x, 16), x), fresh);
                position += used[m];
                uint32x4_t ctx = vshlq_u32(vandq_u32(s, bytes), left);
                vst1q_u32(states + lane, x); vst1q_u32(contexts + lane, ctx);
                vst1q_u32(next + lane, vorrq_u32(ctx, vandq_u32(x, mask)));
            }
        }
        if (order && (i & 31) == 31) nx16_flush_tile(tile, output.data, width, i - 31);
    }
    if (order) {
        for (unsigned row = 0; row < (unsigned)(i & 31); ++row)
            for (unsigned lane = 0; lane < 32; ++lane)
                output.data[lane * width + i - (i & 31) + row] = tile[row][lane];
    }
    *at = position; return i;
}
#endif

static uint32_t nx16_step(uint32_t state, uint32_t slot, unsigned bits) {
    return ((slot >> 20) + 1) * (state >> bits) + ((slot >> 8) & 4095);
}

static uint32_t nx16_renorm(uint32_t state, const uint8_t *input, int64_t *at) {
    uint32_t word = (uint32_t)input[*at] | (uint32_t)input[*at + 1] << 8;
    unsigned needs = state < 32768;
    *at += 2 * needs;
    return needs ? (state << 16) | word : state;
}

static int nx16_tail(uint32_t *state, cram_slice input, int64_t *at) {
    while (*state < 32768) {
        if (input.len - *at < 2) return 0;
        *state = (*state << 16) | (uint32_t)input.data[*at] |
                 (uint32_t)input.data[*at + 1] << 8;
        *at += 2;
    }
    return 1;
}

#if defined(__clang__)
__attribute__((always_inline))
#endif
static inline int64_t nx16_decode(cram_slice input, cram_slice output,
                                  const uint32_t *table, unsigned bits,
                                  unsigned lanes, int order) {
    if (input.len < 4 * (int64_t)lanes) return -1;
    uint32_t states[32], contexts[32] = {0};
    int normalized = 1;
    int64_t at = 0;
    for (unsigned lane = 0; lane < lanes; ++lane, at += 4) {
        states[lane] = (uint32_t)input.data[at] | (uint32_t)input.data[at+1] << 8 |
                       (uint32_t)input.data[at+2] << 16 | (uint32_t)input.data[at+3] << 24;
        normalized &= states[lane] >= 32768;
    }
    uint32_t mask = (1u << bits) - 1;
    int64_t width = output.len / lanes, i = 0;
    if (normalized) {
#ifdef NX16_AVX2
        if (lanes == 32 && __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512bw") &&
            __builtin_cpu_supports("avx512dq") && __builtin_cpu_supports("avx512vl"))
            i = nx16_avx512(input, output, table, bits, order, states, contexts, &at);
        else if (lanes == 32 && __builtin_cpu_supports("avx2"))
            i = nx16_vector(input, output, table, bits, order, states, contexts, &at);
#elif defined(NX16_NEON)
        if (lanes == 32) i = nx16_neon(input, output, table, bits, order, states, contexts, &at);
#endif
        for (; i < width && input.len - at >= 2 * (int64_t)lanes; ++i) {
            for (unsigned lane = 0; lane < lanes; ++lane) {
                uint32_t slot = table[contexts[lane] | (states[lane] & mask)];
                output.data[order ? lane * width + i : i * lanes + lane] = (uint8_t)slot;
                contexts[lane] = order ? (slot & 255) << bits : 0;
                states[lane] = nx16_renorm(nx16_step(states[lane], slot, bits), input.data, &at);
            }
        }
    }
    for (; i < width; ++i) {
        for (unsigned lane = 0; lane < lanes; ++lane) {
            uint32_t slot = table[contexts[lane] | (states[lane] & mask)];
            output.data[order ? lane * width + i : i * lanes + lane] = (uint8_t)slot;
            contexts[lane] = order ? (slot & 255) << bits : 0;
            states[lane] = nx16_step(states[lane], slot, bits);
            if (!nx16_tail(&states[lane], input, &at)) return -1;
        }
    }
    for (int64_t j = width * lanes; j < output.len; ++j) {
        unsigned lane = order ? lanes - 1 : (unsigned)(j - width * lanes);
        uint32_t slot = table[contexts[lane] | (states[lane] & mask)];
        output.data[j] = (uint8_t)slot;
        contexts[lane] = order ? (slot & 255) << bits : 0;
        states[lane] = nx16_step(states[lane], slot, bits);
        if (!nx16_tail(&states[lane], input, &at)) return -1;
    }
    return at;
}

// Keep the four-state registers separate from the address-taken 32-state
// SIMD scratch. Inlining both into the dispatcher spills the scalar loop.
#if defined(__clang__)
__attribute__((noinline))
#endif
static int64_t nx16_four(cram_slice input, cram_slice output, const uint32_t *table, unsigned bits, int order) {
    if (order) return nx16_decode(input, output, table, bits, 4, 1);
    return nx16_decode(input, output, table, bits, 4, 0);
}

int64_t goose_seqair_nx16_packed(cram_slice input, cram_slice output, cram_words table,
                                int32_t bits, int32_t lanes, int32_t order) {
    if (bits < 0 || bits > 12 || input.len < 0 || output.len < 0 ||
        (lanes != 4 && lanes != 32) || (order != 0 && order != 1) ||
        table.len < ((int64_t)(order ? 256 : 1) << bits)) return -1;
    // Constant lane counts/order let Clang specialize the scalar loops.
    if (lanes == 4) return nx16_four(input, output, table.data, bits, order);
    if (order) {
        return nx16_decode(input, output, table.data, bits, 32, 1);
    }
    return nx16_decode(input, output, table.data, bits, 32, 0);
}

void goose_seqair_unpack_nibbles(cram_slice input, cram_slice output, cram_slice mapping) {
    uint8_t table[16] = {0};
    memcpy(table, mapping.data, (size_t)mapping.len);
    // PACK leaves uncovered bytes zero and out-of-map indices zero.
    int64_t covered = output.len;
    if (input.len < covered / 2 + covered % 2) covered = input.len * 2;
    goosimd_decode_nibbles(input.data, output.data, (size_t)covered, table, 1);
}

// RFC 1321 MD5, used only for CRAM reference identity, not security. Process
// complete blocks without copying; the final one or two blocks are local.
static void md5_block(uint32_t h[4], const uint8_t *p) {
    static const uint32_t k[64] = {
        0xd76aa478,0xe8c7b756,0x242070db,0xc1bdceee,0xf57c0faf,0x4787c62a,0xa8304613,0xfd469501,
        0x698098d8,0x8b44f7af,0xffff5bb1,0x895cd7be,0x6b901122,0xfd987193,0xa679438e,0x49b40821,
        0xf61e2562,0xc040b340,0x265e5a51,0xe9b6c7aa,0xd62f105d,0x02441453,0xd8a1e681,0xe7d3fbc8,
        0x21e1cde6,0xc33707d6,0xf4d50d87,0x455a14ed,0xa9e3e905,0xfcefa3f8,0x676f02d9,0x8d2a4c8a,
        0xfffa3942,0x8771f681,0x6d9d6122,0xfde5380c,0xa4beea44,0x4bdecfa9,0xf6bb4b60,0xbebfbc70,
        0x289b7ec6,0xeaa127fa,0xd4ef3085,0x04881d05,0xd9d4d039,0xe6db99e5,0x1fa27cf8,0xc4ac5665,
        0xf4292244,0x432aff97,0xab9423a7,0xfc93a039,0x655b59c3,0x8f0ccc92,0xffeff47d,0x85845dd1,
        0x6fa87e4f,0xfe2ce6e0,0xa3014314,0x4e0811a1,0xf7537e82,0xbd3af235,0x2ad7d2bb,0xeb86d391,
    };
    static const unsigned shifts[16] = {7,12,17,22,5,9,14,20,4,11,16,23,6,10,15,21};
    uint32_t words[16];
    for (unsigned i = 0; i < 16; ++i) {
        words[i] = (uint32_t)p[4*i] | (uint32_t)p[4*i+1] << 8 |
                   (uint32_t)p[4*i+2] << 16 | (uint32_t)p[4*i+3] << 24;
    }
    uint32_t a = h[0], b = h[1], c = h[2], d = h[3];
    for (unsigned i = 0; i < 64; ++i) {
        uint32_t f; unsigned g;
        if (i < 16) { f = (b & c) | (~b & d); g = i; }
        else if (i < 32) { f = (d & b) | (~d & c); g = (5*i+1) & 15; }
        else if (i < 48) { f = b ^ c ^ d; g = (3*i+5) & 15; }
        else { f = c ^ (b | ~d); g = (7*i) & 15; }
        uint32_t x = a + f + k[i] + words[g];
        unsigned shift = shifts[(i >> 4)*4 + (i & 3)];
        a = d; d = c; c = b; b += (x << shift) | (x >> (32-shift));
    }
    h[0] += a; h[1] += b; h[2] += c; h[3] += d;
}

void goose_seqair_reference_md5(cram_slice input, cram_slice output) {
    uint32_t h[4] = {0x67452301,0xefcdab89,0x98badcfe,0x10325476};
    int64_t at = 0;
    while (input.len - at >= 64) { md5_block(h, input.data + at); at += 64; }
    uint8_t tail[128] = {0};
    unsigned left = (unsigned)(input.len - at);
    if (left) memcpy(tail, input.data + at, left);
    tail[left] = 0x80;
    unsigned last = left < 56 ? 64 : 128;
    uint64_t bits = (uint64_t)input.len * 8;
    for (unsigned i = 0; i < 8; ++i) tail[last-8+i] = (uint8_t)(bits >> (8*i));
    md5_block(h, tail);
    if (last == 128) md5_block(h, tail + 64);
    for (unsigned i = 0; i < 16; ++i) output.data[i] = (uint8_t)(h[i/4] >> (8*(i&3)));
}

// Optional standard compression backends. Loading is thread-safe and happens
// only on first use; BAM/SAM and the native Goose CRAM codecs need neither.
// These are the stable libbz2/liblzma C ABIs, not htslib/htscodecs adapters.
#ifndef _WIN32
typedef int (*bz_decode_fn)(char *, unsigned *, char *, unsigned, int, int);
typedef int (*xz_decode_fn)(uint64_t *, uint32_t, const void *, const uint8_t *,
                           size_t *, size_t, uint8_t *, size_t *, size_t);
static bz_decode_fn bz_decode;
static xz_decode_fn xz_decode;
static pthread_once_t compression_once = PTHREAD_ONCE_INIT;

static void load_compression(void) {
#ifdef __APPLE__
    void *bz = dlopen("libbz2.dylib", RTLD_LOCAL | RTLD_LAZY);
    void *xz = dlopen("liblzma.dylib", RTLD_LOCAL | RTLD_LAZY);
#else
    void *bz = dlopen("libbz2.so.1.0", RTLD_LOCAL | RTLD_LAZY);
    void *xz = dlopen("liblzma.so.5", RTLD_LOCAL | RTLD_LAZY);
#endif
    if (bz) bz_decode = (bz_decode_fn)dlsym(bz, "BZ2_bzBuffToBuffDecompress");
    if (xz) xz_decode = (xz_decode_fn)dlsym(xz, "lzma_stream_buffer_decode");
}
#endif

// 1 success, 0 corrupt/oversized data, -1 unavailable codec. Output size is
// exact, including zero bytes. XZ's internal dictionary allocation is bounded.
int32_t goose_seqair_cram_external(int32_t method, cram_slice input,
                                 cram_slice output, int64_t limit) {
#ifdef _WIN32
    (void)method; (void)input; (void)output; (void)limit;
    return -1;
#else
    if (input.len < 0 || output.len < 0 || output.len > limit || limit < 0) return 0;
    pthread_once(&compression_once, load_compression);
    uint8_t empty = 0;
    uint8_t *out = output.len ? output.data : &empty;
    if (method == 2) {
        if (!bz_decode) return -1;
        if ((uint64_t)input.len > UINT32_MAX || (uint64_t)output.len > UINT32_MAX) return 0;
        unsigned size = (unsigned)output.len;
        int result = bz_decode((char *)out, &size, (char *)input.data, (unsigned)input.len, 0, 0);
        return result == 0 && size == (unsigned)output.len;
    }
    if (method == 3) {
        if (!xz_decode) return -1;
        size_t read = 0, written = 0;
        uint64_t memory = (uint64_t)limit;
        int result = xz_decode(&memory, 0, NULL, input.data, &read, (size_t)input.len,
                               out, &written, (size_t)output.len);
        return result == 0 && written == (size_t)output.len;
    }
    return -1;
#endif
}
