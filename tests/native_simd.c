// Native decoder boundary/ISA check. Expected symbols and consumption come
// from the scalar rANS equations, not from another production code path.
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#pragma pack(push, 1)
typedef struct { uint8_t *data; int64_t len; } Bytes;
typedef struct { const uint32_t *data; int64_t len; } Words;
#pragma pack(pop)
extern int64_t goose_seqair_nx16_packed(Bytes, Bytes, Words, int32_t, int32_t, int32_t);

static uint32_t random_word(uint32_t *state) {
    *state ^= *state << 13; *state ^= *state >> 17; *state ^= *state << 5; return *state;
}
static int64_t oracle(Bytes input, Bytes output, Words table, unsigned bits, unsigned lanes, int order) {
    if (input.len < 4 * (int64_t)lanes) return -1;
    uint32_t states[32], context[32] = {0};
    int64_t at = 0, width = output.len / lanes;
    for (unsigned lane = 0; lane < lanes; ++lane, at += 4)
        states[lane] = (uint32_t)input.data[at] | (uint32_t)input.data[at+1] << 8 |
                       (uint32_t)input.data[at+2] << 16 | (uint32_t)input.data[at+3] << 24;
    for (int64_t step = 0; step < output.len; ++step) {
        unsigned lane = step < width * lanes ? step % lanes : order ? lanes - 1 : step - width * lanes;
        uint32_t slot = table.data[(context[lane] << bits) + (states[lane] % (1u << bits))];
        int64_t pos = order && step < width * lanes ? lane * width + step / lanes : step;
        output.data[pos] = slot & 255;
        context[lane] = order ? slot & 255 : 0;
        uint64_t state = ((uint64_t)(slot >> 20) + 1) * (states[lane] / (1u << bits)) + ((slot >> 8) & 4095);
        states[lane] = (uint32_t)state;
        while (states[lane] < 32768) {
            if (input.len - at < 2) return -1;
            states[lane] = states[lane] * 65536u + input.data[at] + 256u * input.data[at+1]; at += 2;
        }
    }
    return at;
}

static const char *isa(void) {
#ifdef GOOSIMD_SCALAR
    return "scalar";
#elif defined(__aarch64__)
    return "neon";
#elif defined(__x86_64__) || defined(__i386__)
    if (__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512bw") &&
        __builtin_cpu_supports("avx512dq") && __builtin_cpu_supports("avx512vl")) return "avx512";
    if (__builtin_cpu_supports("avx2")) return "avx2";
    return "scalar";
#else
    return "scalar";
#endif
}

int main(int argc, char **argv) {
    const char *selected = isa();
    if (argc == 2 && strcmp(argv[1], "auto") && strcmp(argv[1], selected)) {
        fprintf(stderr, "required %s, selected %s; no ISA coverage claimed\n", argv[1], selected); return 2;
    }
    const unsigned lengths[] = {0,1,3,31,32,33,991,992,993,1023,1024,1025,2017,2048,2055,4097};
    const unsigned precisions[] = {0,1,5,8,9,10,11,12};
    uint32_t seed = 0x81726354;
    size_t cases = 0;
    for (unsigned p = 0; p < sizeof precisions / sizeof *precisions; ++p) {
        unsigned bits = precisions[p], width = 1u << bits;
        uint32_t *table = malloc(256 * width * sizeof *table); assert(table);
        for (unsigned row = 0; row < 256; ++row) {
            unsigned begin = 0, symbol = row;
            while (begin < width) {
                unsigned remaining = width - begin;
                unsigned frequency = remaining < 31 ? remaining : 1 + random_word(&seed) % 31;
                for (unsigned bias = 0; bias < frequency; ++bias)
                    table[row * width + begin + bias] = ((frequency - 1) << 20) | (bias << 8) | (symbol & 255);
                begin += frequency; symbol += 53;
            }
        }
        for (int order = 0; order < 2; ++order) for (unsigned lanes = 4; lanes <= 32; lanes *= 8)
        for (unsigned l = 0; l < sizeof lengths / sizeof *lengths; ++l) for (unsigned variant = 0; variant < 4; ++variant) {
            size_t length = lengths[l], available = 4 * lanes + 2 * length + 80;
            uint8_t *bytes = malloc(available + 1), *expected = malloc(length + 2), *actual = malloc(length + 2);
            assert(bytes && expected && actual);
            for (size_t j = 0; j <= available; ++j) bytes[j] = (uint8_t)random_word(&seed);
            // Start unaligned. Normalized states exercise SIMD; small/high-bit
            // states catch signed comparisons and the scalar safety fallback.
            for (unsigned lane = 0; lane < lanes; ++lane) {
                uint32_t x = random_word(&seed) | 32768;
                if (variant == 1) x = lane & 1 ? 0x80000001 : 32768;
                if (variant == 2 && lane == 7 % lanes) x = 7;
                for (unsigned b = 0; b < 4; ++b) bytes[1 + 4*lane + b] = (uint8_t)(x >> (8*b));
            }
            if (variant == 3) {
                available = 4 * lanes + (length % 81);
                bytes = realloc(bytes, available + 1); assert(bytes);
            }
            memset(actual, 0xA5, length + 2); memset(expected, 0xA5, length + 2);
            Bytes input = {bytes + 1, (int64_t)available}, a = {actual + 1, (int64_t)length}, e = {expected + 1, (int64_t)length};
            Words t = {table, (int64_t)(order ? 256 : 1) * width};
            int64_t want = oracle(input, e, t, bits, lanes, order);
            int64_t got = goose_seqair_nx16_packed(input, a, t, bits, lanes, order);
            assert(got == want && !memcmp(actual, expected, length + 2));
            ++cases; free(bytes); free(actual); free(expected);
        }
        free(table);
    }
    printf("PASS %s: %zu rANS decode/consumption/tail cases\n", selected, cases);
    return 0;
}
