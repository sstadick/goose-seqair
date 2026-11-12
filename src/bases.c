#include <stdint.h>
#include "../vendor/goose-simd/include/goosimd_nibbles.h"

// Match Goose's packed byte-slice ABI, like goose-io's native adapters.
#pragma pack(push, 1)
typedef struct { uint8_t *data; int64_t len; } seqair_slice;
#pragma pack(pop)

// The Goose caller checks ceil(output.len / 2) == input.len. Input and output
// belong to different owners (BAM reader and region store). No overread/padding.
void goose_seqair_unpack_bases(seqair_slice input, seqair_slice output) {
    static const uint8_t identity[16] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    goosimd_decode_nibbles(input.data, output.data, (size_t)output.len, identity, 0);
}

void goose_seqair_sequence_text(seqair_slice input, seqair_slice output) {
    static const uint8_t alphabet[16] = "=ACMGRSVTWYHKDBN";
    goosimd_decode_nibbles(input.data, output.data, (size_t)output.len, alphabet, 0);
}

// Seqair's ASCII lookup semantics: printable unknowns become N, U becomes T,
// and 0..3 denote A/C/G/T. The high bit records invalid SAM sequence bytes.
static unsigned sam_code(uint8_t c) {
    static const uint8_t letters[26] = {
        1,14,2,13,15,15,4,11,15,15,12,15,3,15,15,15,15,5,6,8,8,7,9,15,10,15
    };
    unsigned upper = c >= 'a' && c <= 'z' ? c - 32u : c;
    unsigned code = upper >= 'A' && upper <= 'Z' ? letters[upper - 'A'] :
                    c == '=' ? 0 : c >= '0' && c <= '3' ? 1u << (c - '0') : 15;
    return code | ((c < 33 || c > 126) ? 256u : 0u);
}

// Caller supplies ceil(input.len/2) output bytes. No extra decoded sequence
// allocation and no padding reads, including for odd-length sequences.
int32_t goose_seqair_pack_sam(seqair_slice input, seqair_slice output) {
    unsigned invalid = 0;
    int64_t pairs = input.len / 2;
    for (int64_t i = 0; i < pairs; ++i) {
        unsigned a = sam_code(input.data[2*i]), b = sam_code(input.data[2*i+1]);
        invalid |= a | b;
        output.data[i] = (uint8_t)((a << 4) | b);
    }
    if (input.len & 1) {
        unsigned a = sam_code(input.data[input.len - 1]);
        invalid |= a;
        output.data[pairs] = (uint8_t)(a << 4);
    }
    return invalid < 16;
}

int32_t goose_seqair_pack_codes(seqair_slice input, seqair_slice output) {
    unsigned invalid = 0;
    int64_t pairs = input.len / 2;
    for (int64_t i = 0; i < pairs; ++i) {
        unsigned a = input.data[2*i], b = input.data[2*i+1];
        invalid |= a | b;
        output.data[i] = (uint8_t)((a << 4) | b);
    }
    if (input.len & 1) {
        unsigned a = input.data[input.len - 1];
        invalid |= a;
        output.data[pairs] = (uint8_t)(a << 4);
    }
    return invalid < 16;
}
