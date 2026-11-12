#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <locale.h>
#include <pthread.h>

#pragma pack(push, 1)
typedef struct { uint8_t *data; int64_t len; } sam_slice;
#pragma pack(pop)

static pthread_once_t locale_once = PTHREAD_ONCE_INIT;
static locale_t numeric_locale;
static void make_locale(void) { numeric_locale = newlocale(LC_NUMERIC_MASK, "C", (locale_t)0); }

// Seqair sam/format.rs: scalar tags use six-digit %g; B:f uses kputd's
// scaled-double, ties-to-even rounding. No process-wide locale mutation.
int64_t goose_seqair_sam_float(double value, uint8_t array, sam_slice output) {
    pthread_once(&locale_once, make_locale);
    if (!numeric_locale || output.len < 32) return -1;
    locale_t old = uselocale(numeric_locale);
    double magnitude = fabs(value);
    int n;
    if (array && magnitude >= 0.0001 && magnitude <= 999999.0) {
        unsigned decimals = 0;
        double scale = 1.0;
        while (magnitude * scale < 100000.0 && decimals < 9) { scale *= 10.0; ++decimals; }
        double scaled = magnitude * scale, whole = floor(scaled), fraction = scaled - whole;
        uint32_t integer = (uint32_t)whole;
        if (fraction > 0.5 || (fraction == 0.5 && (integer & 1))) ++integer;
        char digits[16]; int count = snprintf(digits, sizeof digits, "%u", integer);
        n = 0;
        if (signbit(value)) output.data[n++] = '-';
        if ((unsigned)count <= decimals) {
            output.data[n++] = '0'; output.data[n++] = '.';
            for (unsigned j = (unsigned)count; j < decimals; ++j) output.data[n++] = '0';
            for (int j = 0; j < count; ++j) output.data[n++] = (uint8_t)digits[j];
        } else {
            for (int j = 0; j < count; ++j) {
                if (decimals && j == count - (int)decimals) output.data[n++] = '.';
                output.data[n++] = (uint8_t)digits[j];
            }
        }
        if (decimals) {
            while (n && output.data[n-1] == '0') --n;
            if (n && output.data[n-1] == '.') --n;
        }
    } else if (isnan(value)) {
        n = snprintf((char *)output.data, (size_t)output.len, "%snan", signbit(value) ? "-" : "");
    } else {
        n = snprintf((char *)output.data, (size_t)output.len, "%g", value);
    }
    uselocale(old);
    return n;
}
