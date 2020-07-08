/*
 * (c) C program code 2020 Kryptonite
 * Other license rights inherited from https://github.com/galtay/hilbertcurve/blob/develop/LICENSE
 */

typedef unsigned long long Index;

Index Z_transposed_to_index(int const d, int const n, Index * const alphas) {
    Index r = 0;
    for (int i = 0; i != n; ++i) {
        for (int j = 0; j != d; ++j) {
            r |= ((alphas[i] >> j) & 1) << (i + j*n);
        }
    }
    return r;
}

void index_to_Z_transposed(int const d, int const n, Index const r, Index * const alphas) {
    for (int i = 0; i != d; ++i) {
        alphas[i] = 0;
    }
    for (int i = 0; i != n; ++i) {
        for (int j = 0; j != d; ++j) {
            alphas[j] |= r >> (i + j*n) & 1 << i;
        }
    }
}

// distance_from_coordinates
Index encode_hilbert(int const d, int const n, Index * const alphas) {
    Index M = 1 << (n - 1);
    // inverse undo excess work
    Index Q = M;
    Index P;
    while (Q > 1) {
        P = Q - 1;
        for (int i = 0; i != n; ++i) {
            if (alphas[i] & Q) {
                alphas[0] ^= P;
            } else {
                Index t = (alphas[0] ^ alphas[i]) & P;
                alphas[0] ^= t;
                alphas[i] ^= t;
            }
        }
        Q >>= 1;
    }
    // Grey encode
    for (int i = 1; i != n; ++i) {
        alphas[i] ^= alphas[i - 1];
    }
    Index t = 0;
    Q = M;
    while (Q > 1) {
        if (alphas[d - 1] & Q) {
            t ^= Q - 1;
        }
        Q >>= 1;
    }
    for (int i = 0; i != n; ++i) {
        alphas[i] ^= t;
    }
    Index h = Z_transposed_to_index(d, n, alphas);
    return h;
}

// coordinates_from_distance
void decode_hilbert(int const d, int const n, Index r, Index * const alphas) {
    Index Z = 1 << n;
    index_to_Z_transposed(d, n, r, alphas);
    // Grey decode
    Index t = alphas[d-1] >> 1;
    for (int i = d - 1; i != -1; --i) {
        alphas[i] ^= alphas[i - 1];
    }
    alphas[0] ^= t;
    // Undo excess work
    Index Q = 2;
    Index P;
    while (Q != Z) {
        P = Q - 1;
        for (int i = d-1; i != -1; --i) {
            if ((alphas[i] & Q) != 0) {
                alphas[0] ^= P; // invert
            } else { // exchange
                Index t = (alphas[0] ^ alphas[i]) & P;
                alphas[0] ^= t;
                alphas[i] ^= t;
            }
        }
        Q <<= 1;
    }
    // alphas is what we "return"
}
