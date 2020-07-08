#include <stdio.h>
#include <stdlib.h>

/**
  * (c) 2020 Kryptonite
  *
  * The code implements the results of https://arxiv.org/pdf/2006.10286.pdf
  */
typedef unsigned long long Index;

Index ** corner_indexes = NULL;

Index parity(Index val) {
    return (val == 0) ? 0 : (val & 1) ^ (val >> 1);
}

Index grey(int const d, Index const val) {
    Index const val2 = val % (1 << d);
    return val2 ^ (val2 >> 1);
}

// May be avoided totally by caching of `grey` values
Index grey_inverse(int const d, Index const val) {
    Index const val2 = val % (1 << d);
    return (val2 == 0 || d < 0) ? 0 : (val2 ^ grey_inverse(d - 1, val2 >> 1));
}

// just stored values
// May be used after the first layer of corner_indexes is computed
Index grey_fast(int const d, Index const r) {
    return corner_indexes[0][r + (1 << d)];
}
// just stored values
// May be used after the first layer of corner_indexes is computed
Index grey_inverse_fast(int const d, Index const r) {
    return corner_indexes[0][r];
}

// Correspondence between usual and Z-order coordinates
void Z_transpose(int const d, int const n, Index const * const coords, Index * const alphas) {
    for (int i = 0; i != d; ++i) {
        alphas[i] = 0;
    }
    for (int i = 0; i != n; ++i) {
        for (int j = 0; j != d; ++j) {
            alphas[j] |= (coords[i] >> j & 1) << i;
        }
    }
}

/**
  * Main functions decode and encode
  * Since encode should return an array, it takes a pointer and writes the results there.
  */
Index encode_h(int const d, int const n, Index const * const alphas) {
    if (n <= 0) return 0;
    Index r = 0; // d*n-bit
    Index r_shift = 0; // d*(n-1)-bit
    Index alpha; // d-bit
    Index r0; // d-bit
    Index const two_power_d = 1 << d; // (d+1)-bit
    Index sub_cell_size; // (d*n + 1)-nit
    Index alpha_inv; // d-bit
    Index need_to_change_last; // 1-bit
    for (int i = n-1; i != -1; --i) {
        alpha = alphas[i];
        r0 = grey_inverse_fast(d, alpha); // contribution of highest bits
        sub_cell_size = 1 << (d* (n - 1 - i));
        alpha_inv = alpha ^ (two_power_d - 1);
        need_to_change_last = 1 ^ parity(alpha_inv); // actually, bool
        r_shift = corner_indexes[n - 1 - i][alpha_inv + two_power_d * need_to_change_last];
        if (d % 2 == 1 && n == 1) { // condition of reversal
            r_shift ^= two_power_d - 1;
            r_shift++;
        }
        r -= r_shift;
        r %= sub_cell_size;
        r += r0 * sub_cell_size;
    }
    return r;
}

void decode_h(int const d, int const n, Index /*const*/ r, Index * const alphas) {
    Index r0; // dn-bit
    Index alpha, alpha_inv; // d-bit
    Index r_shift; // d(n-i)-bit
    Index const two_power_d = 1 << d;
    Index need_to_change_last; // actually, bool
    for (int i = 0; i != n; ++i) {
        r0 = (r >> ((n - 1 - i) * d)) % two_power_d;
        alpha = grey_fast(d, r0);
        alphas[i] = alpha; // here we "return" value
        alpha_inv = alpha ^ (two_power_d - 1);
        need_to_change_last = 1 ^ parity(alpha_inv);
        r_shift = corner_indexes[n - 1 - i][alpha_inv + two_power_d * need_to_change_last];
        if (d % 2 == 1 && n == 1) { // condition of reversal
            r_shift ^= two_power_d - 1;
            r_shift++;
        }
        r += r_shift;
    }
}

/**
  * Preprocessing and evaluating angle indexes to make their evaluation just an access operation by pointer.
  * Also allocation operations.
  */

void evaluate_corner_indexes(int d, int n, Index *** angle_Is) {
    /**
      * Allocation of precomputed corner indexes
      */
    *angle_Is = (Index **)malloc((n + 1) * sizeof(Index *));
    size_t len = sizeof(Index) * (1 << (d + 1));
    Index const two_power_d = 1 << d;
    for (int n2 = 0; n2 != n + 1; ++n2) {
        (*angle_Is)[n2] = (Index *)malloc(len);
    }

    /**
      * evaluation
      * case d=0 separately: store values of g and g^{-1} to avoid its recomputation
      Actually, this is the only place, where 'grey' functions is needed.
      After `grey_fast` and `grey_inverse_fast` may be used, because they need only data access without any computations.
      */
    for (Index i = 0; i != 1 << d; ++ i) {
        Index g = grey(d, i);
        (*angle_Is)[0][g + two_power_d] = i;
        (*angle_Is)[0][i] = g;
    }

    /**
      * evaluation for increasing depth
      * Here the allocation may be flattened.
      * Not all the `alphas` need to be n*d-bit size, actual length linearly depends on n.
      */
    Index * alphas = (Index *)malloc(sizeof(Index) * n);
    for (int n1 = 1; n1 != n + 1; ++n1) {
        for (Index r = 0; r != two_power_d; ++r) {
            for (int n2 = 0; n2 < n; ++n2) {
                alphas[n2] = r;
            }
            (*angle_Is)[n1][r] = encode_h(d, n1, &alphas[n - n1]);
            alphas[n - 1] ^= 1;
            (*angle_Is)[n1][r + two_power_d] = encode_h(d, n1, &alphas[n - n1]);
        }
    }
    free(alphas);
}

void free_corner_indexes(int d, int n, Index *** coordinates) {
    for (int n2 = 0; n2 != n; ++n2) {
        free((*coordinates)[n2]);
    }
    free(*coordinates);
}
