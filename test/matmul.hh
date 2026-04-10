#pragma once

template <typename I, typename T>
void matmul(I m, I n, I k, T const * a, I lda, T const * b, I ldb, T * c, I ldc) {
  for (I j = 0; j < n; ++j) {
    for (I i = 0; i < m; ++i) {
      c[i + ldc*j] = 0;
      for (I l = 0; l < k; ++l) {
        c[i + ldc*j] += a[i + lda*k] * b[k + ldb*j];
      }
    }
  }
}
