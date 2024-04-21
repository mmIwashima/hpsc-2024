#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <x86intrin.h>

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }
  for(int i=0; i<N; i++) {
//    for(int j=0; j<N; j++) {
//      if(i != j) {
//        float rx = x[i] - x[j];
//        float ry = y[i] - y[j];
//        float r = std::sqrt(rx * rx + ry * ry);
//        fx[i] -= rx * m[j] / (r * r * r);
//        fy[i] -= ry * m[j] / (r * r * r);
//      }

    __m512 xvec = _mm512_load_ps(x);
    __m512 yvec = _mm512_load_ps(y);
    __m512 mvec = _mm512_load_ps(m);

    __m512 xivec = _mm512_set1_ps(x[i]);
    __m512 yivec = _mm512_set1_ps(y[i]);

    __m512 rxvec = _mm512_sub_ps(xvec, xivec);
    __m512 ryvec = _mm512_sub_ps(yvec, yivec);

    __m512 rxvec2 = _mm512_mul_ps(rxvec, rxvec);
    __m512 ryvec2 = _mm512_mul_ps(ryvec, ryvec);
    __m512 rvec2 = _mm512_add_ps(rxvec2, ryvec2);

    __m512 r_rev_vec = _mm512_rsqrt14_ps(rvec2);
    __m512 r2_rev_vec = _mm512_mul_ps(r_rev_vec, r_rev_vec);
    __m512 r3_rev_vec = _mm512_mul_ps(r2_rev_vec, r_rev_vec);

    __m512 fx_sub_vec = _mm512_mul_ps(rxvec, r3_rev_vec);
    __m512 fx_vec = _mm512_mul_ps(fx_sub_vec, mvec);
    
    __m512 fy_sub_vec = _mm512_mul_ps(ryvec, r3_rev_vec);
    __m512 fy_vec = _mm512_mul_ps(fy_sub_vec, mvec);

    __m512 fx_masked_vec = _mm512_setzero_ps();
    __m512 fy_masked_vec = _mm512_setzero_ps();
    __m512 threshold = _mm512_set1_ps(10e+5);
    __mmask16 mask_x = _mm512_cmp_ps_mask(fx_vec, threshold, _MM_CMPINT_GT);
    fx_masked_vec = _mm512_mask_blend_ps(mask_x, fx_vec, fx_masked_vec);
    __mmask16 mask_y = _mm512_cmp_ps_mask(fy_vec, threshold, _MM_CMPINT_GT);
    fy_masked_vec = _mm512_mask_blend_ps(mask_y, fy_vec, fy_masked_vec);

    fx[i] = _mm512_reduce_add_ps(fx_masked_vec);
    fy[i] = _mm512_reduce_add_ps(fy_masked_vec);

//    }
    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}
