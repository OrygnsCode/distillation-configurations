// brute force: try every set of submixtures and keep the ones that are basic configurations.
// also counts the thermally coupled ones.
//
// components 1..n, most volatile first. a stream is an interval [i,j]. the feed [1,n] and the pure
// products [c,c] are always there. a stream [i,j] with i<j goes into a column: its top product is the
// next stream to the left in row i, its bottom product the next one down in column j.
//   - nothing can get lost in a split: the bottom product starts at most one past where the top ends
//   - every stream except the feed has to come out of some column, so there is a longer stream that
//     starts where it starts or ends where it ends
// a submixture that comes out of exactly one column has its own condenser or reboiler, and that one
// can be swapped for a thermal coupling. m counts them, so each config gives 2^m - 1 coupled ones.
//
// gcc -O3 -fopenmp -o tc_bf tc_bf.c
// ./tc_bf 8        (n up to 9; 9 is 2^35 sets and takes a while)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 10
static int n, nfree, fi[64], fj[64];

int main(int argc, char **argv) {
    n = atoi(argv[1]);
    if (n < 2 || n > 9) { fprintf(stderr, "n 2..9\n"); return 1; }
    nfree = 0;
    for (int i = 1; i <= n; i++)
        for (int j = i + 1; j <= n; j++)
            if (!(i == 1 && j == n)) { fi[nfree] = i; fj[nfree] = j; nfree++; }
    long long total = 1LL << nfree;
    unsigned long long basic = 0, hist[64] = {0};
    unsigned __int128 tc = 0;

    #pragma omp parallel
    {
        unsigned long long lb = 0, lh[64] = {0};
        unsigned __int128 ltc = 0;
        #pragma omp for schedule(dynamic, 4096)
        for (long long s = 0; s < total; s++) {
            int in[MAXN + 2][MAXN + 2];
            memset(in, 0, sizeof in);
            for (int c = 1; c <= n; c++) in[c][c] = 1;
            in[1][n] = 1;
            for (int b = 0; b < nfree; b++) if (s >> b & 1) in[fi[b]][fj[b]] = 1;
            int ok = 1;
            // nothing lost
            for (int i = 1; i <= n && ok; i++)
                for (int j = i + 1; j <= n && ok; j++) {
                    if (!in[i][j]) continue;
                    int k = j - 1; while (!in[i][k]) k--;      // top product [i,k]
                    int l = i + 1; while (!in[l][j]) l++;      // bottom product [l,j]
                    if (l > k + 1) ok = 0;
                }
            if (!ok) continue;
            // everything comes out of a column, and count m on the way
            int m = 0;
            for (int i = 1; i <= n && ok; i++)
                for (int j = i; j <= n && ok; j++) {
                    if (!in[i][j] || (i == 1 && j == n)) continue;
                    int up_row = 0, up_col = 0;
                    for (int jj = j + 1; jj <= n; jj++) if (in[i][jj]) { up_row = 1; break; }
                    for (int ii = i - 1; ii >= 1; ii--) if (in[ii][j]) { up_col = 1; break; }
                    if (!up_row && !up_col) ok = 0;
                    else if (i < j && (up_row ^ up_col)) m++;
                }
            if (!ok) continue;
            lb++; lh[m]++; ltc += ((unsigned __int128)1 << m) - 1;
        }
        #pragma omp critical
        { basic += lb; tc += ltc; for (int t = 0; t < 64; t++) hist[t] += lh[t]; }
    }
    unsigned long long tclo = (unsigned long long)tc;   // fits for n <= 9
    printf("n=%d basic=%llu tc=%llu hist:", n, basic, tclo);
    for (int t = 0; t < 64; t++) if (hist[t]) printf(" m%d=%llu", t, hist[t]);
    printf("\n");
    return 0;
}
