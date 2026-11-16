// the same count done a different way, as a check on the row sweep: streams get added shortest first.
// for each component it remembers the longest stream so far that starts there (a) and the longest that
// ends there (b), 4 bits each, plus which rows and columns still need a longer stream for their longest
// one to come out of.
// [i,j] can only go in if its two products (the longest shorter streams in row i and in column j) don't
// lose anything between them. after length d nothing new can start at n - d or end at d + 1, so that row
// and column get checked and dropped. a layer is done a chunk of states at a time to keep memory down.
// counts the basic configs and the thermally coupled ones, like the row version.
//
//   g++ -O3 -o lendp2.exe lendp2.cpp
//   ./lendp2.exe n [chunk_states]      (n = 12 needs about 8 GB)
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <chrono>

typedef unsigned __int128 u128;

#pragma pack(push, 1)
struct E {
    uint64_t ka, kb;   // a and b, 4 bits per component
    uint32_t kf;       // rows still waiting | columns still waiting << 16
    u128 cnt, w;
};
#pragma pack(pop)

static int n;
static inline int geta(const E& e, int i) { return (e.ka >> (4 * i)) & 15; }
static inline int getb(const E& e, int j) { return (e.kb >> (4 * j)) & 15; }
static inline void seta(E& e, int i, int v) { e.ka = (e.ka & ~(15ULL << (4 * i))) | ((uint64_t)v << (4 * i)); }
static inline void setb(E& e, int j, int v) { e.kb = (e.kb & ~(15ULL << (4 * j))) | ((uint64_t)v << (4 * j)); }
static inline bool kl(const E& x, const E& y) {
    if (x.ka != y.ka) return x.ka < y.ka;
    if (x.kb != y.kb) return x.kb < y.kb;
    return x.kf < y.kf;
}
static inline bool ke(const E& x, const E& y) { return x.ka == y.ka && x.kb == y.kb && x.kf == y.kf; }
static void merge(std::vector<E>& v) {
    std::sort(v.begin(), v.end(), kl);
    size_t o = 0;
    for (size_t i = 0; i < v.size(); i++) {
        if (o && ke(v[o - 1], v[i])) { v[o - 1].cnt += v[i].cnt; v[o - 1].w += v[i].w; }
        else v[o++] = v[i];
    }
    v.resize(o);
    v.shrink_to_fit();
}
static void p128(u128 x) {
    char b[64]; int p = 63; b[p] = 0;
    if (!x) b[--p] = '0';
    while (x) { b[--p] = '0' + (int)(x % 10); x /= 10; }
    printf("%s", b + p);
}

// row r and column cc are done after length d: check them and drop them. false means the config can't work
static bool finalize(E& e, int r, int cc) {
    uint32_t nr = e.kf & 0xFFFF, nc = e.kf >> 16;
    int ar = geta(e, r), ee = r + ar;
    if (nr >> r & 1) return false;
    if (ee <= cc) { if (ee == cc && !(getb(e, ee) > ar)) return false; }
    else if (getb(e, ee) == ar && ee != n) nc |= 1u << ee;
    if (nc >> cc & 1) return false;
    int bc = getb(e, cc), i0 = cc - bc;
    if (i0 < r && geta(e, i0) == bc && i0 != 1) nr |= 1u << i0;
    int mult = 0;
    if (r >= 2 && r <= n - 1 && ar >= 1) mult++;
    if (cc >= 2 && cc <= n - 1 && bc >= 1) mult++;
    e.w <<= mult;
    seta(e, r, 0); setb(e, cc, 0);
    nr &= ~(1u << r); nc &= ~(1u << cc);
    e.kf = nr | (nc << 16);
    return true;
}

int main(int argc, char** argv) {
    n = atoi(argv[1]);
    size_t chunk = argc > 2 ? strtoull(argv[2], 0, 10) : 1000000;
    if (n < 3 || n > 15) { fprintf(stderr, "3 <= n <= 15\n"); return 1; }
    auto t0 = std::chrono::steady_clock::now();
    std::vector<E> st(1), out, cur, nxt;
    st[0].ka = st[0].kb = 0; st[0].kf = 0; st[0].cnt = 1; st[0].w = 1;
    size_t peak_between = 1, peak_live = 1;
    for (int d = 1; d < n; d++) {
        bool forced = (d == n - 1);
        out.clear();
        for (size_t c0 = 0; c0 < st.size(); c0 += chunk) {
            size_t c1 = std::min(st.size(), c0 + chunk);
            cur.assign(st.begin() + c0, st.begin() + c1);
            for (int i = 1; i <= n - d; i++) {
                int j = i + d;
                nxt.clear();
                for (const E& e : cur) {
                    if (!forced) nxt.push_back(e);
                    if (d <= geta(e, i) + getb(e, j) + 1) {
                        E f = e;
                        seta(f, i, d); setb(f, j, d);
                        uint32_t nr = f.kf & 0xFFFF, nc = f.kf >> 16;
                        nr &= ~(1u << i); nc &= ~(1u << j);
                        f.kf = nr | (nc << 16);
                        nxt.push_back(f);
                    }
                }
                merge(nxt);
                std::swap(cur, nxt);
                size_t live = st.size() + out.size() + cur.size();
                if (live > peak_live) peak_live = live;
            }
            if (d < n - 1) {
                int r = n - d, cc = d + 1;
                for (E e : cur) if (finalize(e, r, cc)) out.push_back(e);
            } else {
                for (const E& e : cur) out.push_back(e);
            }
            if (out.size() > 4 * chunk) merge(out);
            std::vector<E>().swap(nxt);
        }
        std::vector<E>().swap(cur);
        merge(out);
        std::swap(st, out);
        std::vector<E>().swap(out);
        if (st.size() > peak_between) peak_between = st.size();
        fprintf(stderr, "layer %d between %zu\n", d, st.size());
    }
    u128 basic = 0, w = 0;
    for (const E& e : st) {
        uint32_t nr = e.kf & 0xFFFF, nc = e.kf >> 16;
        if ((nr & ~(1u << 1)) == 0 && (nc & ~(1u << n)) == 0) { basic += e.cnt; w += e.w; }
    }
    double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    printf("%d ", n); p128(basic); printf(" "); p128(w - basic);
    printf("   (between-layer peak %zu, live peak %zu, %zu-byte records, %.1fs)\n", peak_between, peak_live, sizeof(E), s);
    return 0;
}
