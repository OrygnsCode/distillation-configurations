// counts the basic configurations and the thermally coupled ones row by row. it's a second program for
// the same sweep as sadp4.cpp, written separately, and it got n = 16 too: the two agree on both counts
// and on the number of states after every row.
//
// row i is the streams that start at component i. between rows all that's left of a config is a 5-bit
// field per column that is still open (a state), and configs with the same state get counted together.
//
// what's different from sadp4:
//   - the counts are exact, no CRT. the basic count is carried in 128 bits and the weighted one (basic +
//     thermally coupled) in 192, as three 64-bit words.
//   - a row's states are summed in a hash table instead of sorted runs. a row is done in P passes over
//     the previous row's files. pass p keeps only the states whose key hashes to p mod P, sums them in
//     the table, sorts them and writes one file. so memory is one table whatever n is, and the price is
//     scanning the row P times. P doubles, and the row starts over, when the table gets over 80% full.
//
// it can be stopped any time. done.txt lists the files of the last finished row. inside a row,
// rowR_progress.txt lists the passes that are finished, each file with its size and a checksum, and a
// restart keeps those passes if the files still check out. if they don't, it redoes the row.
//
//   g++ -O2 -std=c++17 -static -o rowdp_disk.exe rowdp_disk.cpp -lpsapi      (Windows, MinGW)
//   rowdp_disk.exe n M log2_slots workdir headroom_GB maxmem_GB
//     M           states are scanned in groups that agree on their last M columns (4 for n = 16)
//     log2_slots  the table has 2^log2_slots slots of 56 bytes, so 26 is 3.76 GB
//     it stops before a write that would leave less than headroom_GB free on the disk, and if its
//     memory goes over maxmem_GB.
//   three more arguments are only for testing: a shift (every weight starts at 2^shift, to get past 128
//   bits at small n), and a row and a pass to stop after, as if it had been killed there.
//   n = 16 was  rowdp_disk.exe 16 4 26 workdir 130 9  : about 19 hours, 6 GB of memory, 20 GB of disk.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <chrono>
typedef unsigned __int128 K;
typedef unsigned __int128 U;
typedef uint64_t u64;
struct W { u64 a[3]; };   // 192-bit integer, low word first
static bool w_overflow = false;
static inline void wadd(W& x, const W& y) {
    unsigned __int128 t = (unsigned __int128)x.a[0] + y.a[0]; x.a[0] = (u64)t;
    t = (unsigned __int128)x.a[1] + y.a[1] + (u64)(t >> 64); x.a[1] = (u64)t;
    t = (unsigned __int128)x.a[2] + y.a[2] + (u64)(t >> 64); x.a[2] = (u64)t;
    if (t >> 64) w_overflow = true;
}
static inline void wdbl(W& x) {
    if (x.a[2] >> 63) w_overflow = true;
    x.a[2] = (x.a[2] << 1) | (x.a[1] >> 63); x.a[1] = (x.a[1] << 1) | (x.a[0] >> 63); x.a[0] <<= 1;
}
static void wsub128(W& x, U y) {   // x -= y, for y below 2^128 and not above x
    u64 ylo = (u64)y, yhi = (u64)(y >> 64);
    unsigned __int128 t = (unsigned __int128)x.a[0] - ylo; x.a[0] = (u64)t; u64 b = (t >> 64) ? 1 : 0;
    t = (unsigned __int128)x.a[1] - yhi - b; x.a[1] = (u64)t; b = (t >> 64) ? 1 : 0;
    x.a[2] -= b;
}
static void wprint(W x) {
    char buf[80]; int p = 79; buf[p] = 0;
    bool zero = !(x.a[0] | x.a[1] | x.a[2]);
    if (zero) buf[--p] = '0';
    while (x.a[0] | x.a[1] | x.a[2]) {
        unsigned __int128 rem = 0;
        for (int i = 2; i >= 0; i--) { unsigned __int128 cur = (rem << 64) | x.a[i]; x.a[i] = (u64)(cur / 10); rem = cur % 10; }
        buf[--p] = '0' + (int)rem;
    }
    printf("%s", buf + p);
}
struct E { u64 klo, khi, clo, chi; W w; };   // a state: key, basic count, weighted count. 56 bytes
static inline K key(const E& e) { return ((K)e.khi << 64) | e.klo; }
static inline U cnt(const E& e) { return ((U)e.chi << 64) | e.clo; }
static inline E mk(K k, U c, const W& w) { E e; e.klo = (u64)k; e.khi = (u64)(k >> 64); e.clo = (u64)c; e.chi = (u64)(c >> 64); e.w = w; return e; }
static bool c_overflow = false;
// sort by key and add up the states that are equal
static void merge(std::vector<E>& v) {
    std::sort(v.begin(), v.end(), [](const E& a, const E& b) { return a.khi != b.khi ? a.khi < b.khi : a.klo < b.klo; });
    size_t o = 0;
    for (size_t r = 0; r < v.size();) {
        K k = key(v[r]); U c = 0; W w = {{0, 0, 0}};
        while (r < v.size() && key(v[r]) == k) {
            U x = cnt(v[r]); if (c + x < c) c_overflow = true; c += x;
            wadd(w, v[r].w); r++;
        }
        v[o++] = mk(k, c, w);
    }
    v.resize(o);
}
// the key: 5 bits per column, then at bit 100 the column of the last stream put in this row, and one flag bit
static inline int col(K k, int j) { return (int)((k >> (5 * (j - 1))) & 31); }
static inline K setcol(K k, int j, int v) { K m = (K)31 << (5 * (j - 1)); return (k & ~m) | ((K)v << (5 * (j - 1))); }
static inline int getk(K k) { return (int)((k >> 100) & 31); }
static inline int gett(K k) { return (int)((k >> 105) & 1); }
static inline K setkt(K k, int kk, int t) { k &= ~((K)63 << 100); return k | ((K)kk << 100) | ((K)t << 105); }
static void print128(U x) { char b[64]; int p = 63; b[p] = 0; if (!x) b[--p] = '0'; while (x) { b[--p] = '0' + (int)(x % 10); x /= 10; } printf("%s", b + p); }

static int n;
// push a batch of states through row r, one column at a time. returns the most states it held at once
static size_t scan_row(int r, const E* b, const E* e, std::vector<E>& out, std::vector<E>& cur, std::vector<E>& nxt) {
    cur.clear();
    for (const E* p = b; p != e; ++p) {
        K k0 = key(*p);
        int top = col(k0, r) == 0;   // the cell (r,r) on the diagonal
        cur.push_back(mk(setkt(setcol(k0, r, 0), r, top), cnt(*p), p->w));
    }
    size_t peak = cur.size();
    for (int j = r + 1; j <= n; j++) {
        bool must = (r == 1 && j == n);
        nxt.clear();
        for (auto& s : cur) {
            K sk = key(s);
            int v = col(sk, j), k = getk(sk), t = gett(sk);
            bool due = (v >= 2 && v - 1 == r);   // column j has to get a stream in this row at the latest
            if (!must && !due) nxt.push_back(s);
            int nv = (k + 1 >= j) ? 1 : 1 + (k + 1);
            W w = s.w;
            if (k > r && t) wdbl(w);   // the stream before this one comes out of one split only: it can be coupled or not
            nxt.push_back(mk(setkt(setcol(sk, j, nv), j, v == 0), cnt(s), w));
        }
        merge(nxt); cur.swap(nxt);
        peak = std::max(peak, cur.size());
    }
    for (auto& s : cur) {
        K sk = key(s);
        int k = getk(sk), t = gett(sk);
        bool feed = (r == 1 && k == n);
        if (t && !feed) continue;   // the last stream of the row would have nothing that makes it
        W w = s.w;
        if (k > r && !feed) wdbl(w);   // the last stream of the row comes out of one split only too
        out.push_back(mk(setkt(sk, 0, 0), cnt(s), w));
    }
    return peak;
}

#include <string>
#include <cstring>
#include <filesystem>
#include <thread>
#include <io.h>
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <psapi.h>
namespace fs = std::filesystem;

// stops the program if its memory goes over the limit. looks every 50 ms
static void watchdog(double limit_gb) {
    for (;;) {
        PROCESS_MEMORY_COUNTERS pmc;
        if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof pmc)) {
            double gb = pmc.PagefileUsage / 1e9;
            if (gb > limit_gb) {
                printf("STOPPED: memory %.2f GB passed the limit of %.2f GB\n", gb, limit_gb);
                fflush(stdout);
                _exit(7);
            }
        }
        Sleep(50);
    }
}

static const u64 EMPTY = ~0ULL;   // marks an empty slot. real keys are below 2^107

static inline u64 hmix(u64 lo, u64 hi) {
    u64 x = lo ^ (hi * 0x9e3779b97f4a7c15ULL);
    x ^= x >> 31; x *= 0xbf58476d1ce4e5b9ULL; x ^= x >> 29; x *= 0x94d049bb133111ebULL; return x ^ (x >> 32);
}
// which pass a key belongs to. a second hash, so the keys of one pass still spread over the whole table
static inline u64 hpart(u64 lo, u64 hi) { return hmix(hi ^ 0x632be59bd9b4e019ULL, lo + 0x2545f4914f6cdd1dULL); }

struct Table {
    size_t cap, mask, used = 0;
    std::vector<E> slot;
    explicit Table(size_t c) : cap(c), mask(c - 1), slot(c) { clear(); }
    void clear() { for (auto& s : slot) s.khi = EMPTY; used = 0; }
    bool add(const E& e) {   // false means the table is too full
        size_t i = hmix(e.klo, e.khi) & mask;
        for (;;) {
            E& s = slot[i];
            if (s.khi == EMPTY) {
                if ((used + 1) * 10 > cap * 8) return false;
                s = e; used++; return true;
            }
            if (s.khi == e.khi && s.klo == e.klo) {
                U a = cnt(s), b = cnt(e), c = a + b;
                if (c < a) c_overflow = true;
                s.clo = (u64)c; s.chi = (u64)(c >> 64);
                wadd(s.w, e.w);
                return true;
            }
            i = (i + 1) & mask;
        }
    }
    // used slots to the front, sorted by key. the table can't be used after this
    size_t compact_sorted() {
        size_t o = 0;
        for (size_t i = 0; i < cap; i++) if (slot[i].khi != EMPTY) slot[o++] = slot[i];
        std::sort(slot.begin(), slot.begin() + o, [](const E& a, const E& b) { return a.khi != b.khi ? a.khi < b.khi : a.klo < b.klo; });
        return o;
    }
};

struct Reader {
    FILE* f = nullptr; std::vector<E> buf; size_t pos = 0, len = 0; bool eof = false;
    bool open(const std::string& p) { f = fopen(p.c_str(), "rb"); buf.resize(1 << 18); pos = len = 0; eof = false; return f != nullptr; }
    bool fill() { if (eof) return false; len = fread(buf.data(), sizeof(E), buf.size(), f); pos = 0; if (len == 0) { eof = true; return false; } return true; }
    bool peek(E& e) { if (pos == len && !fill()) return false; e = buf[pos]; return true; }
    bool group(std::vector<E>& g, K gmask) {   // the next run of states that agree on their last M columns
        g.clear(); E e;
        if (!peek(e)) return false;
        K gk = key(e) & gmask;
        while (peek(e) && (key(e) & gmask) == gk) { g.push_back(e); pos++; }
        return true;
    }
    void close() { if (f) fclose(f); f = nullptr; }
};

static void write_file(const std::string& path, const E* p, size_t m) {
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) { printf("CANNOT WRITE %s\n", path.c_str()); exit(4); }
    if (m && fwrite(p, sizeof(E), m, f) != m) { printf("SHORT WRITE %s\n", path.c_str()); exit(4); }
    // on the disk before anything lists the file as finished
    if (fflush(f) != 0 || _commit(_fileno(f)) != 0) { printf("CANNOT FLUSH %s\n", path.c_str()); exit(4); }
    fclose(f);
}

static std::string part_name(int r, int p) { return "row" + std::to_string(r) + "_part" + std::to_string(p) + ".bin"; }

// checksum of a pass file, one state at a time
static inline u64 rec_sum(u64 h, const E& e) {
    h = hmix(h ^ e.klo, e.khi); h = hmix(h ^ e.clo, e.chi); h = hmix(h ^ e.w.a[0], e.w.a[1]);
    return hmix(h, e.w.a[2]);
}
// the checksum starts from n, the row, P, the pass and the shift, so it only fits the file it was made for.
// without that a progress file with P cut in half would get through: the files of 16 passes are halves of the files of 8
static inline u64 part_seed(int r, int P, int p, int shift) {
    return hmix((u64)n * 65536 + (u64)r, ((u64)P * 65536 + (u64)p) * 256 + (u64)shift);
}

// the finished passes of row r: P, pass 0's counters, then size, checksum and name of each file
static void write_progress(const std::string& dir, int r, int P, size_t groups, size_t gpeak, size_t appended,
                           const std::vector<size_t>& sizes, const std::vector<u64>& sums) {
    std::string pf = dir + "/row" + std::to_string(r) + "_progress.txt";
    FILE* d = fopen((pf + ".tmp").c_str(), "w");
    if (!d) { printf("CANNOT WRITE %s.tmp\n", pf.c_str()); exit(4); }
    fprintf(d, "n %d row %d P %d groups %llu gpeak %llu appended %llu\n", n, r, P,
            (unsigned long long)groups, (unsigned long long)gpeak, (unsigned long long)appended);
    for (size_t i = 0; i < sizes.size(); i++)
        fprintf(d, "%llu %016llx %s\n", (unsigned long long)sizes[i], (unsigned long long)sums[i], part_name(r, (int)i).c_str());
    fclose(d);
    fs::rename(pf + ".tmp", pf);
}

// a pass file left by a stopped run is only used again if it has the right size and checksum, its keys go
// up, and every key belongs to that pass
static bool check_part(const std::string& path, size_t m, u64 want, u64 seed, int p, int P) {
    std::error_code ec;
    auto sz = fs::file_size(path, ec);
    if (ec || sz != (uintmax_t)m * sizeof(E)) return false;
    Reader rd; E e, prev = {}; bool first = true; size_t seen = 0; u64 sum = seed;
    if (!rd.open(path)) return false;
    while (rd.peek(e)) {
        bool good = (P == 1 || (int)(hpart(e.klo, e.khi) % (u64)P) == p) && (first || key(prev) < key(e));
        if (!good) { rd.close(); return false; }
        sum = rec_sum(sum, e);
        prev = e; first = false; seen++; rd.pos++;
    }
    rd.close();
    return seen == m && sum == want;
}

int main(int argc, char** argv) {
    // a number that doesn't parse all the way is an error, never 0
    auto num = [&](int i, double lo, double hi) {
        char* end = nullptr;
        double v = strtod(argv[i], &end);
        if (end == argv[i] || *end || v < lo || v > hi) { fprintf(stderr, "argument %d (%s) must be a number in [%g, %g]\n", i, argv[i], lo, hi); exit(1); }
        return v;
    };
    if (argc < 7 || argc > 10) {
        fprintf(stderr, "rowdp_disk.exe N M CAPLOG2 DIR HEADROOM_GB MAXMEM_GB [SHIFT] [STOPAFTER] [STOPPART]\n");
        return 1;
    }
    n = (int)num(1, 2, 20);
    int M = (int)num(2, 1, n), caplog = (int)num(3, 8, 31);
    std::string dir = argv[4];
    double headroom = num(5, 0, 10000);
    double maxmem = num(6, 0.01, 1000);
    int shift = argc > 7 ? (int)num(7, 0, 150) : 0;
    int stopafter = argc > 8 ? (int)num(8, 0, 20) : 0;   // testing: stop after this row
    int stoppart = argc > 9 ? (int)num(9, 0, 65536) : -1;   // testing: stop inside that row, after this pass
    std::thread(watchdog, maxmem).detach();
    fs::create_directories(dir);
    auto t0 = std::chrono::steady_clock::now();
    K gmask = ~(((K)1 << (5 * (n - M))) - 1) & (((K)1 << 100) - 1);   // the last M columns
    Table T((size_t)1 << caplog);

    std::vector<std::string> inputs;
    int start = 1, P = 1;
    std::string donef = dir + "/done.txt";
    if (FILE* d = fopen(donef.c_str(), "r")) {   // pick up after the last finished row
        int rd, pd, sd = 0; char line[4096];   // no shift field means shift 0
        if (fgets(line, sizeof line, d) && sscanf(line, "row %d P %d shift %d", &rd, &pd, &sd) >= 2) {
            if (sd != shift) {   // the weights in those files are 2^sd times too big for this run
                fprintf(stderr, "%s was made with shift %d, this run has shift %d: not resumed\n", donef.c_str(), sd, shift);
                return 1;
            }
            while (fgets(line, sizeof line, d)) { line[strcspn(line, "\r\n")] = 0; if (*line) inputs.push_back(line); }
            start = rd + 1; P = pd;
            printf("resuming after row %d with %zu files\n", rd, inputs.size());
        }
        fclose(d);
    }
    if (start == 1) {
        W one = {{1, 0, 0}};
        for (int i = 0; i < shift; i++) wdbl(one);
        E e = mk(0, 1, one);
        std::string p0 = dir + "/row0_part0.bin";
        write_file(p0, &e, 1);
        inputs = {p0};
    }
    std::vector<E> g, out, cur, nxt;
    for (int r = start; r <= n; r++) {
        auto tr = std::chrono::steady_clock::now();
        size_t states_in = 0;
        for (auto& f : inputs) states_in += (size_t)fs::file_size(f) / sizeof(E);
        std::vector<std::string> outs;
        std::vector<size_t> sizes;
        std::vector<u64> sums;
        size_t groups = 0, gpeak = 0, appended = 0, states_out = 0;
        int first = 0;   // passes a stopped run already finished
        fs::remove(dir + "/row" + std::to_string(r - 1) + "_progress.txt");
        if (FILE* d = fopen((dir + "/row" + std::to_string(r) + "_progress.txt").c_str(), "r")) {
            int nn, rr, pp; unsigned long long g0, k0, a0; char line[4096];
            bool good = fscanf(d, "n %d row %d P %d groups %llu gpeak %llu appended %llu\n", &nn, &rr, &pp, &g0, &k0, &a0) == 6
                        && nn == n && rr == r && pp >= P && pp <= (1 << 16);
            std::vector<size_t> ks; std::vector<u64> kh;
            while (good && fgets(line, sizeof line, d)) {
                line[strcspn(line, "\r\n")] = 0;
                if (!*line) continue;
                char *end = nullptr, *end2 = nullptr;
                unsigned long long m = strtoull(line, &end, 10);
                unsigned long long h = (end != line && *end == ' ') ? strtoull(end + 1, &end2, 16) : 0;
                int p = (int)ks.size();
                good = end2 && end2 != end + 1 && *end2 == ' ' && part_name(r, p) == end2 + 1 && p < pp
                       && check_part(dir + "/" + part_name(r, p), (size_t)m, (u64)h, part_seed(r, pp, p, shift), p, pp);
                if (good) { ks.push_back((size_t)m); kh.push_back((u64)h); }
            }
            fclose(d);
            if (good) {
                P = pp; sizes = ks; sums = kh; first = (int)sizes.size();
                for (int p = 0; p < first; p++) { outs.push_back(dir + "/" + part_name(r, p)); states_out += sizes[p]; }
                if (first > 0) { groups = (size_t)g0; gpeak = (size_t)k0; appended = (size_t)a0; }
                printf("row %2d: resuming at P = %d, %d parts already written and checked\n", r, P, first);
            } else
                printf("row %2d: the progress file does not match its files, the row starts again\n", r);
            fflush(stdout);
        }
        for (;;) {
            bool ok = true;
            // anything else lying around for this row goes: half-written files, files made with another P
            std::string pre = "row" + std::to_string(r) + "_part";
            std::vector<fs::path> stale;
            for (auto& de : fs::directory_iterator(dir)) {
                std::string nm = de.path().filename().string();
                if (nm.rfind(pre, 0) != 0) continue;
                bool keep = false;
                for (int p = 0; p < first && !keep; p++) keep = nm == part_name(r, p);
                if (!keep) stale.push_back(de.path());
            }
            for (auto& f : stale) fs::remove(f);
            write_progress(dir, r, P, groups, gpeak, appended, sizes, sums);
            for (int p = first; p < P && ok; p++) {
                T.clear();
                for (auto& f : inputs) {
                    Reader rd;
                    if (!rd.open(f)) { printf("CANNOT READ %s\n", f.c_str()); return 4; }
                    while (ok && rd.group(g, gmask)) {
                        out.clear();
                        size_t pk = scan_row(r, g.data(), g.data() + g.size(), out, cur, nxt);
                        if (p == 0) { groups++; gpeak = std::max(gpeak, pk); appended += out.size(); }
                        for (auto& o : out) {
                            if (P > 1 && (int)(hpart(o.klo, o.khi) % (u64)P) != p) continue;
                            if (!T.add(o)) { ok = false; break; }
                        }
                    }
                    rd.close();
                    if (!ok) break;
                }
                if (!ok) break;
                size_t m = T.compact_sorted();
                double need = (double)m * sizeof(E) / 1e9;
                double avail = (double)fs::space(dir).available / 1e9;
                if (avail - need < headroom) {
                    printf("STOPPED: writing row %d part %d (%.1f GB) would leave %.1f GB free, under %.0f GB. Resume later from done.txt.\n",
                           r, p, need, avail - need, headroom);
                    return 5;
                }
                std::string path = dir + "/" + part_name(r, p);
                write_file(path, T.slot.data(), m);
                outs.push_back(path); states_out += m;
                u64 sum = part_seed(r, P, p, shift);
                for (size_t i = 0; i < m; i++) sum = rec_sum(sum, T.slot[i]);
                sizes.push_back(m); sums.push_back(sum);
                write_progress(dir, r, P, groups, gpeak, appended, sizes, sums);
                if (P > 1) {
                    printf("row %2d part %d of %d written: %zu states, %.0fs into this run of the row\n", r, p + 1, P, m,
                           std::chrono::duration<double>(std::chrono::steady_clock::now() - tr).count());
                    fflush(stdout);
                }
                if (r == stopafter && p == stoppart) { printf("stopped in row %d after pass %d of %d (test)\n", r, p + 1, P); return 6; }
            }
            if (ok) break;
            for (auto& f : outs) fs::remove(f);
            outs.clear(); sizes.clear(); sums.clear(); states_out = 0; first = 0;
            if (P >= (1 << 16)) { printf("STOPPED: row %d needs more than %d partitions\n", r, P); return 2; }
            P *= 2;
            groups = gpeak = appended = 0;
            printf("row %2d: table over 80%% full, restarting the row with %d partitions\n", r, P);
            fflush(stdout);
        }
        int Pused = P;
        if (P > 1 && states_out * 10 < T.cap * 2) P = std::max(1, P / 2);   // the next row may get by with fewer
        FILE* d = fopen((donef + ".tmp").c_str(), "w");
        if (!d) { printf("CANNOT WRITE %s.tmp\n", donef.c_str()); return 4; }
        fprintf(d, "row %d P %d shift %d\n", r, P, shift);   // the P the next row starts with
        for (auto& f : outs) fprintf(d, "%s\n", f.c_str());
        fclose(d);
        fs::rename(donef + ".tmp", donef);
        fs::remove(dir + "/row" + std::to_string(r) + "_progress.txt");
        for (auto& f : inputs) fs::remove(f);
        inputs = outs;
        double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - tr).count();
        printf("row %2d: %zu states in, %zu groups, largest group scan %zu, %zu appended, %zu states out, %d partitions, %.1fs\n",
               r, states_in, groups, gpeak, appended, states_out, Pused, s);
        fflush(stdout);
        if (r == stopafter) { printf("stopped after row %d (test)\n", r); return 6; }
    }
    U a = 0; W w = {{0, 0, 0}};
    for (auto& f : inputs) {
        Reader rd; E e;
        if (!rd.open(f)) { printf("CANNOT READ %s\n", f.c_str()); return 4; }
        while (rd.peek(e)) { U x = cnt(e); if (a + x < a) c_overflow = true; a += x; wadd(w, e.w); rd.pos++; }
        rd.close();
    }
    double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    printf("n=%d basic=", n); print128(a);
    if (shift == 0) { W tc = w; wsub128(tc, a); printf(" TC="); wprint(tc); }
    printf(" W%s=", shift ? "_shifted" : ""); wprint(w);
    printf("   (table %zu slots = %.2f GB, %.1fs)%s%s\n", T.cap, T.cap * sizeof(E) / 1e9, sec,
           c_overflow ? " BASIC COUNT OVERFLOWED 128 BITS" : "", w_overflow ? " WEIGHTED COUNT OVERFLOWED 192 BITS" : "");
    return (c_overflow || w_overflow) ? 3 : 0;
}
