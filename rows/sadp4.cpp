// counts the basic configurations and the thermally coupled ones row by row, with the states between
// rows kept on disk. this is the one that got n = 16.
//
// row i is the set of streams that start at component i. a config gets filled in one row at a time and
// all that has to be remembered between rows is a 5-bit field per column that is still open, so configs
// with the same fields get summed into one state. the cost is the number of states, not configs.
// counts are kept mod 2^128 plus a copy mod 2^61 - 1, crt.py puts them back together.
//
// a row goes chunk by chunk: push a chunk of states through the row in memory, sort, write it out as a
// sorted run. then the runs get merged (equal states summed) into the next row's file.
// it picks up again after the last finished row (row_<i>.done) and stops if the disk gets too full.
//
//   g++ -O3 -std=c++17 -static -o sadp4.exe sadp4.cpp
//   (no -march=native, that build crashed for me from n = 10 on. -static because git bash was picking
//    up an old libstdc++ without std::filesystem)
//   sadp4.exe n chunk_states workdir [headroom_GB] [--break-merge]
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <queue>
#include <chrono>
#include <filesystem>

namespace fs = std::filesystem;
typedef unsigned __int128 u128;
static const uint64_t P = (1ULL << 61) - 1;
static inline uint64_t addp(uint64_t a, uint64_t b) { uint64_t s = a + b; return s >= P ? s - P : s; }

// D comes out to 64 bytes, no packing needed
struct E {          // state inside a row
    u128 cols;
    uint8_t r, pend;
    u128 cnt, w;
    uint64_t cp, wp;
};
struct D {          // state between rows, on disk
    u128 cols;
    u128 cnt, w;
    uint64_t cp, wp;
};


static int n;
static std::string dir;
static double headroom_gb = 130.0;
static bool BREAK_MERGE = false;   // test switch: the merge stops summing equal states, so the answer has to come out wrong

static inline int getc_(u128 c, int k) { return (int)((c >> (5 * k)) & 31); }
static inline u128 setc_(u128 c, int k, int v) { return (c & ~((u128)31 << (5 * k))) | ((u128)v << (5 * k)); }
static inline bool kl(const E& a, const E& b) {
    if (a.cols != b.cols) return a.cols < b.cols;
    if (a.r != b.r) return a.r < b.r;
    return a.pend < b.pend;
}
static inline bool ke(const E& a, const E& b) { return a.cols == b.cols && a.r == b.r && a.pend == b.pend; }
static void mergeE(std::vector<E>& v) {
    std::sort(v.begin(), v.end(), kl);
    size_t o = 0;
    for (size_t i = 0; i < v.size(); i++) {
        if (o && ke(v[o - 1], v[i])) {
            v[o - 1].cnt += v[i].cnt; v[o - 1].w += v[i].w;
            v[o - 1].cp = addp(v[o - 1].cp, v[i].cp); v[o - 1].wp = addp(v[o - 1].wp, v[i].wp);
        } else v[o++] = v[i];
    }
    v.resize(o);
}
static void p128(FILE* f, u128 x) {
    char b[64]; int p = 63; b[p] = 0;
    if (!x) b[--p] = '0';
    while (x) { b[--p] = '0' + (int)(x % 10); x /= 10; }
    fprintf(f, "%s", b + p);
}
static void die(const char* m) { fprintf(stderr, "FATAL: %s\n", m); exit(2); }
static void check_space(uint64_t need_bytes) {
    fs::space_info s = fs::space(dir);
    double free_gb = (double)s.available / 1e9, need_gb = (double)need_bytes / 1e9;
    if (free_gb - need_gb < headroom_gb) {
        fprintf(stderr, "FATAL: disk would drop to %.1f GB free (headroom %.0f)\n", free_gb - need_gb, headroom_gb);
        exit(4);
    }
}
static std::string rowfile(int i) { return dir + "/row_" + std::to_string(i) + ".bin"; }
static std::string runfile(int i, int k) { return dir + "/run_" + std::to_string(i) + "_" + std::to_string(k) + ".bin"; }

// push a chunk of states through row i. what comes out is sorted and summed but still has the full key,
// the column that closes gets dropped in the last merge of the row
static void expand_chunk(std::vector<E>& cur, int i, std::vector<E>& outv) {
    std::vector<E> nxt;
    for (E& e : cur) { e.r = (uint8_t)i; e.pend = 0; }
    for (int j = i + 1; j <= n; j++) {
        int k = j - i - 1;
        bool forced = (i == 1 && j == n);
        bool interior = (i >= 2 && j <= n - 1);
        nxt.clear();
        nxt.reserve(cur.size() * 2);
        for (const E& e : cur) {
            int cs = getc_(e.cols, k);
            if (!forced && cs != i + 2) nxt.push_back(e);
            int d = e.r + 1;
            int ncs = (d >= j) ? 1 : d + 2;
            E f = e;
            f.cols = setc_(e.cols, k, ncs);
            f.r = (uint8_t)j;
            f.pend = (interior && cs == 0) ? 1 : 0;
            if (!forced && cs == 0) { f.w <<= 1; f.wp = addp(f.wp, f.wp); }
            nxt.push_back(f);
        }
        mergeE(nxt);
        std::swap(cur, nxt);
    }
    outv.clear();
    for (E e : cur) {
        if (e.pend) continue;
        if (e.r > i && !(i == 1 && e.r == n)) { e.w <<= 1; e.wp = addp(e.wp, e.wp); }
        e.r = 0; e.pend = 0;
        outv.push_back(e);
    }
    mergeE(outv);
}

struct Reader {
    FILE* f = nullptr; std::vector<D> buf; size_t pos = 0, len = 0; bool eof = false;
    void open(const std::string& p) { f = fopen(p.c_str(), "rb"); if (!f) die("open run"); buf.resize(8192); fill(); }
    void fill() { len = fread(buf.data(), sizeof(D), buf.size(), f); pos = 0; if (len == 0) eof = true; }
    bool has() const { return !eof; }
    const D& cur() const { return buf[pos]; }
    void next() { if (++pos >= len) fill(); }
    void close() { if (f) fclose(f); f = nullptr; }
};
struct Writer {
    FILE* f = nullptr; std::vector<D> buf; size_t len = 0; uint64_t count = 0;
    void open(const std::string& p) { f = fopen(p.c_str(), "wb"); if (!f) die("open out"); buf.resize(8192); }
    void put(const D& d) { buf[len++] = d; count++; if (len == buf.size()) flush(); }
    void flush() { if (len && fwrite(buf.data(), sizeof(D), len, f) != len) die("write"); len = 0; }
    void close() { flush(); if (f) fclose(f); f = nullptr; }
};

static inline void addD(D& a, const D& b) {
    a.cnt += b.cnt; a.w += b.w; a.cp = addp(a.cp, b.cp); a.wp = addp(a.wp, b.wp);
}

// merge sorted files, summing equal keys. on the last merge of a row it also counts the distinct states,
// then drops the column that just closed (the lowest 5 bits) and sums again. the states that become
// equal are already next to each other since that field is the lowest one.
static uint64_t kmerge(const std::vector<std::string>& ins, const std::string& out, bool final_stage, int i) {
    std::vector<Reader> rd(ins.size());
    for (size_t k = 0; k < ins.size(); k++) rd[k].open(ins[k]);
    auto cmp = [&](int a, int b) { return rd[a].cur().cols > rd[b].cur().cols; };
    std::priority_queue<int, std::vector<int>, decltype(cmp)> pq(cmp);
    for (size_t k = 0; k < ins.size(); k++) if (rd[k].has()) pq.push((int)k);
    bool shift = final_stage && (i + 1 <= n - 1);
    Writer w; w.open(out);
    uint64_t distinct = 0;
    bool have = false, shave = false; D acc, sacc;
    memset(&acc, 0, sizeof acc); memset(&sacc, 0, sizeof sacc);
    auto emit = [&](const D& x) {
        distinct++;
        if (!shift) { w.put(x); return; }
        D s = x; s.cols >>= 5;
        if (shave && s.cols == sacc.cols) addD(sacc, s);
        else { if (shave) w.put(sacc); sacc = s; shave = true; }
    };
    while (!pq.empty()) {
        int k = pq.top(); pq.pop();
        const D& d = rd[k].cur();
        if (have && d.cols == acc.cols) { if (!BREAK_MERGE) addD(acc, d); }
        else { if (have) emit(acc); acc = d; have = true; }
        rd[k].next();
        if (rd[k].has()) pq.push(k);
    }
    if (have) emit(acc);
    if (shave) w.put(sacc);
    w.close();
    for (auto& r : rd) r.close();
    return distinct;
}

int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: sadp4 n chunk dir [headroom_GB] [--break-merge]\n"); return 1; }
    n = atoi(argv[1]);
    size_t chunk = strtoull(argv[2], 0, 10);
    dir = argv[3];
    for (int a = 4; a < argc; a++) {
        // anything that isn't the test switch or a plain number gets refused. a typo shouldn't turn into
        // headroom 0 (what you get from junk) and switch the disk check off
        if (!strcmp(argv[a], "--break-merge")) { BREAK_MERGE = true; continue; }
        char* end = nullptr;
        double v = strtod(argv[a], &end);
        if (end == argv[a] || *end != 0 || v < 0) { fprintf(stderr, "FATAL: unknown argument '%s'\n", argv[a]); return 1; }
        headroom_gb = v;
    }
    if (n < 3 || n > 24) die("3 <= n <= 24");
    fs::create_directories(dir);
    auto t0 = std::chrono::steady_clock::now();

    // start from one empty state, or from the last finished row
    int start = 1;
    for (int i = n - 1; i >= 1; i--) if (fs::exists(dir + "/row_" + std::to_string(i + 1) + ".done")) { start = i + 1; break; }
    if (start == 1) {
        Writer w; w.open(rowfile(1));
        D d; memset(&d, 0, sizeof d); d.cnt = 1; d.w = 1; d.cp = 1; d.wp = 1;
        w.put(d); w.close();
        FILE* f = fopen((dir + "/row_1.done").c_str(), "w"); fclose(f);
    } else fprintf(stderr, "resuming at row %d\n", start);

    for (int i = start; i < n; i++) {
        // leftovers from a killed run of this row never get read, they just take up disk
        {
            std::string pre[3] = {"run_" + std::to_string(i) + "_", "cmp_" + std::to_string(i) + "_",
                                  "mrg_" + std::to_string(i) + "_"};
            std::vector<fs::path> stale;
            for (auto& de : fs::directory_iterator(dir)) {
                std::string nm = de.path().filename().string();
                for (auto& p : pre) if (nm.rfind(p, 0) == 0) stale.push_back(de.path());
            }
            for (auto& p : stale) fs::remove(p);
            fs::remove(rowfile(i + 1) + ".part");
        }
        // pass 1: chunks -> sorted runs
        FILE* in = fopen(rowfile(i).c_str(), "rb"); if (!in) die("open row");
        std::vector<D> dbuf(chunk);
        std::vector<E> cur, outv;
        int nruns = 0, ncompact = 0;
        // every 16 runs at one level get merged into one run at the next level, so copies of the same
        // state from different chunks collapse early and the disk only holds about the distinct states
        const size_t BATCH = 16;
        std::vector<std::vector<std::string>> lv(1);
        size_t got;
        while ((got = fread(dbuf.data(), sizeof(D), chunk, in)) > 0) {
            cur.resize(got);
            for (size_t q = 0; q < got; q++) {
                cur[q].cols = dbuf[q].cols; cur[q].r = 0; cur[q].pend = 0;
                cur[q].cnt = dbuf[q].cnt; cur[q].w = dbuf[q].w; cur[q].cp = dbuf[q].cp; cur[q].wp = dbuf[q].wp;
            }
            expand_chunk(cur, i, outv);
            check_space((uint64_t)outv.size() * sizeof(D));
            std::string rf = runfile(i, nruns++);
            Writer w; w.open(rf);
            for (const E& e : outv) { D d; d.cols = e.cols; d.cnt = e.cnt; d.w = e.w; d.cp = e.cp; d.wp = e.wp; w.put(d); }
            w.close();
            std::vector<E>().swap(cur); std::vector<E>().swap(outv);
            lv[0].push_back(rf);
            for (size_t L = 0; L < lv.size() && lv[L].size() >= BATCH; L++) {
                uint64_t bytes = 0;
                for (auto& s : lv[L]) bytes += fs::file_size(s);
                check_space(bytes);
                std::string o = dir + "/cmp_" + std::to_string(i) + "_" + std::to_string(ncompact++) + ".bin";
                kmerge(lv[L], o, false, i);
                for (auto& s : lv[L]) fs::remove(s);
                lv[L].clear();
                if (L + 1 == lv.size()) lv.emplace_back();
                lv[L + 1].push_back(o);
            }
        }
        fclose(in);
        std::vector<D>().swap(dbuf);

        // pass 2: merge what's left, 200 files at a time at most (the C runtime caps open files)
        const size_t MAXOPEN = 200;
        std::vector<std::string> ins;
        for (auto& l : lv) for (auto& s : l) ins.push_back(s);
        int level = 0;
        while (ins.size() > MAXOPEN) {
            std::vector<std::string> next;
            uint64_t bytes = 0;
            for (auto& s : ins) bytes += fs::file_size(s);
            check_space(bytes);
            for (size_t g = 0; g < ins.size(); g += MAXOPEN) {
                std::vector<std::string> grp(ins.begin() + g, ins.begin() + std::min(ins.size(), g + MAXOPEN));
                std::string o = dir + "/mrg_" + std::to_string(i) + "_" + std::to_string(level) + "_" + std::to_string(g / MAXOPEN) + ".bin";
                kmerge(grp, o, false, i);
                for (auto& s : grp) fs::remove(s);
                next.push_back(o);
            }
            ins.swap(next);
            level++;
        }
        uint64_t nout = kmerge(ins, rowfile(i + 1) + ".part", true, i);
        for (auto& s : ins) fs::remove(s);
        fs::rename(rowfile(i + 1) + ".part", rowfile(i + 1));
        FILE* f = fopen((dir + "/row_" + std::to_string(i + 1) + ".done").c_str(), "w"); fclose(f);
        if (i >= 2) fs::remove(rowfile(i));                 // only row 1 and the newest row are kept
        double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        // "between" is the number of distinct states after row i, counted before the closed column is
        // dropped. the file ends up with fewer records than that, so don't check it against the file size
        fprintf(stderr, "row %d between %llu   (runs %d, compactions %d, final merge levels %d, %.0fs)\n", i,
                (unsigned long long)nout, nruns, ncompact, level + 1, s);
        fflush(stderr);
    }

    // the last row's file is the finished configs
    u128 basic = 0, wsum = 0; uint64_t bp = 0, wpp = 0;
    {
        Reader r; r.open(rowfile(n));
        while (r.has()) { const D& d = r.cur(); basic += d.cnt; wsum += d.w; bp = addp(bp, d.cp); wpp = addp(wpp, d.wp); r.next(); }
        r.close();
    }
    double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    printf("RES n=%d basic_m128=", n); p128(stdout, basic); printf(" basic_mP=%llu sum_m128=", (unsigned long long)bp);
    p128(stdout, wsum); printf(" sum_mP=%llu\n", (unsigned long long)wpp);
    printf("%d ", n); p128(stdout, basic); printf(" "); p128(stdout, wsum - basic); printf("   (%.1fs, mod 2^128 values; rebuild with crt.py)\n", s);
    return 0;
}
