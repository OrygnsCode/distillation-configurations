# the three checks exactly as the matrix method prints them, done by brute force.
#
#   python printed_rules.py
#
# this is the tie between the counts here and the published method. it builds every candidate
# matrix and keeps the ones that pass the three checks as written, with no dynamic programming
# and nothing clever, so it can be read line by line against the quoted text. it stops at n = 7,
# which is 2^20 candidates; n = 8 is 2^27 and wants the C programs.
#
# the model, in the method's own terms. components 1..n, most volatile first. the upper triangular
# matrix has one spot per contiguous submixture, so a spot is an interval [i..j]. the feed [1..n]
# and the pure products [i..i] are always present and everything else is present or absent, which
# leaves (n-2)(n+1)/2 free spots, 9 at n = 5, matching "the number of submixture streams is 9 ...
# 2^9=512 candidate matrices".
#   "the only possible top products ... are those found on a right-traveling horizontal line"
#       so the top products of [i..j] are its proper prefixes [i..k], k < j
#   "the only bottom products possible ... found on a downward/right traveling along the diagonal"
#       so the bottom products of [i..j] are its proper suffixes [k..j], k > i
#   "the only feed nodes that can produce a given stream are those positioned horizontally to the
#    left or diagonally to the left and up"
#       so the feeds of [i..j] are [i..m] with m > j, and [m..j] with m < i
#
# the three checks, quoted:
#   check 1  "For every stream that exists, ensure that at least one stream that can act as its
#             feed also exists. This excludes the feed node itself."
#   check 2  "For every node with a top and bottom product, ensure that all components contained
#             in the node are also contained in one of the products."
#   check 3  "Ensure that at least n-2 submixtures are transferred."
#
# check 2 needs one reading and it does not matter which. a split of [i..j] may take any present
# proper prefix as its top and any present proper suffix as its bottom, and the pair covering the
# most is the longest of each, so "some legal split exists" is the same predicate as
# (longest prefix) + (longest suffix) >= length.
#
# check 3 turns out to remove nothing: checks 1 and 2 already force at least n-2 submixtures at
# every n tested here. the run prints that, and the controls at the end are set up for it.

TARGET = {3: 3, 4: 18, 5: 203, 6: 4373, 7: 185421}      # the published basic counts


def setup(n):
    ivs = [(i, j) for i in range(1, n + 1) for j in range(i, n + 1)]
    idx = {v: k for k, v in enumerate(ivs)}
    feed = idx[(1, n)]
    fixed = [feed] + [idx[(i, i)] for i in range(1, n + 1)]
    free = [k for k in range(len(ivs)) if k not in set(fixed)]
    pref = [[idx[(i, k)] for k in range(i, j)] for (i, j) in ivs]        # proper prefixes
    suff = [[idx[(k, j)] for k in range(i + 1, j + 1)] for (i, j) in ivs]
    feeds = [[idx[(i, m)] for m in range(j + 1, n + 1)] + [idx[(m, j)] for m in range(1, i)]
             for (i, j) in ivs]
    return ivs, idx, feed, fixed, free, pref, suff, feeds


def count(n, use1=True, use2=True, use3=True, mutate=None):
    ivs, idx, feed, fixed, free, pref, suff, feeds = setup(n)
    base = 0
    for k in fixed:
        base |= 1 << k
    L = [j - i + 1 for (i, j) in ivs]
    total = 0
    for m in range(1 << len(free)):
        present = base
        got = m
        while got:
            b = got & -got
            present |= 1 << free[b.bit_length() - 1]
            got ^= b
        if use3 and bin(m).count("1") < n - 2:
            continue
        ok = True
        for k in range(len(ivs)):
            if not (present >> k) & 1:
                continue
            if use1 and k != feed:
                if not any((present >> f) & 1 for f in feeds[k]):
                    ok = False
                    break
            if use2 and L[k] >= 2:
                p = max((L[x] for x in pref[k] if (present >> x) & 1), default=0)
                q = max((L[x] for x in suff[k] if (present >> x) & 1), default=0)
                if mutate == "drop_cover":
                    pass                      # control: skip the covering test entirely
                elif p == 0 or q == 0 or p + q < L[k]:
                    ok = False
                    break
        if ok:
            total += 1
    return total


print("n   free spots   2^free      checks 1+2+3   published   checks 1+2 only   1+3 only   2+3 only")
allok = True
for n in range(3, 8):
    free = (n - 2) * (n + 1) // 2
    c123 = count(n)
    c12 = count(n, use3=False)
    c13 = count(n, use2=False)
    c23 = count(n, use1=False)
    good = c123 == TARGET[n]
    allok &= good
    print("%d   %-10d   %-10d  %-14d %-11s %-17d %-10d %d"
          % (n, free, 1 << free, c123, "%d %s" % (TARGET[n], "ok" if good else "***"), c12, c13, c23))

print()
print("all five published terms reproduced from the printed checks:", allok)

print("\nis Check 3 doing any work?")
for n in range(3, 8):
    c123, c12 = count(n), count(n, use3=False)
    print("   n=%d  with Check 3 %-8d  without %-8d  %s"
          % (n, c123, c12, "Check 3 is REDUNDANT here" if c123 == c12 else "Check 3 removes %d" % (c12 - c123)))

# controls. checks 1 and 2 each have to be load-bearing: drop either and the published counts stop
# coming out. check 3 is the other way round, since it is implied by the first two, so dropping it
# has to change nothing. a control that expected all three to break would pass on a bug here.
print("\ncontrols:")
ok = True
pub = [TARGET[n] for n in range(3, 7)]
for label, kw, must_break in [
        ("Check 1 dropped", dict(use1=False), True),
        ("Check 2 dropped", dict(use2=False), True),
        ("Check 2 kept but its covering test gutted", dict(mutate="drop_cover"), True),
        ("Check 3 dropped (redundant, so must NOT break)", dict(use3=False), False)]:
    vals = [count(n, **kw) for n in range(3, 7)]
    broke = vals != pub
    good = broke == must_break
    ok &= good
    if must_break:
        verdict = "differs, as it must" if broke else "*** STILL MATCHES, the test is worthless"
    else:
        verdict = "unchanged, as it must be" if not broke else "*** CHANGED, Check 3 is not redundant after all"
    print("   %-46s %s  %s" % (label, vals, verdict))

print("\n" + ("ALL AS EXPECTED" if (allok and ok) else "*** SOMETHING IS NOT AS EXPECTED"))
raise SystemExit(0 if (allok and ok) else 1)
