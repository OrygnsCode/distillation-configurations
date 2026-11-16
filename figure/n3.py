# draws figure/n3.svg, the three basic configurations for n = 3, after checking them against the rules.
# standard library only. run it from the figure folder: python n3.py
from itertools import combinations

n = 3
color = {1: "#d55e00", 2: "#009e73", 3: "#0072b2"}   # A, B, C

# every set of streams that follows the two rules (see the README)
free = [(i, j) for i in range(1, n + 1) for j in range(i + 1, n + 1) if (i, j) != (1, n)]
fixed = [(1, n)] + [(c, c) for c in range(1, n + 1)]
configs = []
for k in range(len(free) + 1):
    for pick in combinations(free, k):
        s = set(fixed) | set(pick)
        ok = True
        for (i, j) in s:
            if i < j:
                top = max(b for (a, b) in s if a == i and b < j)
                bot = min(a for (a, b) in s if b == j and a > i)
                ok &= bot <= top + 1
            if (i, j) != (1, n):
                ok &= any((a == i and b > j) or (b == j and a < i) for (a, b) in s)
        if ok:
            configs.append(s)
assert len(configs) == 3

# each split: feed -> longest stream below it in its row (top) and in its column (bottom).
# splits that make a common product share a column
def columns(s):
    splits = []
    for (i, j) in sorted(s, key=lambda t: t[0] - t[1]):
        if i < j:
            top = (i, max(b for (a, b) in s if a == i and b < j))
            bot = (min(a for (a, b) in s if b == j and a > i), j)
            splits.append(((i, j), top, bot))
    cols = []
    for sp in splits:
        for c in cols:
            if any({sp[1], sp[2]} & {q[1], q[2]} for q in c):
                c.append(sp)
                break
        else:
            cols.append([sp])
    return [[sp[0] for sp in c] for c in cols]

# what gets drawn, left to right: the feeds of column 1 and of column 2
drawn = [
    [[(1, 3)], [(2, 3)]],            # A off first, then B / C
    [[(1, 3)], [(1, 2)]],            # C off first, then A / B
    [[(1, 3)], [(1, 2), (2, 3)]],    # AB and BC both into column 2, B from the middle
]
norm = lambda cfg: sorted(sorted(c) for c in cfg)
got = [columns(s) for s in configs]
assert sorted(map(norm, got)) == sorted(map(norm, drawn)), got
for d in drawn:
    assert len(d) == n - 1            # n - 1 columns

# drawing
W, H, PW = 960, 300, 320
out = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d">' % (W, H, W, H),
       '<rect width="%d" height="%d" fill="white"/>' % (W, H)]


def stream(x1, x2, y, comps):
    # one stripe per component, centered on y, with an arrowhead
    w = 6
    for k, c in enumerate(comps):
        yy = y + (k - (len(comps) - 1) / 2) * (w + 2)
        out.append('<line x1="%g" y1="%g" x2="%g" y2="%g" stroke="%s" stroke-width="%d"/>' % (x1, yy, x2 - 8, yy, color[c], w))
    h = len(comps) * (w + 2) / 2 + 3
    out.append('<polygon points="%g,%g %g,%g %g,%g" fill="#555"/>' % (x2 - 9, y - h, x2, y, x2 - 9, y + h))


def elbow(x, y1, y2, x2, c):
    # a single product leaving the top or bottom of a column and going right
    out.append('<polyline points="%g,%g %g,%g %g,%g" fill="none" stroke="%s" stroke-width="6"/>' % (x, y1, x, y2, x2 - 8, y2, color[c]))
    out.append('<polygon points="%g,%g %g,%g %g,%g" fill="#555"/>' % (x2 - 9, y2 - 6, x2, y2, x2 - 9, y2 + 6))


def column(x):
    out.append('<rect x="%g" y="50" width="44" height="200" rx="20" fill="#f2f2f2" stroke="#555" stroke-width="2"/>' % x)


comp = lambda iv: list(range(iv[0], iv[1] + 1))
for p, d in enumerate(drawn):
    ox = p * PW
    c1, c2 = ox + 80, ox + 200
    column(c1)
    column(c2)
    stream(ox + 15, c1, 150, comp((1, 3)))
    if d[1] == [(2, 3)]:
        elbow(c1 + 22, 50, 25, ox + 305, 1)
        stream(c1 + 44, c2, 210, comp((2, 3)))
        stream(c2 + 44, ox + 305, 80, [2])
        stream(c2 + 44, ox + 305, 220, [3])
    elif d[1] == [(1, 2)]:
        stream(c1 + 44, c2, 90, comp((1, 2)))
        elbow(c1 + 22, 250, 275, ox + 305, 3)
        stream(c2 + 44, ox + 305, 80, [1])
        stream(c2 + 44, ox + 305, 220, [2])
    else:
        stream(c1 + 44, c2, 90, comp((1, 2)))
        stream(c1 + 44, c2, 210, comp((2, 3)))
        stream(c2 + 44, ox + 305, 70, [1])
        stream(c2 + 44, ox + 305, 150, [2])
        stream(c2 + 44, ox + 305, 230, [3])
out.append('</svg>')
open("n3.svg", "w", newline="\n").write("\n".join(out) + "\n")
print("checked: 3 configurations, n - 1 = 2 columns each; wrote n3.svg")
