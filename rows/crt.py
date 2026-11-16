"""the thermally coupled count passes 128 bits at n = 16 (the basic one at 17), so the program keeps each
count mod 2^128 and mod 2^61 - 1 and this puts them back together (chinese remainder theorem, exact below 2^189).
python crt.py out.txt   (or pipe the program's output in)
"""
import re
import sys

P = (1 << 61) - 1
M = 1 << 128
INV = pow(M % P, P - 2, P)


def crt(r128, rp):
    t = ((rp - r128) % P) * INV % P
    return r128 + M * t


src = open(sys.argv[1]).read() if len(sys.argv) > 1 else sys.stdin.read()
for line in src.splitlines():
    m = re.match(r"RES n=(\d+) basic_m128=(\d+) basic_mP=(\d+) sum_m128=(\d+) sum_mP=(\d+)", line)
    if not m:
        continue
    n, b128, bp, s128, sp = map(int, m.groups())
    basic = crt(b128, bp)
    total = crt(s128, sp)
    print(n, basic, total - basic)
