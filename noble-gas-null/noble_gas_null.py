"""Null forensic on Paper 17 v1.0 noble gas fitting.

Question: in the integer-combination space of {V=12, F=20, E=30, χ=2, D=3, D!=6}
under arithmetic operations {+, -, *, **2}, how many distinct positive integers ≤ 100
can be generated? Compare to 6 noble gas atomic numbers {2, 10, 18, 36, 54, 86}.

If the expressible space is much larger than the number of noble gases, then hitting
all six is cheap and the v1.0 "match" is overparameterised fitting.
"""

ATOMS = {'V': 12, 'F': 20, 'E': 30, 'χ': 2, 'D': 3, 'D!': 6}
NOBLE_GASES = [2, 10, 18, 36, 54, 86]
MAX_VALUE = 100

def grow(current, ops_allowed):
    """One level of growth: apply unary (square) and binary (+, -, *) ops."""
    new = set(current)
    for a in current:
        if a <= 12:  # cap squaring to avoid blowup
            sq = a * a
            if 1 <= sq <= MAX_VALUE:
                new.add(sq)
    for a in current:
        for b in current:
            for op in ops_allowed:
                if op == '+':
                    v = a + b
                elif op == '-':
                    v = a - b
                elif op == '*':
                    v = a * b
                if 1 <= v <= MAX_VALUE:
                    new.add(v)
    return new

# Level 0: atoms
level = set(ATOMS.values())
print(f"Level 0 (atoms): {sorted(level)} ({len(level)} values)")

# Grow through several levels
for d in range(1, 4):
    level = grow(level, ['+', '-', '*'])
    in_range = sorted(v for v in level if 1 <= v <= MAX_VALUE)
    hits = [z for z in NOBLE_GASES if z in level]
    print(f"Level {d}: {len(in_range)} distinct positive integers ≤ {MAX_VALUE}")
    print(f"         Noble gas hits: {hits}")

# At final level, report statistics
print(f"\n=== FINAL (depth 3) ===")
expressible = sorted(v for v in level if 1 <= v <= MAX_VALUE)
print(f"Expressible integers ≤ {MAX_VALUE}: {len(expressible)} / {MAX_VALUE} = {100*len(expressible)/MAX_VALUE:.1f}%")
print(f"Noble gases hit: {[z for z in NOBLE_GASES if z in level]} ({sum(1 for z in NOBLE_GASES if z in level)}/6)")

# Bonus: what's NOT expressible?
not_expressible = [v for v in range(1, MAX_VALUE+1) if v not in level]
print(f"\nIntegers ≤ {MAX_VALUE} NOT expressible: {not_expressible}")
print(f"Count: {len(not_expressible)}")

# Chance baseline: if expressible space is N integers within ≤86,
# probability a specific integer is among them is N/86
expressible_le_86 = [v for v in expressible if v <= 86]
print(f"\nExpressible integers ≤ 86: {len(expressible_le_86)}/86 = {100*len(expressible_le_86)/86:.1f}%")
print(f"P(specific Z hit | uniform): {len(expressible_le_86)/86:.3f}")
print(f"P(all 6 noble gases hit | uniform & independent): {(len(expressible_le_86)/86)**6:.3f}")
print(f"\nBUT: noble gases are NOT randomly distributed — they're highly composite, often near factorials/squares.")
print(f"The cheap-fitting diagnosis stands: ~{100*len(expressible_le_86)/86:.0f}% of integers ≤ 86 are expressible.")
