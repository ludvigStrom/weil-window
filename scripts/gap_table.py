"""Match the lowest ordinary eigenstates of the odd half of the Weil form to the gaps between zeta zeros.
Reads flat_track.json (written by flat_track.py) and zeta_zeros_to_2300.json (written by zeros_rs.py) from the current directory.
For each window c: the three lowest ordinary states, the two zeros on either side of the centre of each state's frequency band,
the gap between them, and that gap's rank among all gaps between consecutive zeros lying above 2 pi c and within the basis.
A state counts as a newborn tiny level, and is skipped, if its eigenvalue is below 0.1 and more than 2% of its weight lies
within 0.02 of the window's ends.
Usage: python3 gap_table.py"""
import json, math
import numpy as np
zeros = np.array(json.load(open("zeta_zeros_to_2300.json"))); gaps = np.diff(zeros)
best = {}
for r in json.load(open("flat_track.json")):
    if r["c"] not in best or r["ntiny"] < best[r["c"]]["ntiny"]: best[r["c"]] = r      # prefer the run with the lower threshold
print("    c |  2 pi c | 2/sqrt(c) | eigenvalue | edge share | zeros on either side | gap  | rank above 2 pi c")
for c in sorted(best):
    r = best[c]; W = r["W"]; top = 2 * math.pi * r["N"] / r["L"]; usable = (zeros[:-1] > W) & (zeros[1:] < 0.9 * top); shown = 0
    for s in r["states"]:
        if s["eps"] < 0.1 and s["edge"] > 0.02: continue
        mid = s["band"][1]; i = int(np.searchsorted(zeros, mid)) - 1
        if i < 0 or i >= len(gaps): continue
        rank = int((np.where(usable, gaps, 0) > gaps[i]).sum()) + 1 if usable[i] else None
        where = f"{rank}" if rank else ("below 2 pi c" if zeros[i] < W else "near the top of the basis")
        print(f"{c:5.1f} | {W:6.0f} | {2 / math.sqrt(c):9.3f} | {s['eps']:10.4f} | {s['edge']:10.4f} | {zeros[i]:8.2f}, {zeros[i + 1]:8.2f}   | {gaps[i]:.2f} | {where}")
        shown += 1
        if shown == 3: break
