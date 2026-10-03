import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(11.5, 13.5), dpi=200)
ax.set_xlim(0, 100); ax.set_ylim(0, 135); ax.axis('off')

C_MAIN   = '#DCE9F7'
C_SIDE   = '#F2F2F2'
C_GATE   = '#FDEBD0'
C_ELIM   = '#FADBD8'
C_FINAL  = '#D5F5E3'
EDGE     = '#2C3E50'
MARGIN   = 2.0
GAP      = 4.5          # vertical gap between stacked funnel boxes

def measure(text, fs, bold):
    t = ax.text(0, 0, text, ha='center', va='center', fontsize=fs,
                fontweight='bold' if bold else 'normal')
    fig.canvas.draw()
    bb = t.get_window_extent(fig.canvas.get_renderer())
    inv = ax.transData.inverted()
    (x0, y0), (x1, y1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
    t.remove()
    return (x1 - x0), (y1 - y0)

def draw_box(x, y, w, h, text, fc, fs, bold):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6,rounding_size=1.2",
                       linewidth=1.4, edgecolor=EDGE, facecolor=fc, zorder=2)
    ax.add_patch(p)
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=fs,
            fontweight='bold' if bold else 'normal', color='#1B2631', zorder=3)
    return (x, y, w, h)

def tbox(cx, cy, text, fc, fs=10.5, bold=False, pad_x=3.2, pad_y=1.8):
    tw, th = measure(text, fs, bold)
    w = tw + 2 * pad_x
    h = th + 2 * pad_y
    x = cx - w / 2
    y = cy - h / 2
    if x < MARGIN: x = MARGIN
    if x + w > 100 - MARGIN: x = 100 - MARGIN - w
    if y < MARGIN: y = MARGIN
    if y + h > 135 - MARGIN: y = 135 - MARGIN - h
    return draw_box(x, y, w, h, text, fc, fs, bold)

def tbox_top(cx, top_y, text, fc, fs=10.5, bold=False, pad_x=3.2, pad_y=1.8):
    """Box whose TOP edge sits at top_y (sequential stacking)."""
    tw, th = measure(text, fs, bold)
    w = tw + 2 * pad_x
    h = th + 2 * pad_y
    x = cx - w / 2
    y = top_y - h
    if x < MARGIN: x = MARGIN
    if x + w > 100 - MARGIN: x = 100 - MARGIN - w
    return draw_box(x, y, w, h, text, fc, fs, bold)

def cxb(b): return b[0] + b[2] / 2
def top(b): return b[1] + b[3]
def bot(b): return b[1]

def arrow(x1, y1, x2, y2, color=EDGE, lw=1.6, ls='-'):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>', mutation_scale=16,
                                 linewidth=lw, color=color, linestyle=ls, zorder=1))

ax.text(50, 132.5, 'Candidate discovery and background-triage funnel', ha='center',
        fontsize=14, fontweight='bold', color='#1B2631')

# Top source
b0 = tbox(50, 126.5, 'TCGA-UCEC 450K re-analysis\n23,099 DMRs (7,901 hypermethylated)', C_MAIN, bold=True)

# Two tracks
b1 = tbox(24, 112, 'DMR track (parallel)\n50 candidate regions\ncross-cohort consistency required\n(same-direction Δβ in GSE67116)', C_SIDE)
b2 = tbox(74, 112, 'Single-probe track (candidate-generating)\n476 elite probes\nFDR<0.05, Δβ≥0.50, tumour positivity ≥90%,\npremenopausal positivity ≥85%', C_MAIN, bold=True)
arrow(44, bot(b0), cxb(b1), top(b1)); arrow(56, bot(b0), cxb(b2), top(b2))

# DMR outcome note
b3 = tbox(24, 99.5, 'Regional context and auxiliary evidence only\nno DMR-track region carried into the final pool', C_SIDE, fs=9.5)
arrow(cxb(b1), bot(b1), cxb(b3), top(b3), ls='--')

# Funnel: gates -> elimination -> score, stacked with explicit GAP so arrows show
b4 = tbox_top(64, 94.5, 'Background triage — four P95 gates (β ≤ 0.15, union rule)\nhealthy premenopausal endometrium · benign surgical endometrium\nHPV-negative cervix · blood-reference fractions\nn/a = dimension not assessable (missing, not failed)\n+ menstrual-cycle dimension: descriptive, no gate', C_GATE, fs=10)
arrow(cxb(b2), bot(b2), cxb(b4), top(b4))
# dashed "context" feed from DMR note into the gate box (corner to corner, no overlap)
arrow(b3[0] + b3[2], bot(b3), b4[0] + 2.5, top(b4), ls='--', color='#7F8C8D')
ax.text((b3[0] + b3[2] + b4[0]) / 2 + 1.5, (bot(b3) + top(b4)) / 2 + 1.0, 'context',
        fontsize=8.5, color='#7F8C8D', ha='center')

b5 = tbox_top(64, bot(b4) - GAP, '461/476 eliminated by background triage\n(65/476 pass healthy-endometrium screen;\n15/65 pass four-gate triage — scoring of survivors only)', C_ELIM, fs=10)
arrow(cxb(b4), bot(b4), cxb(b5), top(b5))

b6 = tbox_top(64, bot(b5) - GAP, 'Three-term effective score\n0.25·s_bg + 0.25·premenopausal positivity + 0.20·overall positivity\n(two constant historical terms cannot affect ranking)', C_MAIN, fs=10)
arrow(cxb(b5), bot(b5), cxb(b6), top(b6))

# Final pool
b7 = tbox_top(50, bot(b6) - GAP - 1.5, '23-probe experimental pool\nheterogeneous library — three retention rationales', C_FINAL, bold=True, fs=11)
arrow(cxb(b6), bot(b6), cxb(b7), top(b7))

# Four branches — auto layout across the full width with font-size fallback
branch_texts = [
    'Tier A leads (2)\nBOLL anchor: 4/4 gates\nZSCAN12 anchor: 3/3\n(benign-surgical n/a)',
    'Tier B comparators (5)\nretained on specimen-relevant\nbackground or clinical evidence\ndespite gate failures',
    'BOLL partner probe (1)\nregional characterisation\n(not prespecified for\npositivity calls)',
    'Region-coverage probes (15)\ngreedy incremental coverage;\ndescriptive post-assembly audit\n(14/15 cervix+blood ≤0.30)',
]
fs_try = 9.5
while True:
    widths = [measure(t, fs_try, False)[0] + 2 * 3.0 for t in branch_texts]
    total = sum(widths)
    if total + 3 * 1.5 <= 100 - 2 * MARGIN or fs_try <= 7.5:
        break
    fs_try -= 0.5
gap = (100 - 2 * MARGIN - total) / 3
boxes = []
x_cur = MARGIN
branch_top = bot(b7) - 4.0
for t, w in zip(branch_texts, widths):
    tw, th = measure(t, fs_try, False)
    h = th + 2 * 1.8
    y = branch_top - h
    p = FancyBboxPatch((x_cur, y), w, h, boxstyle="round,pad=0.6,rounding_size=1.2",
                       linewidth=1.4, edgecolor=EDGE, facecolor=C_FINAL, zorder=2)
    ax.add_patch(p)
    ax.text(x_cur + w / 2, y + h / 2, t, ha='center', va='center', fontsize=fs_try,
            color='#1B2631', zorder=3)
    boxes.append((x_cur, y, w, h))
    x_cur += w + gap
c1, c2, c3, c4 = boxes
for c in boxes:
    arrow(cxb(b7), bot(b7), cxb(c), top(c))

# Core set
b8 = tbox(50, 13.2, 'First-stage core set — 4 CpGs across 3 genes\nZSCAN12 (cg27577527) · CDO1 (cg23180938)\nBOLL region (anchor cg24589459 + partner cg07495363)\nexpansion into the pool conditional on core-set performance', C_MAIN, bold=True, fs=10.5)
for c in boxes:
    arrow(cxb(c), bot(c), cxb(b8), top(b8))

ax.text(50, 2.4, 'Counts as in Table 2. All 23 experimental-pool probes entered via the single-probe track;\nthe DMR track provided regional context but contributed no candidate to the final pool.',
        ha='center', fontsize=9, color='#566573', style='italic')

out = r'C:\Users\ld\WorkBuddy\2026-09-25-20-09-18\fig1_v6.png'
plt.savefig(out, bbox_inches='tight', facecolor='white', dpi=200)
print('saved, branch fs =', fs_try)
for name, b in [('b0', b0), ('b1', b1), ('b2', b2), ('b3', b3), ('b4', b4), ('b5', b5),
                ('b6', b6), ('b7', b7), ('c1', c1), ('c2', c2), ('c3', c3), ('c4', c4), ('b8', b8)]:
    ok = b[0] >= MARGIN - 0.01 and b[0] + b[2] <= 100 - MARGIN + 0.01 and b[1] >= 0
    print(f'{name}: x=[{b[0]:.1f},{b[0]+b[2]:.1f}] y=[{b[1]:.1f},{b[1]+b[3]:.1f}] {"OK" if ok else "CLIPPED!"}')
# gaps between the three funnel boxes
print('gap b4->b5:', bot(b4) - top(b5))
print('gap b5->b6:', bot(b5) - top(b6))
print('gap b6->b7:', bot(b6) - top(b7))
