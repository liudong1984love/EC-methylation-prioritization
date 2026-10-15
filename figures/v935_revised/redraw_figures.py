#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""确定性重绘图1与图3；不改动源文件。
依赖：Python 3、matplotlib、numpy、Pillow。
运行：python redraw_figures.py
可选：--output-dir 输出目录 --source-dir 正文目录 --dpi 600
图内按英文稿要求使用英文；SVG保留可编辑文字和矢量图元。
"""
from pathlib import Path
import argparse
import hashlib
import re
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch
from matplotlib.text import Text
import numpy as np
from PIL import Image

INK = '#203242'
MUTED = '#526573'
BLUE = '#EAF2F9'
BLUE_EDGE = '#557D9B'
GREEN_BG = '#E9F4EF'
GREEN_EDGE = '#388265'
GREY = '#F3F5F6'
AMBER_BG = '#FFF5E4'
AMBER_EDGE = '#B58637'
BAR_GREEN = '#37956A'
BAR_YELLOW = '#EAC34A'
BAR_RED = '#D55351'
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 10,
    'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': INK,
    'ytick.color': INK, 'axes.edgecolor': INK,
    'svg.fonttype': 'none', 'svg.hashsalt': 'ec-revised-figures-v1',
    'hatch.linewidth': 0.7, 'savefig.facecolor': 'white',
})

EXPECTED = [
    ('ARL5C', .041, .052), ('ZSCAN12', .044, .058),
    ('CDO1', .055, .080), ('PCDHGA', .080, .116),
    ('CYP1B1', .108, .122), ('BOLL', .141, .112),
    ('GHSR', .072, .149), ('AJAP1', .088, .295),
    ('GALR1', .226, .435), ('GYPC', .680, .838),
]


def verify_source(source):
    """从正文表6逐项读取，校验用户指定顺序中的20个数值。"""
    raw = (source / 'manuscript.txt').read_text(encoding='utf-8')
    lines = [re.sub(r'^\[P\d+\]\s*', '', s).strip() for s in raw.splitlines()]
    first = next(i for i, s in enumerate(lines) if s.startswith('Table 6.'))
    last = next(i for i, s in enumerate(lines[first + 1:], first + 1) if s.startswith('Table 7.'))
    table = lines[first:last]
    result = []
    for gene, cervix, blood in EXPECTED:
        i = table.index(gene)
        actual = (float(table[i + 1]), float(table[i + 2]))
        assert actual == (cervix, blood), f'表6数据不一致：{gene}'
        result.append((gene, *actual))
    for source_image, size in [('image1.png', (1822, 2119)), ('image3.png', (2520, 1170))]:
        with Image.open(source / source_image) as im:
            assert im.size == size, f'原图尺寸发生变化：{source_image}'
    return result


def box(ax, x, y, w, h, face=BLUE, edge=BLUE_EDGE, lw=1.15):
    p = FancyBboxPatch((x, y), w, h,
                       boxstyle='round,pad=0,rounding_size=0.9',
                       facecolor=face, edgecolor=edge, linewidth=lw, zorder=2)
    ax.add_patch(p)
    return (x, y, w, h)


def label(ax, checks, x, y, text, size=10, weight='normal', color=INK,
          ha='center', container=None, **kwargs):
    t = ax.text(x, y, text, fontsize=size, fontweight=weight,
                color=color, ha=ha, va='center', linespacing=1.35, zorder=5, **kwargs)
    checks.append((t, ax, container))
    return t


def arrow(ax, start, end, dashed=False, color=BLUE_EDGE):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle='-|>',
                 mutation_scale=11, linewidth=1.25,
                 linestyle=(0, (4, 3)) if dashed else '-', color=color,
                 shrinkA=1.5, shrinkB=1.5, zorder=1))


def check_layout(fig, checks):
    """检查全部文字是否超出画布、流程框，及显式标签是否互相碰撞。"""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    frame = fig.bbox
    for t in fig.findobj(match=Text):
        if not t.get_visible() or not t.get_text():
            continue
        b = t.get_window_extent(renderer)
        assert b.x0 >= frame.x0 - 1 and b.y0 >= frame.y0 - 1 and b.x1 <= frame.x1 + 1 and b.y1 <= frame.y1 + 1, f'文字超出画布：{t.get_text()}'
    bounds = []
    for t, ax, container in checks:
        b = t.get_window_extent(renderer)
        if container is not None:
            x, y, w, h = container
            lo = ax.transData.transform((x, y))
            hi = ax.transData.transform((x + w, y + h))
            assert b.x0 >= lo[0] + 1 and b.y0 >= lo[1] + 1 and b.x1 <= hi[0] - 1 and b.y1 <= hi[1] - 1, f'文字超出流程框：{t.get_text()}'
        bounds.append((t, b))
    for i, (t1, b1) in enumerate(bounds):
        for t2, b2 in bounds[i + 1:]:
            if b1.overlaps(b2):
                raise AssertionError(f'文字相互重叠：{t1.get_text()} / {t2.get_text()}')


def make_figure1():
    # 纵横比与原图1822×2119保持一致；全部文字采用独立可编辑对象。
    fig = plt.figure(figsize=(10.32, 10.32 * 2119 / 1822), dpi=150)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, 100), ylim=(0, 116.3))
    ax.axis('off')
    checks = []
    def text(x, y, s, **kw):
        return label(ax, checks, x, y, s, **kw)

    text(50, 113.6, 'Candidate discovery and background triage', size=18, weight='bold')
    root = box(ax, 4, 101, 92, 8)
    text(50, 106.2, 'TCGA-UCEC 450K re-analysis', size=13.2, weight='bold', container=root)
    text(50, 103.2, '23,099 DMRs  |  7,901 hypermethylated', size=12, container=root)
    arrow(ax, (23, 101), (18, 97))
    arrow(ax, (66, 101), (65.5, 97))

    dmr = box(ax, 4, 84, 28, 13, GREY, '#9AA7AE')
    text(18, 94.2, 'Parallel DMR track', size=11.4, weight='bold', container=dmr)
    text(18, 90.7, '50 candidate regions', size=12, weight='bold', container=dmr)
    text(18, 86.8, 'Cross-cohort consistency\nand regional entry criteria', size=9.2, container=dmr)
    arrow(ax, (18, 84), (18, 80), color='#89969F')
    context = box(ax, 4, 72, 28, 8, GREY, '#9AA7AE')
    text(18, 77.4, 'Regional context only', size=10.7, weight='bold', container=context)
    text(18, 74.4, 'No DMR-track region\nentered the final pool', size=9.3, container=context)

    single = box(ax, 36, 84, 59, 13)
    text(65.5, 94.8, '476 single-probe candidates', size=14, weight='bold', container=single)
    text(65.5, 89.3,
         'FDR < 0.05; Δβ ≥ 0.50; tumour positivity ≥ 90%\n'
         'Premenopausal positivity ≥ 85%; adjacent-normal P95 β ≤ 0.40\n'
         'Concordant GSE67116 direction; present in GSE73949',
         size=9.4, container=single)
    arrow(ax, (65.5, 84), (65.5, 80))
    healthy = box(ax, 36, 72, 59, 8)
    text(65.5, 77.1, '65 pass the healthy-endometrium gate', size=12.4, weight='bold', container=healthy)
    text(65.5, 74.1, 'Healthy premenopausal endometrium: P95 β ≤ 0.15', size=9.8, container=healthy)
    arrow(ax, (65.5, 72), (65.5, 68))
    text(48.4, 70.0, 'Add cervix + blood', size=8.6, color=MUTED)
    three = box(ax, 36, 60, 59, 8)
    text(65.5, 65.1, '15 pass three background dimensions', size=12.3, weight='bold', container=three)
    text(65.5, 62.0, 'Healthy endometrium + HPV-negative cervix + blood', size=9.4, container=three)
    arrow(ax, (65.5, 60), (65.5, 55))
    text(65.5, 57.5, 'Apply benign-surgical gate where assessable', size=9,
         bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': 1.5})

    rules = box(ax, 4, 40, 28, 27, GREY, '#C0CBD1')
    text(18, 63.8, 'STRICT BACKGROUND RULE', size=9.5, weight='bold', container=rules)
    text(18, 58.8, 'P95 β ≤ 0.15 in every\nassessable dimension', size=10, container=rules)
    text(18, 51.7, 'n/a = not assessable;\nnot a pass or a fail', size=9.7, container=rules)
    text(18, 44.5, 'Menstrual-cycle background:\ndescriptive only, no gate', size=9.1, container=rules)

    tier = box(ax, 36, 40, 59, 15, GREEN_BG, GREEN_EDGE, 1.45)
    text(65.5, 52.5, '2 Tier-A survivors', size=14.5, weight='bold', container=tier)
    text(65.5, 49.8, 'All assessable gates passed by point estimate', size=10, container=tier)
    text(65.5, 47.0, 'BOLL: 4/4 — borderline cervical gate', size=10.5, weight='bold', container=tier)
    text(65.5, 44.4, 'ZSCAN12: 3/3; benign-surgical n/a (absent from EPIC)', size=9.5, container=tier)
    text(65.5, 41.8, 'Score ranks Tier-A survivors only', size=10.4, weight='bold', color=GREEN_EDGE, container=tier)
    arrow(ax, (65.5, 40), (65.5, 35), color=GREEN_EDGE)
    text(48.5, 37.5, '2 strict-gate leads', size=9, color=GREEN_EDGE)

    # 虚线从476候选直接绕过严格门槛，进入异质池；不连接15通过者。
    ax.plot([95, 98.4, 98.4], [90, 90, 31.8], color=AMBER_EDGE,
            linewidth=1.3, linestyle=(0, (4, 3)), zorder=1)
    arrow(ax, (98.4, 31.8), (96, 31.8), dashed=True, color=AMBER_EDGE)
    text(97.0, 63.5, 'Different retention paths: Tier B / partner / region coverage',
         size=8.7, rotation=90, color='#896222',
         bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': 1.3})

    pool = box(ax, 4, 17, 92, 18, '#FAFCFD', BLUE_EDGE, 1.45)
    text(50, 32.7, '23-probe experimental pool  |  2 + 5 + 1 + 15', size=14, weight='bold', container=pool)
    text(50, 30.2, 'Heterogeneous retention — not an all-gates-passed set', size=10.4, color=MUTED, container=pool)
    cards = [(6, 21, GREEN_BG, GREEN_EDGE), (28.3, 21.4, AMBER_BG, AMBER_EDGE),
             (51, 21, AMBER_BG, AMBER_EDGE), (73.3, 20.7, AMBER_BG, AMBER_EDGE)]
    contents = [
        ('2 Tier-A leads', 'All assessable\nbackground gates passed'),
        ('5 Tier-B comparators', 'Specimen-relevant background\nor clinical evidence,\ndespite gate failures'),
        ('1 BOLL partner', 'Regional characterisation;\nnot prespecified\nfor positivity calls'),
        ('15 region-coverage probes', 'Coverage assembly /\nregional characterisation;\nnot strict-gate compliance'),
    ]
    for (x, w, face, edge), (heading, body) in zip(cards, contents):
        card = box(ax, x, 19.3, w, 9, face, edge, .9)
        text(x + w/2, 26.4, heading, size=9.05, weight='bold', container=card)
        text(x + w/2, 22.8, body, size=8.45, container=card)
    text(50, 18.1, 'Tier B, partner and region-coverage probes were retained on separate grounds, not by score.', size=8.4, color=MUTED, container=pool)
    arrow(ax, (50, 17), (50, 13), color=BLUE_EDGE)
    text(70, 15.0, 'Staged wet-lab evaluation', size=9, color=MUTED)

    core = box(ax, 4, 3, 92, 10)
    text(50, 10.9, 'First-stage core set: 4 CpGs across 3 genes', size=12.6, weight='bold', container=core)
    text(50, 8.2, 'ZSCAN12  cg27577527     |     CDO1  cg23180938', size=11.2, weight='bold', container=core)
    text(50, 5.4, 'BOLL  cg24589459 (anchor) + cg07495363 (partner)', size=11.1, weight='bold', container=core)
    text(50, 1.3, 'Sources: Table 2; Supplementary Notes S2.3, S2.6 and S2.14. All 23 pool probes entered via the single-probe track.', size=8.2, color=MUTED)
    return fig, checks


def category(v):
    return BAR_GREEN if v <= .15 else BAR_YELLOW if v <= .30 else BAR_RED


def make_figure3(data):
    fig = plt.figure(figsize=(14, 6.5), dpi=150)
    checks = []
    title = fig.text(.067, .941, 'Reference methylation background by compartment',
                     fontsize=17, weight='bold', va='center')
    checks.append((title, None, None))
    series = [Patch(facecolor='#B8C3CC', edgecolor=INK, label='Cervix (solid)'),
              Patch(facecolor='#B8C3CC', edgecolor=INK, hatch='////', label='Blood (hatched)')]
    fig.legend(handles=series, loc='center left', bbox_to_anchor=(.063, .866),
               frameon=False, ncol=2, fontsize=10.6, handlelength=2.0, columnspacing=1.6)
    categories = [Patch(facecolor=BAR_GREEN, label='≤ 0.15'),
                  Patch(facecolor=BAR_YELLOW, label='> 0.15 to ≤ 0.30'),
                  Patch(facecolor=BAR_RED, label='> 0.30')]
    fig.legend(handles=categories, title='Colour: P95 β category (per bar)',
               loc='center right', bbox_to_anchor=(.942, .866), frameon=False,
               ncol=3, fontsize=10.2, title_fontsize=10.2, columnspacing=1.6)
    ax = fig.add_axes([.075, .25, .855, .53])
    ax.set_ylim(0, .94)
    ax.set_xlim(-.65, 9.65)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color='#E6EBEE', linewidth=.7)
    ax.set_yticks(np.arange(0, 1, .2))
    ax.set_ylabel('Background P95 β', fontsize=12, labelpad=10)
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(.85)
    ax.tick_params(axis='both', labelsize=10.5, length=3, pad=6)
    xs = np.arange(len(data))
    width = .32
    for j, hatch in enumerate(['', '////']):
        values = [row[j + 1] for row in data]
        pos = xs + (j - .5) * .37
        bars = ax.bar(pos, values, width=width,
                      color=[category(v) for v in values], edgecolor=INK,
                      linewidth=.75, hatch=hatch, zorder=3)
        for i, (bar, v) in enumerate(zip(bars, values)):
            starred = (data[i][0] == 'BOLL' and j == 0) or (data[i][0] == 'GHSR' and j == 1)
            s = f'{v:.3f}' + ('*' if starred else '')
            t = ax.text(bar.get_x() + width/2, v + .015, s,
                        rotation=90 if v <= .15 else 0,
                        fontsize=9.4, weight='bold' if starred else 'normal',
                        ha='center', va='bottom', zorder=6,
                        bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': .4})
            checks.append((t, ax, None))
    for value, linestyle, lab in [(.15, (0, (5, 3)), '0.15 gate'), (.30, (0, (1.5, 2.5)), '0.30')]:
        ax.axhline(value, color='#586976', linewidth=1.1, linestyle=linestyle, zorder=4)
        t = ax.text(1.008, value, lab, transform=ax.get_yaxis_transform(),
                    fontsize=9, color=MUTED, va='center', ha='left')
        checks.append((t, ax, None))
    ax.set_xticks(xs, [row[0] for row in data], rotation=30, ha='right', rotation_mode='anchor')
    footnotes = [
        (.104, '* Point-estimate passes with bootstrap 95% CIs crossing 0.15: BOLL Cervix, 0.100–0.158; GHSR Blood, 0.101–0.160.'),
        (.065, 'Source: Table 6 (P95 values and colour categories); CI details: main text and Supplementary Note S2.7.'),
        (.029, 'Reference-background P95 estimates only — not diagnostic performance in actual sampled specimens.'),
    ]
    for y, s in footnotes:
        t = fig.text(.075, y, s, fontsize=9.5, color=MUTED, va='center')
        checks.append((t, None, None))
    # 图3精确核验：20根柱、顺序、20个三位小数标签及两个星号。
    assert len(ax.patches) == 20
    numbers = [t.get_text() for t, _, _ in checks if re.fullmatch(r'0\.\d{3}\*?', t.get_text())]
    assert len(numbers) == 20 and sum(s.endswith('*') for s in numbers) == 2
    return fig, checks


def save_figure(fig, checks, target, dpi):
    check_layout(fig, checks)
    fig.set_dpi(dpi)
    check_layout(fig, checks)
    # 不使用tight裁切，从而保留原纵横比和预留边距。
    fig.savefig(target.with_suffix('.png'), dpi=dpi)
    fig.savefig(target.with_suffix('.svg'), metadata={'Date': None})
    with Image.open(target.with_suffix('.png')) as im:
        im.verify()
    with Image.open(target.with_suffix('.png')) as im:
        print(f'已生成：{target.with_suffix(".png")}；{im.size[0]}×{im.size[1]}像素；{dpi} dpi')
    tree = ET.parse(target.with_suffix('.svg'))
    assert tree.findall('.//{http://www.w3.org/2000/svg}text'), 'SVG缺少可编辑文字'
    assert not tree.findall('.//{http://www.w3.org/2000/svg}image'), 'SVG不应包含栅格图片'
    print(f'已生成：{target.with_suffix(".svg")}；可编辑文字和矢量图元；文字边界检查通过')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description='重绘医学稿图1、图3，不修改正文或Word。')
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--source-dir', type=Path, default=None, help='可选：原始抽取文本和图片目录；不指定时使用已核对的表6数值。')
    parser.add_argument('--dpi', type=int, default=600)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    originals = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in args.source_dir.iterdir() if p.is_file()} if args.source_dir else {}
    if args.source_dir:
        data = verify_source(args.source_dir)
        print('表6全部20个数值及原图尺寸校验通过。')
    else:
        data = list(EXPECTED)
        print('使用随程序保存、在本次编辑时已与表6逐项核对的20个数值；未重新读取原始数据。')
    fig, checks = make_figure1()
    save_figure(fig, checks, args.output_dir / 'Figure1_revised', args.dpi)
    fig, checks = make_figure3(data)
    save_figure(fig, checks, args.output_dir / 'Figure3_revised', args.dpi)
    assert all(hashlib.sha256(p.read_bytes()).hexdigest() == digest for p, digest in originals.items()), '源文件发生变化'
    print('校验完成：源目录所有文件内容未改动。')


if __name__ == '__main__':
    main()
