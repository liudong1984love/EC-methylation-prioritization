# fix_zh_v932.py — Chinese manuscript v9.31 -> v9.32 (mirror of fix_ms_v932)
import zipfile, os

SRC = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Manuscript_EC_methylation_panel_v9.31_中文全文翻译.docx'
DST = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\Manuscript_EC_methylation_panel_v9.32_中文全文翻译.docx'
FIG4 = r'C:\Users\ld\WorkBuddy\2026-08-19-15-50-19\results\v932_figure4_crosscohort.png'

zin = zipfile.ZipFile(SRC)
xml = zin.read('word/document.xml').decode('utf-8')

MD = '\u2014'

REPL = [
 # 1) Methods: donor-level blood resampling
 ('GSE35069 在细胞类型组分内）',
  'GSE35069 按供者水平：六名经核实的供者有放回重抽样，每名被抽中的供者贡献全部十种血液参考组分'
  + MD + '逐组分 P95 再取最大值' + MD + '保留供者跨组分相关结构；供者身份已对照 ArrayExpress '
  '样本元数据核实，补充注 S2.18）'),
 # 2) Methods: cross-cohort consistency correction
 ('且单探针路径的纳入不要求跨队列一致性（补充数据 S2b）',
  '且单探针路径的纳入还要求跨队列方向一致性（GSE67116 中肿瘤' + MD + '增生 Δβ>0 方向一致）'
  '并存在于 GSE73949 noob 矩阵（补充数据 S2b）'),
 # 3) Results: v3 stability numbers
 ('分别在 53.1% 和 48.3% 的重复中重新进入 A 级',
  '分别在 52.8% 和 48.7% 的重复中重新进入 A 级'),
 ('跌破 ≥0.90 的发现标准（45.1% 的重复）',
  '跌破 ≥0.90 的发现标准（45.9% 的重复）'),
 ('23.1% 的重复中总体阳性率低于 ≥0.90（观测值 0.909），16.1% 中绝经前阳性率低于 ≥0.85，'
  '28.4% 中临界的宫颈门槛不通过（71.6% 合规）',
  '21.7% 的重复中总体阳性率低于 ≥0.90（观测值 0.909），16.8% 中绝经前阳性率低于 ≥0.85，'
  '26.6% 中临界的宫颈门槛不通过（73.5% 合规）'),
 ('cg27577527 在全部 1,061 次保留重复中排名 1',
  'cg27577527 在全部 1,055 次保留重复中排名 1'),
 ('cg24589459 在 99.8% 中排名 ≤5',
  'cg24589459 在 99.7% 中排名 ≤5'),
 # Results stability preamble: donor-level note
 ('（2,000 次嵌套 bootstrap 重复；绝经前肿瘤与其他肿瘤作为独立分层分别重抽样后合并，'
  '使两个阳性率指标来自同一次抽取；补充注 S2.9）',
  '（2,000 次嵌套 bootstrap 重复；绝经前肿瘤与其他肿瘤作为独立分层分别重抽样后合并，'
  '使两个阳性率指标来自同一次抽取，血液参考组分按供者水平重抽样；补充注 S2.9）'),
 # 4) purified blood cells
 ('宫颈刮片和纯化血细胞）',
  '宫颈刮片和血液参考组分）'),
 # 5) Abstract OR rule
 ('简单的双探针 OR 组合在良性手术子宫内膜中频繁阳性',
  'BOLL 区域双探针的 OR 规则在良性手术子宫内膜中频繁阳性'),
 # 6) Intro softening
 ('由于没有任何公共队列与该人群匹配',
  '由于我们的检索未发现任何符合条件的队列与该人群匹配'),
 # 7) BOLL expression wording
 ('无论甲基化状态如何，BOLL 表达在子宫内膜组织中基本缺失',
  '无论甲基化状态如何，BOLL 表达在大多数子宫内膜组织样本中很低'),
 # 8) CDO1 dual-use (Discussion)
 ('与其在宫颈（CISENDO）和宫腔内两种形式中的双重应用一致',
  '支持其宫颈形式的应用（CISENDO）；其宫腔内应用依据的是既往临床研究，因为我们的队列中'
  '其良性子宫内膜背景未达到门槛'),
 # 8b) CDO1 dual-use (Results)
 ('CDO1 显示双部位谱（宫颈 0.055、血液 0.080），与其在宫颈和宫腔内两种形式中的应用一致。',
  'CDO1 显示低宫颈（0.055）和血液（0.080）背景，支持宫颈路径应用；其宫腔内应用依据既往'
  '临床研究，而非本研究评估的良性子宫内膜背景。'),
 # 9) Limitations 15->10 simplified
 ('第二阶段分流通过集合' + MD + MD + '位于严格四门 A 级认定上游的辅助阶段，并非 A 级存活者' +
  MD + MD + '由 15 减至 10',
  '辅助性的第二阶段分流通过集合（位于严格 A 级认定的上游，并非 A 级存活者）由 15 减至 10'),
 # 10) Figure 4 caption
 ('图 4. BOLL 和 ZSCAN12 区域的跨队列表现。',
  '图 4. BOLL 和 ZSCAN12 区域的跨队列甲基化背景与肿瘤组织阳性率。'),
]

for old, new in REPL:
    n = xml.count(old)
    print(n, '|', old[:50])
    assert n == 1, f'anchor count {n} != 1: {old[:40]}'
    xml = xml.replace(old, new, 1)

for bad in ['在细胞类型组分内', '纯化血细胞', '不要求跨队列一致性', '53.1%', '48.3%', '45.1%',
            '23.1% 的重复', '1,061', '基本缺失', '双重应用一致', '双部位谱', '没有任何公共队列',
            '双探针 OR 组合', '跨队列表现。']:
    c = xml.count(bad)
    if c: print('RESIDUAL', bad, c)

fig4_bytes = open(FIG4, 'rb').read()
with zipfile.ZipFile(DST + '.tmp', 'w', zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == 'word/document.xml':
            data = xml.encode('utf-8')
        elif item.filename == 'word/media/image4.png':
            data = fig4_bytes
        zout.writestr(item, data)
zin.close()
os.replace(DST + '.tmp', DST)

import xml.etree.ElementTree as ET
z = zipfile.ZipFile(DST)
assert z.testzip() is None
ET.fromstring(z.read('word/document.xml'))
print('image4 size:', z.getinfo('word/media/image4.png').file_size)
print('OK ->', DST, os.path.getsize(DST))
