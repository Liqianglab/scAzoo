from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

ROOT = Path('/mnt/data/supfig4_bubble_redraw')
DATA_IN = ROOT / 'input' / 'supfig4_source_data'
OUT = Path('/mnt/data/SupFig4_redrawn_ST_pathway_bubbles')
SRC = OUT / 'supfig4_source_data'
CODE = OUT / 'supfig4_code'
PDF = OUT / 'panels_pdf'
SVG = OUT / 'panels_svg'
PNG = OUT / 'panels_png'
PREV = OUT / 'preview'
for d in [SRC, CODE, PDF, SVG, PNG, PREV]:
    d.mkdir(parents=True, exist_ok=True)

group_order = ['AZFc_Del', 'iNOA_B', 'iNOA_S', 'KS']

# Pathways selected to keep each panel readable while preserving the main GSEA signal.
stage_a_paths = [
    'Chemokine Signaling Pathway',
    'Endocytosis',
    'Focal Adhesion',
    'Mapk Signaling Pathway',
    'Regulation Of Actin Cytoskeleton',
    'Wnt Signaling Pathway',
    'Insulin Signaling Pathway',
    'Tight Junction',
    'Oxidative Phosphorylation',
    'Spliceosome',
]
stage_c_paths = [
    'Apoptosis',
    'Cell Adhesion Molecules Cams',
    'Chemokine Signaling Pathway',
    'Endocytosis',
    'Insulin Signaling Pathway',
    'Mapk Signaling Pathway',
    'Oxidative Phosphorylation',
    'Regulation Of Actin Cytoskeleton',
    'Spliceosome',
    'Wnt Signaling Pathway',
]

label_map = {
    'Chemokine Signaling Pathway': 'Chemokine signaling',
    'Endocytosis': 'Endocytosis',
    'Focal Adhesion': 'Focal adhesion',
    'Mapk Signaling Pathway': 'MAPK signaling',
    'Regulation Of Actin Cytoskeleton': 'Actin cytoskeleton regulation',
    'Wnt Signaling Pathway': 'Wnt signaling',
    'Insulin Signaling Pathway': 'Insulin signaling',
    'Tight Junction': 'Tight junction',
    'Oxidative Phosphorylation': 'Oxidative phosphorylation',
    'Spliceosome': 'Spliceosome',
    'Apoptosis': 'Apoptosis',
    'Cell Adhesion Molecules Cams': 'Cell adhesion molecules',
}


def load_panel(filename, selected_paths):
    df = pd.read_csv(DATA_IN / filename)
    # Keep only terms selected for the clean plot. Deduplicate if needed, keeping the lowest adjusted P value.
    df = df[df['pathway_readable'].isin(selected_paths)].copy()
    df['group'] = pd.Categorical(df['group'], categories=group_order, ordered=True)
    df['pathway_readable'] = pd.Categorical(df['pathway_readable'], categories=selected_paths, ordered=True)
    df = df.sort_values(['pathway_readable', 'group', 'p.adjust'])
    df = df.drop_duplicates(['pathway_readable', 'group'], keep='first')
    df['pathway_label'] = df['pathway_readable'].map(label_map).fillna(df['pathway_readable'].astype(str))
    # Save selected long data in plotting order.
    return df

stage_a = load_panel('SupFig4D_ST1_disease_vs_ctrl_top10.csv', stage_a_paths)
stage_c = load_panel('SupFig4F_ST3_disease_vs_ctrl_top10.csv', stage_c_paths)

# Also keep the original source tables for traceability.
for fn in ['SupFig4D_ST1_disease_vs_ctrl_top10.csv', 'SupFig4D_ST1_disease_vs_ctrl_gsea_all.csv',
           'SupFig4F_ST3_disease_vs_ctrl_top10.csv', 'SupFig4F_ST3_disease_vs_ctrl_gsea_all.csv']:
    shutil.copy2(DATA_IN / fn, SRC / fn)
stage_a.to_csv(SRC / 'SupFig4C_Stage_a_immature_stress_selected_bubble_data.csv', index=False)
stage_c.to_csv(SRC / 'SupFig4D_Stage_c_mature_supportive_selected_bubble_data.csv', index=False)

all_nes = pd.concat([stage_a['NES'], stage_c['NES']], ignore_index=True)
nes_lim = max(3.0, float(np.ceil(np.nanmax(np.abs(all_nes)))))
norm = TwoSlopeNorm(vmin=-nes_lim, vcenter=0, vmax=nes_lim)
cmap = plt.get_cmap('RdBu_r')
all_neglog = pd.concat([stage_a['neglog10_fdr'], stage_c['neglog10_fdr']], ignore_index=True)
neg_min = float(np.nanmin(all_neglog))
neg_max = float(np.nanmax(all_neglog))
if neg_max == neg_min:
    neg_max += 0.1


def size_from_fdr(x):
    # Matplotlib scatter area in pt^2. Keep differences visible even when FDR range is narrow.
    return 45 + (np.asarray(x, dtype=float) - neg_min) / (neg_max - neg_min) * 240


def plot_bubble(df, selected_paths, title, panel_label, out_stem, fig_width=8.2, fig_height=4.2):
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    x_map = {g: i for i, g in enumerate(group_order)}
    # Put first selected pathway at top.
    y_order = selected_paths[::-1]
    y_map = {p: i for i, p in enumerate(y_order)}
    for _, r in df.iterrows():
        x = x_map[str(r['group'])]
        y = y_map[str(r['pathway_readable'])]
        ax.scatter(x, y, s=size_from_fdr(r['neglog10_fdr']),
                   c=[cmap(norm(r['NES']))], edgecolor='0.25', linewidth=0.35, alpha=0.95)
    ax.set_xlim(-0.55, len(group_order)-0.45)
    ax.set_ylim(-0.6, len(y_order)-0.4)
    ax.set_xticks(range(len(group_order)))
    ax.set_xticklabels(group_order, rotation=45, ha='right', fontsize=8)
    ax.set_yticks(range(len(y_order)))
    ax.set_yticklabels([label_map.get(p, p) for p in y_order], fontsize=8)
    ax.grid(True, color='0.88', linewidth=0.7)
    ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_linewidth(0.8)
    ax.set_title(title, fontsize=10, weight='bold', pad=7)
    ax.text(-0.22, 1.10, panel_label, transform=ax.transAxes, fontsize=14, weight='bold')
    # Colorbar
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.025)
    cbar.set_label('NES', fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    # Size legend
    legend_vals = [round(neg_min, 1), round((neg_min + neg_max)/2, 1), round(neg_max, 1)]
    # Avoid duplicates
    legend_vals = sorted(set(legend_vals))
    handles = [plt.scatter([], [], s=size_from_fdr(v), color='0.2', edgecolor='0.2') for v in legend_vals]
    labels = [f'{v:g}' for v in legend_vals]
    leg = ax.legend(handles, labels, title='-log10(FDR)', frameon=False,
                    loc='center left', bbox_to_anchor=(1.17, 0.52), fontsize=7, title_fontsize=8,
                    labelspacing=1.0, handletextpad=1.0)
    fig.subplots_adjust(left=0.31, right=0.80, top=0.84, bottom=0.20)
    # Save all formats.
    for folder, ext in [(PDF, 'pdf'), (SVG, 'svg'), (PNG, 'png')]:
        fig.savefig(folder / f'{out_stem}.{ext}', dpi=300, transparent=False)
    plt.close(fig)

plot_bubble(
    stage_a, stage_a_paths,
    'Stage_a / immature-like Sertoli cells',
    'c',
    'SupFig4C_Stage_a_immature_like_ST_pathway_bubble',
)
plot_bubble(
    stage_c, stage_c_paths,
    'Stage_c / mature-supportive Sertoli cells',
    'd',
    'SupFig4D_Stage_c_mature_supportive_ST_pathway_bubble',
)

# Combined two-panel preview/PDF for quick Illustrator placement.
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8.4, 8.6))
for ax, df, selected_paths, title, label in [
    (axes[0], stage_a, stage_a_paths, 'Stage_a / immature-like Sertoli cells', 'c'),
    (axes[1], stage_c, stage_c_paths, 'Stage_c / mature-supportive Sertoli cells', 'd')]:
    x_map = {g: i for i, g in enumerate(group_order)}
    y_order = selected_paths[::-1]
    y_map = {p: i for i, p in enumerate(y_order)}
    for _, r in df.iterrows():
        ax.scatter(x_map[str(r['group'])], y_map[str(r['pathway_readable'])],
                   s=size_from_fdr(r['neglog10_fdr']), c=[cmap(norm(r['NES']))],
                   edgecolor='0.25', linewidth=0.35, alpha=0.95)
    ax.set_xlim(-0.55, len(group_order)-0.45)
    ax.set_ylim(-0.6, len(y_order)-0.4)
    ax.set_xticks(range(len(group_order)))
    ax.set_xticklabels(group_order, rotation=45, ha='right', fontsize=8)
    ax.set_yticks(range(len(y_order)))
    ax.set_yticklabels([label_map.get(p, p) for p in y_order], fontsize=8)
    ax.grid(True, color='0.88', linewidth=0.7)
    ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_linewidth(0.8)
    ax.set_title(title, fontsize=10, weight='bold', pad=7)
    ax.text(-0.22, 1.10, label, transform=ax.transAxes, fontsize=14, weight='bold')
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02)
cbar.set_label('NES', fontsize=8)
cbar.ax.tick_params(labelsize=7)
legend_vals = [round(neg_min, 1), round((neg_min + neg_max)/2, 1), round(neg_max, 1)]
legend_vals = sorted(set(legend_vals))
handles = [plt.scatter([], [], s=size_from_fdr(v), color='0.2', edgecolor='0.2') for v in legend_vals]
axes[1].legend(handles, [f'{v:g}' for v in legend_vals], title='-log10(FDR)', frameon=False,
               loc='center left', bbox_to_anchor=(1.20, 0.88), fontsize=7, title_fontsize=8)
fig.subplots_adjust(left=0.30, right=0.80, top=0.94, bottom=0.08, hspace=0.55)
for folder, ext in [(PDF, 'pdf'), (SVG, 'svg'), (PNG, 'png')]:
    fig.savefig(folder / f'SupFig4CD_redrawn_pathway_bubbles_combined.{ext}', dpi=300)
plt.close(fig)

# Write README and copy code.
readme = '''SupFig4 redrawn Sertoli pathway bubble plots\n\nFiles included:\n- panels_pdf/: redrawn PDF panels for Illustrator\n- panels_svg/: editable SVG panels\n- panels_png/: high-resolution PNG previews\n- supfig4_source_data/: original GSEA tables plus selected plotting tables\n- supfig4_code/: plotting code\n\nPlot encoding:\n- x-axis: disease group compared with control\n- y-axis: selected, readable top pathways\n- color: normalized enrichment score (NES), blue negative and red positive\n- size: -log10(FDR)\n\nRecommended usage:\nUse the individual panel PDFs for Supplementary Figure 4 panels c and d.\nThe combined PDF is only a quick preview / layout helper.\n'''
(OUT / 'README.txt').write_text(readme, encoding='utf-8')
shutil.copy2(Path(__file__).resolve(), CODE / 'plot_supfig4_redrawn_bubbles.py')
# Zip package
zip_path = Path('/mnt/data/SupFig4_redrawn_ST_pathway_bubbles.zip')
if zip_path.exists():
    zip_path.unlink()
shutil.make_archive(str(zip_path.with_suffix('')), 'zip', OUT)
print('Wrote:', OUT)
print('Zip:', zip_path)
