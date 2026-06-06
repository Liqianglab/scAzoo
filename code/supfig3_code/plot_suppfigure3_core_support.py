#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'supfig3_source_data'
OUT = ROOT / 'supfig3_panels'
for sub in ['pdf', 'svg', 'png']:
    (OUT / sub).mkdir(parents=True, exist_ok=True)
(ROOT / 'preview').mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'axes.linewidth': 0.8,
})

GROUP_ORDER = ['Ctrl', 'OA', 'AZFc_Del', 'iNOA_B', 'iNOA_S', 'KS']
CELLTYPE_LABEL = {'Immune': 'Immune cells', 'Myoid': 'Myoid cells', 'EC': 'Endothelial cells'}


def save_panel(fig, name: str) -> None:
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT / ext / f'{name}.{ext}', bbox_inches='tight', dpi=300)
    plt.close(fig)


def plot_panel_a() -> None:
    df = pd.read_csv(SRC / 'SupFig3A_contextual_somatic_pathways_long.csv')
    selected = []
    for ct in ['Immune', 'Myoid', 'EC']:
        sub = df[df['panel_celltype'] == ct].copy()
        top_pathways = (
            sub.groupby('pathway_readable')['neglog10_fdr']
            .max()
            .sort_values(ascending=False)
            .head(7)
            .index.tolist()
        )
        selected.append(sub[sub['pathway_readable'].isin(top_pathways)])
    plot_df = pd.concat(selected, ignore_index=True)

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.5), sharex=True)
    norm = Normalize(vmin=-3, vmax=3)
    cmap = plt.cm.RdBu_r
    for ax, ct in zip(axes, ['Immune', 'Myoid', 'EC']):
        sub = plot_df[plot_df['panel_celltype'] == ct].copy()
        pathways = (
            sub.groupby('pathway_readable')['neglog10_fdr']
            .max()
            .sort_values(ascending=True)
            .index.tolist()
        )
        ymap = {p: i for i, p in enumerate(pathways)}
        xmap = {g: i for i, g in enumerate(GROUP_ORDER)}
        sizes = np.clip(sub['neglog10_fdr'].values, 0, 8) * 18 + 10
        colors = cmap(norm(sub['NES'].values))
        ax.scatter(sub['group'].map(xmap), sub['pathway_readable'].map(ymap),
                   s=sizes, c=colors, edgecolor='none')
        ax.set_yticks(range(len(pathways)))
        ax.set_yticklabels(pathways, fontsize=8)
        ax.set_xticks(range(len(GROUP_ORDER)))
        ax.set_xticklabels(GROUP_ORDER, rotation=45, ha='right', fontsize=8)
        ax.set_title(CELLTYPE_LABEL[ct], fontsize=11, fontweight='bold')
        ax.grid(axis='x', color='#dddddd', linewidth=0.5)
        ax.set_xlim(-0.5, len(GROUP_ORDER) - 0.5)
    axes[0].set_ylabel('Pathway')
    fig.suptitle('a  Contextual somatic pathway remodeling', x=0.01, ha='left', fontsize=14, fontweight='bold')
    sm = ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label('NES', fontsize=9)
    for s, lab in [(30, '1'), (80, '3'), (150, '6')]:
        axes[-1].scatter([], [], s=s, color='#777777', edgecolor='none', label=lab)
    axes[-1].legend(title='-log10(FDR)', loc='lower right', bbox_to_anchor=(1.38, -0.03), frameon=False, fontsize=8, title_fontsize=8)
    fig.tight_layout(rect=[0, 0, 0.94, 0.92])
    save_panel(fig, 'SupFig3A_contextual_somatic_pathways')


def plot_panel_b() -> None:
    mat = pd.read_csv(SRC / 'SupFig3B_core_LC_ST_germ_correlation_matrix.csv', index_col=0)
    txt = pd.read_csv(SRC / 'SupFig3B_core_LC_ST_germ_correlation_text_matrix.csv', index_col=0)
    labels = mat.index.tolist()
    fig, ax = plt.subplots(figsize=(7.5, 6.7))
    im = ax.imshow(mat.values, cmap='RdBu_r', vmin=-1, vmax=1)
    ax.set_title('b  Core LC-ST-germ donor-level correlations', loc='left', fontsize=14, fontweight='bold')
    ax.set_xticks(np.arange(len(labels)))
    ax.set_yticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xticks(np.arange(-0.5, len(labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(labels), 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=1)
    ax.tick_params(which='minor', bottom=False, left=False)
    for i in range(len(labels)):
        for j in range(len(labels)):
            s = str(txt.iloc[i, j])
            if s and s != 'nan':
                val = mat.iloc[i, j]
                color = 'white' if abs(val) > 0.65 else 'black'
                ax.text(j, i, s, ha='center', va='center', fontsize=7, color=color)
    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cbar.set_label('Spearman rho', fontsize=9)
    fig.tight_layout()
    save_panel(fig, 'SupFig3B_core_LC_ST_germ_correlations')


def plot_panel_c() -> None:
    mat = pd.read_csv(SRC / 'SupFig3C_core_LC_ST_germ_group_summary_matrix.csv', index_col=0)
    groups = [g for g in GROUP_ORDER if g in mat.columns]
    mat = mat[groups]
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    im = ax.imshow(mat.values, cmap='RdBu_r', vmin=-2, vmax=2, aspect='auto')
    ax.set_title('c  Group-level core somatic niche score summary', loc='left', fontsize=14, fontweight='bold')
    ax.set_yticks(np.arange(mat.shape[0]))
    ax.set_yticklabels(mat.index.tolist(), fontsize=9)
    ax.set_xticks(np.arange(len(groups)))
    ax.set_xticklabels(groups, rotation=45, ha='right', fontsize=9)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iloc[i, j]
            ax.text(j, i, f'{val:.2f}', ha='center', va='center', fontsize=7, color='black')
    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cbar.set_label('Mean row-wise z-score', fontsize=9)
    fig.tight_layout()
    save_panel(fig, 'SupFig3C_core_group_level_summary')


def make_preview() -> None:
    panel_paths = [
        OUT / 'png/SupFig3A_contextual_somatic_pathways.png',
        OUT / 'png/SupFig3B_core_LC_ST_germ_correlations.png',
        OUT / 'png/SupFig3C_core_group_level_summary.png',
    ]
    imgs = [Image.open(p).convert('RGB') for p in panel_paths]

    def resize(img, w):
        scale = w / img.width
        return img.resize((int(img.width * scale), int(img.height * scale)))

    img_a = resize(imgs[0], 1700)
    img_b = resize(imgs[1], 820)
    img_c = resize(imgs[2], 820)
    pad = 40
    title_h = 90
    width = 1700 + 2 * pad
    height = title_h + img_a.height + pad + max(img_b.height, img_c.height) + 2 * pad
    canvas = Image.new('RGB', (width, height), 'white')
    draw = ImageDraw.Draw(canvas)
    try:
        title_font = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf', 32)
        sub_font = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf', 18)
    except Exception:
        title_font = ImageFont.load_default()
        sub_font = ImageFont.load_default()
    draw.text((pad, 22), 'Supplementary Figure 3 - core support for revised Figure 3', font=title_font, fill='black')
    draw.text((pad, 60), 'Panels retained: contextual pathways + core LC-ST-germ correlations + core group-level niche summary', font=sub_font, fill='#555555')
    y = title_h
    canvas.paste(img_a, (pad, y))
    y += img_a.height + pad
    canvas.paste(img_b, (pad, y))
    canvas.paste(img_c, (pad + 860, y))
    canvas.save(ROOT / 'preview/Supplementary_Figure3_core_support_preview.png')


def main() -> None:
    plot_panel_a()
    plot_panel_b()
    plot_panel_c()
    make_preview()


if __name__ == '__main__':
    main()
