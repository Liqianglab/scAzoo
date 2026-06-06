from pathlib import Path
import os, zipfile, textwrap, math, shutil
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

ZIP = Path('/mnt/data/归档 5.zip')
WORK = Path('/mnt/data/supfig5_highlevel_work')
OUT = Path('/mnt/data/SupFig5_highlevel_GO_program_matrices')
if WORK.exists(): shutil.rmtree(WORK)
if OUT.exists(): shutil.rmtree(OUT)
WORK.mkdir(parents=True)
OUT.mkdir(parents=True)
for sub in ['panels_pdf','panels_svg','panels_png','source_data','code','preview']:
    (OUT/sub).mkdir(parents=True, exist_ok=True)
with zipfile.ZipFile(ZIP) as z:
    z.extractall(WORK)
BASE = WORK/'supfig5_source_data'

# Read full GO enrichment tables when available. They include GeneRatio and adjusted P values.
up_raw = pd.read_csv(BASE/'SupFig5B_GO_up_long.csv').dropna(subset=['group']).copy()
down_raw = pd.read_csv(BASE/'SupFig5B_GO_down_long.csv').dropna(subset=['group']).copy()

group_order = ['AZFc','iNOA_B','iNOA_S','KS']
group_display = {'AZFc':'AZFc_Del','iNOA_B':'iNOA_B','iNOA_S':'iNOA_S','KS':'KS'}

UP_MODULES = {
    'Collagen / ECM remodeling':[
        'collagen-containing extracellular matrix', 'extracellular matrix structural constituent',
        'extracellular matrix organization', 'extracellular structure organization',
        'external encapsulating structure organization'],
    'Adhesion / cell-substrate junction':['focal adhesion','cell-substrate junction','collagen binding','cell leading edge'],
    'Vesicle / secretory lumen':['vesicle lumen','cytoplasmic vesicle lumen','secretory granule lumen','endoplasmic reticulum lumen'],
    'Interstitial / tissue development':['reproductive system development','reproductive structure development','gland development','connective tissue development','muscle tissue development','renal system vasculature development','kidney vasculature development','metanephric nephron development','ossification','regeneration'],
    'Actin / contractile cytoskeleton':['actin binding','actin filament binding','contractile fiber','myofibril','structural constituent of muscle'],
    'Growth-factor / TGF-beta receptor axis':['transmembrane receptor protein serine/threonine kinase signaling pathway','growth factor binding','glycosaminoglycan binding','peptidase regulator activity'],
    'Wnt / transcriptional co-regulation':['wnt signaling','cell-cell signaling by wnt','DNA-binding transcription factor binding','transcription coregulator activity','RNA polymerase II-specific DNA-binding transcription factor binding'],
    'Oxidative-stress response':['response to oxidative stress'],
}

DOWN_MODULES = {
    'Steroid biosynthesis':['steroid biosynthetic process','sterol biosynthetic process'],
    'Steroid / cholesterol metabolism':['steroid metabolic process','sterol metabolic process','cholesterol metabolic process','secondary alcohol metabolic process','alcohol metabolic process'],
    'Mitochondrial compartment':['mitochondrial matrix','mitochondrial inner membrane','mitochondrial outer membrane','organelle outer membrane','outer membrane'],
    'OXPHOS / electron transport':['oxidative phosphorylation','electron transport chain','electron transfer activity'],
    'Respiration / ATP generation':['generation of precursor metabolites and energy','energy derivation by oxidation of organic compounds','aerobic respiration','cellular respiration','ATP metabolic process'],
    'Fatty-acid / lipid metabolism':['fatty acid metabolic process','lipid catabolic process'],
    'Vesicle / lysosomal lumen':['vacuolar lumen','secretory granule lumen','cytoplasmic vesicle lumen','vesicle lumen','primary lysosome','azurophil granule','ficolin-1-rich granule','ficolin-1-rich granule lumen'],
    'Proteostasis / chaperone':['protein folding chaperone','chaperone-mediated protein folding','endoplasmic reticulum protein-containing complex'],
    'Adhesion / junction':['focal adhesion','cell-substrate junction','cadherin binding'],
    'Oxidoreductase / antioxidant':['oxidoreductase activity','antioxidant activity'],
}

def ratio_float(s):
    if isinstance(s, str) and '/' in s:
        a,b = s.split('/')[:2]
        try:
            return float(a) / float(b)
        except Exception:
            return np.nan
    return np.nan

def select_representatives(df, modules, direction):
    out=[]
    for module, keywords in modules.items():
        for g in group_order:
            dfg = df[df['group'].eq(g)].copy()
            desc = dfg['Description'].astype(str).str.lower()
            mask = desc.apply(lambda x: any(k.lower() in x for k in keywords))
            sub = dfg[mask].copy()
            if sub.empty:
                out.append({
                    'direction':direction,'group':g,'group_display':group_display[g],
                    'module':module,'representative_GO_term':'', 'GO_ID':'',
                    'Count':0,'GeneRatio':'','gene_ratio':np.nan,
                    'p.adjust':np.nan,'minus_log10_p_adjust':np.nan,
                    'geneID':''
                })
            else:
                # Use the most significant representative term within each program for each disease group.
                r = sub.sort_values(['minus_log10_p_adjust','Count'], ascending=[False, False]).iloc[0]
                out.append({
                    'direction':direction,'group':g,'group_display':group_display[g],
                    'module':module,'representative_GO_term':r.get('Description',''), 'GO_ID':r.get('ID',''),
                    'Count':r.get('Count',np.nan),'GeneRatio':r.get('GeneRatio',''), 'gene_ratio':ratio_float(r.get('GeneRatio','')),
                    'p.adjust':r.get('p.adjust',np.nan),'minus_log10_p_adjust':r.get('minus_log10_p_adjust',np.nan),
                    'geneID':r.get('geneID','')
                })
    return pd.DataFrame(out)

up = select_representatives(up_raw, UP_MODULES, 'Up')
down = select_representatives(down_raw, DOWN_MODULES, 'Down')
up.to_csv(OUT/'source_data'/'SupFig5B_highlevel_up_GO_program_matrix_source.csv', index=False)
down.to_csv(OUT/'source_data'/'SupFig5C_highlevel_down_GO_program_matrix_source.csv', index=False)
# Also copy raw input tables for traceability.
for fn in ['SupFig5B_GO_up_long.csv','SupFig5B_GO_down_long.csv','SupFig5B_GO_up_top15_by_group.csv','SupFig5B_GO_down_top15_by_group.csv']:
    if (BASE/fn).exists():
        shutil.copy(BASE/fn, OUT/'source_data'/fn)

# Plotting utilities
up_cmap = LinearSegmentedColormap.from_list('up_red', ['#fff5eb','#fdd0a2','#fb6a4a','#cb181d','#67000d'])
down_cmap = LinearSegmentedColormap.from_list('down_blue', ['#f7fbff','#c6dbef','#6baed6','#2171b5','#08306b'])

def wrap_label(s, width=28):
    return '\n'.join(textwrap.wrap(s, width=width, break_long_words=False))

def draw_dot_matrix(df, modules, title, letter, cmap, outfile_stem, vmax=None):
    module_order = list(modules.keys())
    df = df.copy()
    df['module'] = pd.Categorical(df['module'], categories=module_order, ordered=True)
    df['group'] = pd.Categorical(df['group'], categories=group_order, ordered=True)
    valid = df.dropna(subset=['minus_log10_p_adjust'])
    if vmax is None:
        vmax = max(3, float(np.nanmax(valid['minus_log10_p_adjust'])) if not valid.empty else 3)
    norm = Normalize(vmin=0, vmax=vmax)
    # Determine size from gene ratio; fallback count if missing. Size is area.
    ratios = df['gene_ratio'].dropna()
    min_ratio = max(0.03, float(ratios.min()) if len(ratios) else 0.05)
    max_ratio = max(0.20, float(ratios.max()) if len(ratios) else 0.20)
    def size_map(r):
        if pd.isna(r) or r<=0:
            return 0
        return 80 + 650*(r-min_ratio)/(max_ratio-min_ratio+1e-9)
    df['size'] = df['gene_ratio'].apply(size_map)
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    # background grid
    for x in range(len(group_order)):
        ax.axvline(x, color='#e6e6e6', lw=0.8, zorder=0)
    for y in range(len(module_order)):
        ax.axhline(y, color='#e6e6e6', lw=0.8, zorder=0)
    # scatter
    for _, r in df.dropna(subset=['minus_log10_p_adjust']).iterrows():
        x = group_order.index(str(r['group']))
        y = module_order.index(str(r['module']))
        ax.scatter(x, y, s=r['size'], c=[cmap(norm(r['minus_log10_p_adjust']))], edgecolor='black', linewidth=0.35, zorder=3)
    ax.set_xlim(-0.5, len(group_order)-0.5)
    ax.set_ylim(-0.5, len(module_order)-0.5)
    ax.invert_yaxis()
    ax.set_xticks(range(len(group_order)))
    ax.set_xticklabels([group_display[g] for g in group_order], rotation=35, ha='right', fontsize=10)
    ax.set_yticks(range(len(module_order)))
    ax.set_yticklabels([wrap_label(m, 25) for m in module_order], fontsize=9)
    ax.tick_params(axis='both', length=0)
    for spine in ['top','right']:
        ax.spines[spine].set_visible(False)
    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
    ax.text(-0.13, 1.08, letter, transform=ax.transAxes, fontsize=16, fontweight='bold', va='bottom')
    # colorbar
    sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.025)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(r'$-\log_{10}$(adj. P)', fontsize=9)
    # Size legend with representative gene ratios.
    # choose nice sizes based on actual range
    ratio_vals = [0.05,0.10,0.15,0.20]
    ratio_vals = [r for r in ratio_vals if r >= min_ratio-0.01 and r <= max_ratio+0.03]
    if len(ratio_vals)<3:
        ratio_vals = np.linspace(min_ratio, max_ratio, 3)
    handles=[]
    labels=[]
    for rv in ratio_vals:
        handles.append(ax.scatter([],[],s=size_map(rv), c='lightgrey', edgecolor='black', linewidth=0.35))
        labels.append(f'{rv:.2f}')
    leg = ax.legend(handles, labels, title='Gene ratio', loc='upper left', bbox_to_anchor=(1.14, 0.98), frameon=False, fontsize=8, title_fontsize=9, borderaxespad=0)
    # save
    fig.savefig(OUT/'panels_pdf'/f'{outfile_stem}.pdf', bbox_inches='tight')
    fig.savefig(OUT/'panels_svg'/f'{outfile_stem}.svg', bbox_inches='tight')
    fig.savefig(OUT/'panels_png'/f'{outfile_stem}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    return vmax

vmax_up = draw_dot_matrix(up, UP_MODULES, 'Upregulated GO programs in disease Leydig cells', 'b', up_cmap, 'SupFig5B_highlevel_up_GO_program_matrix')
vmax_down = draw_dot_matrix(down, DOWN_MODULES, 'Downregulated GO programs in disease Leydig cells', 'c', down_cmap, 'SupFig5C_highlevel_down_GO_program_matrix')

# Combined side-by-side layout for AI preview.
fig, axes = plt.subplots(1,2, figsize=(13.2,5.0), gridspec_kw={'wspace':0.85})
# re-draw on provided axes simplified with shared functions inline

def draw_on_ax(ax, df, modules, title, letter, cmap, vmax):
    module_order=list(modules.keys())
    norm=Normalize(vmin=0, vmax=vmax)
    ratios=df['gene_ratio'].dropna()
    min_ratio=max(0.03,float(ratios.min()) if len(ratios) else 0.05)
    max_ratio=max(0.20,float(ratios.max()) if len(ratios) else 0.20)
    def size_map(r):
        if pd.isna(r) or r<=0: return 0
        return 60+520*(r-min_ratio)/(max_ratio-min_ratio+1e-9)
    for x in range(len(group_order)): ax.axvline(x, color='#e6e6e6', lw=0.8, zorder=0)
    for y in range(len(module_order)): ax.axhline(y, color='#e6e6e6', lw=0.8, zorder=0)
    for _,r in df.dropna(subset=['minus_log10_p_adjust']).iterrows():
        x=group_order.index(str(r['group'])); y=module_order.index(str(r['module']))
        ax.scatter(x,y,s=size_map(r['gene_ratio']), c=[cmap(norm(r['minus_log10_p_adjust']))], edgecolor='black', lw=0.3, zorder=3)
    ax.set_xlim(-0.5,len(group_order)-0.5); ax.set_ylim(-0.5,len(module_order)-0.5); ax.invert_yaxis()
    ax.set_xticks(range(len(group_order))); ax.set_xticklabels([group_display[g] for g in group_order], rotation=35, ha='right', fontsize=9)
    ax.set_yticks(range(len(module_order))); ax.set_yticklabels([wrap_label(m, 23) for m in module_order], fontsize=8)
    ax.tick_params(length=0)
    for spine in ['top','right']: ax.spines[spine].set_visible(False)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.text(-0.17, 1.07, letter, transform=ax.transAxes, fontsize=16, fontweight='bold')

draw_on_ax(axes[0], up, UP_MODULES, 'Upregulated GO programs', 'b', up_cmap, vmax_up)
draw_on_ax(axes[1], down, DOWN_MODULES, 'Downregulated GO programs', 'c', down_cmap, vmax_down)
# add independent colorbars
for ax, cmap, vmax in zip(axes,[up_cmap,down_cmap],[vmax_up,vmax_down]):
    sm=ScalarMappable(norm=Normalize(vmin=0, vmax=vmax), cmap=cmap); sm.set_array([])
    cb=fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.02)
    cb.ax.tick_params(labelsize=7); cb.set_label(r'$-\log_{10}$(adj. P)', fontsize=8)
fig.savefig(OUT/'panels_pdf'/'SupFig5BC_highlevel_GO_program_matrices_combined.pdf', bbox_inches='tight')
fig.savefig(OUT/'panels_svg'/'SupFig5BC_highlevel_GO_program_matrices_combined.svg', bbox_inches='tight')
fig.savefig(OUT/'panels_png'/'SupFig5BC_highlevel_GO_program_matrices_combined.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# Write reproducible plotting script
script_path = OUT/'code'/'draw_supfig5_highlevel_GO_program_matrices.py'
with open(script_path, 'w') as f:
    f.write(Path(__file__).read_text() if '__file__' in globals() and Path(__file__).exists() else '# Reproducible script was generated by ChatGPT. See source_data and panel files.\n')
# Write README
with open(OUT/'README.txt','w') as f:
    f.write('SupFig5 high-level GO program matrices\n')
    f.write('Input: SupFig5B_GO_up_long.csv and SupFig5B_GO_down_long.csv from the uploaded SupFig5 data package.\n')
    f.write('Method: redundant GO terms were grouped into curated Leydig-cell functional programs. For each disease group and program, the most significant matching GO term was selected as the representative term. Dot color encodes -log10(adjusted P); dot size encodes GeneRatio.\n')
    f.write('Outputs: PDF, SVG and PNG files in panels_* folders; source_data contains the selected representative terms.\n')
# zip
zip_out=Path('/mnt/data/SupFig5_highlevel_GO_program_matrices.zip')
if zip_out.exists(): zip_out.unlink()
with zipfile.ZipFile(zip_out,'w',zipfile.ZIP_DEFLATED) as z:
    for root, dirs, files in os.walk(OUT):
        for file in files:
            full=Path(root)/file
            z.write(full, full.relative_to(OUT.parent))
print(zip_out)
