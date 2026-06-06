from __future__ import annotations
from pathlib import Path
import shutil, zipfile, textwrap
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D
from PIL import Image, ImageDraw, ImageFont

BASE=Path('/mnt/data')
OLD_ZIP=BASE/'归档.zip'
NICHE_ZIP=BASE/'归档 3.zip'
BCD_ZIP=BASE/'Figure3_BCD_options_package.zip'
OUT=BASE/'Figure3_final_realdata_AI_panels'
TMP=OUT/'_tmp'
DATA=OUT/'figure3_source_data'
CODE=OUT/'figure3_code'
PDF=OUT/'panels_pdf'
SVG=OUT/'panels_svg'
PNG=OUT/'panels_png'
PREV=OUT/'preview'
if OUT.exists(): shutil.rmtree(OUT)
for d in [TMP,DATA,CODE,PDF,SVG,PNG,PREV]: d.mkdir(parents=True, exist_ok=True)
with zipfile.ZipFile(OLD_ZIP) as z: z.extractall(TMP/'old')
with zipfile.ZipFile(NICHE_ZIP) as z: z.extractall(TMP/'niche')
with zipfile.ZipFile(BCD_ZIP) as z: z.extract('B_panel_priority_metrics.csv', TMP/'bmetrics')
OLD_SRC=TMP/'old'/'figure3_source_data'
NICHE_SRC=TMP/'niche'/'supfig8_source_data'
B_METRICS=TMP/'bmetrics'/'B_panel_priority_metrics.csv'

mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42
mpl.rcParams['svg.fonttype']='none'
mpl.rcParams['font.family']='DejaVu Sans'
mpl.rcParams['axes.linewidth']=0.8

cell_cols={'Endothelial cells':'#1f77b4','Leydig cells':'#ff7f0e','Immune cells':'#2ca02c','Myoid cells':'#d62728','Sertoli cells':'#8c564b'}
group_cols={'Ctrl':'#1f77b4','OA':'#d62728','AZFc_Del':'#9467bd','iNOA_B':'#ff7f0e','iNOA_S':'#2ca02c','KS':'#8c564b'}
group_order=['Ctrl','OA','AZFc_Del','iNOA_B','iNOA_S','KS']
disease_order=['OA','AZFc_Del','iNOA_B','iNOA_S','KS']

# Load raw data
umap=pd.read_csv(OLD_SRC/'Fig3A_somatic_umap_cells.csv')
pathways=pd.read_csv(OLD_SRC/'Fig3C_bubble_pathways_long.csv')
priority=pd.read_csv(B_METRICS)
scores_z=pd.read_csv(NICHE_SRC/'SupFig8A_somatic_scores_heatmap_z_matrix.csv')
scores_long=pd.read_csv(NICHE_SRC/'SupFig8A_somatic_scores_long.csv')
sample_order=pd.read_csv(NICHE_SRC/'SupFig8A_sample_order_by_PC1v2.csv')
group_bar=pd.read_csv(NICHE_SRC/'SupFig8A_group_bar_annotation.csv')
corr_long=pd.read_csv(NICHE_SRC/'SupFig8E_score_correlation_long.csv')
corr_mat=pd.read_csv(NICHE_SRC/'SupFig8E_score_correlation_spearman_matrix.csv', index_col=0)
pca=pd.read_csv(NICHE_SRC/'SupFig8B_pca_points.csv')
germ_mat=pd.read_csv(NICHE_SRC/'SupFig8C_niche_vs_germ_maturity.csv')
germ_latest=pd.read_csv(NICHE_SRC/'SupFig8D_niche_vs_germ_latest_stage.csv')
sum_mat=pd.read_csv(NICHE_SRC/'SupFig8C_spearman_summary.csv')
sum_latest=pd.read_csv(NICHE_SRC/'SupFig8D_spearman_summary.csv')

# Copy source files / selected data
umap.to_csv(DATA/'Fig3A_somatic_umap_cells.csv', index=False)
priority.to_csv(DATA/'Fig3B_integrated_priority_metrics.csv', index=False)
# Original remodeling table for traceability
pd.read_csv(OLD_SRC/'Fig3B_radar_matrix_long.csv').to_csv(DATA/'Fig3B_original_remodeling_signal.csv', index=False)

lc_map={'Cholesterol\nhomeostasis':'Cholesterol Homeostasis','Fatty acid\nmetabolism':'Fatty Acid Metabolism','OXPHOS':'Oxidative Phosphorylation','TNF-alpha/\nNF-kB':'Tnfa Signaling Via Nfkb','ECM / adhesion':'Epithelial Mesenchymal Transition','Hypoxia':'Hypoxia','Apoptosis':'Apoptosis'}
st_map={'ECM / adhesion':'Epithelial Mesenchymal Transition','TNF-alpha/\nNF-kB':'Tnfa Signaling Via Nfkb','IFN-gamma\nresponse':'Interferon Gamma Response','Apoptosis':'Apoptosis','p53 stress':'P53 Pathway','Glycolysis':'Glycolysis','OXPHOS':'Oxidative Phosphorylation'}
sel=set(lc_map.values())|set(st_map.values())
panelC_src=pathways[(pathways.panel_celltype.isin(['LC','ST'])) & (pathways.pathway_readable.isin(sel))].copy()
panelC_src.to_csv(DATA/'Fig3C_LC_ST_pathway_enrichment_selected.csv', index=False)

ordered_samples=sample_order.sort_values('sample_rank_pc1v2')['Sample_ID'].tolist()
core_map={'ST Maturity':'Sertoli maturity','BTB Integrity':'BTB integrity','ST Stress IEG':'Sertoli stress','ST SASP':'Sertoli SASP','LC Mature':'Leydig maturity','LC Immature':'Leydig immaturity'}
score_core=scores_z.set_index('score_display').loc[list(core_map.keys()), ordered_samples]
score_core.index=[core_map[x] for x in score_core.index]
score_core.reset_index(names='score_display').to_csv(DATA/'Fig3D_core_somatic_niche_scores_z_matrix.csv', index=False)
sample_order.to_csv(DATA/'Fig3D_sample_order_by_PC1v2.csv', index=False)
group_bar.to_csv(DATA/'Fig3D_group_bar_annotation.csv', index=False)
# Optional extended modules with HH/lactate for supplement or later decisions
scores_z.to_csv(DATA/'Fig3D_optional_extended_all_niche_scores_z_matrix.csv', index=False)

# Correlation edges for panel E
def fetch_corr(a,b):
    r=corr_long[(corr_long.row_score==a)&(corr_long.col_score==b)]
    if r.empty: r=corr_long[(corr_long.row_score==b)&(corr_long.col_score==a)]
    if r.empty: return None
    return r.iloc[0]
edge_pairs=[('LC Mature','ST Maturity'),('LC Mature','BTB Integrity'),('LC Immature','ST Maturity'),('LC Immature','BTB Integrity'),('LC Immature','Germ Maturity Index'),('LC Immature','Germ Latest Stage'),('ST Maturity','Germ Maturity Index'),('ST Maturity','Germ Latest Stage'),('BTB Integrity','Germ Maturity Index'),('BTB Integrity','Germ Latest Stage')]
edges=[]
for s,t in edge_pairs:
    r=fetch_corr(s,t)
    if r is not None:
        edges.append({'source':s,'target':t,'spearman_rho':float(r['spearman_rho']),'p_value':float(r['p_value']),'fdr_bh_q':float(r['fdr_bh_q']) if 'fdr_bh_q' in r else np.nan,'stars':r.get('stars','')})
edges=pd.DataFrame(edges)
edges.to_csv(DATA/'Fig3E_core_LC_ST_germ_association_edges.csv', index=False)
corr_mat.to_csv(DATA/'Fig3E_full_score_correlation_matrix.csv')
corr_long.to_csv(DATA/'Fig3E_full_score_correlation_long.csv', index=False)
pca.to_csv(DATA/'Fig3F_PCA_niche_axis_points.csv', index=False)
germ_mat.to_csv(DATA/'Fig3G_niche_axis_vs_germ_maturity.csv', index=False)
germ_latest.to_csv(DATA/'Fig3G_niche_axis_vs_latest_stage.csv', index=False)
sum_mat.to_csv(DATA/'Fig3G_germ_maturity_spearman_summary.csv', index=False)
sum_latest.to_csv(DATA/'Fig3G_latest_stage_spearman_summary.csv', index=False)

# Helpers
def savefig(fig,name):
    fig.savefig(PDF/f'{name}.pdf', bbox_inches='tight')
    fig.savefig(SVG/f'{name}.svg', bbox_inches='tight')
    fig.savefig(PNG/f'{name}.png', dpi=500, bbox_inches='tight')
    plt.close(fig)

def clean(ax):
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

def label(ax, letter, x=-0.12, y=1.12):
    ax.text(x,y,letter,transform=ax.transAxes,ha='left',va='top',fontsize=16,fontweight='bold')

# Panel A: Somatic UMAP
fig, ax=plt.subplots(figsize=(4.2,3.4))
for ct in ['Endothelial cells','Leydig cells','Immune cells','Myoid cells','Sertoli cells']:
    sub=umap[umap.somatic_celltype==ct]
    ax.scatter(sub.UMAP1, sub.UMAP2, s=0.45, c=cell_cols[ct], alpha=0.9, linewidths=0, label=ct)
ax.set_title('Somatic UMAP', fontsize=11, pad=4)
ax.set_xlabel('UMAP1', fontsize=9); ax.set_ylabel('UMAP2', fontsize=9)
ax.set_xticks([]); ax.set_yticks([])
# legend to right, away from cells
ax.legend(frameon=False, fontsize=7, markerscale=6, handlelength=0.8, loc='lower right', bbox_to_anchor=(1.03,0.02))
clean(ax); label(ax,'a')
savefig(fig,'Fig3A_somatic_UMAP')

# Panel B: prioritization
fig, ax=plt.subplots(figsize=(7.4,3.15))
role={'Leydig':'endocrine / interstitial support','Sertoli':'tubular germline support','Immune':'inflammatory niche context','Myoid':'peritubular structural support','Endothelial':'vascular niche support'}
col={'Leydig':cell_cols['Leydig cells'],'Sertoli':cell_cols['Sertoli cells'],'Immune':'#BDBDBD','Myoid':'#BDBDBD','Endothelial':'#BDBDBD'}
df=priority.copy().sort_values('integrated_priority', ascending=True)
y=np.arange(len(df))
ax.barh(y, df.integrated_priority, color=[col[d] for d in df.display], height=0.62)
ax.scatter(df.integrated_priority, y, s=95, color=[col[d] for d in df.display], edgecolor='black', zorder=5)
ranks=df.integrated_priority.rank(ascending=False, method='first').astype(int)
for yi,(_,r),rk in zip(y, df.iterrows(), ranks):
    ax.text(-0.038, yi, str(rk), va='center', ha='center', fontsize=8, fontweight='bold', bbox=dict(boxstyle='circle,pad=0.18', fc='white', ec='black', lw=0.8))
    weight='bold' if r.display in ['Leydig','Sertoli'] else 'normal'
    ax.text(r.integrated_priority+0.02, yi, role[r.display], va='center', fontsize=8.4, fontweight=weight, color='black' if weight=='bold' else '#444')
ax.set_yticks(y); ax.set_yticklabels([d+' cells' for d in df.display], fontsize=9)
ax.set_xlabel('Integrated niche-driver priority score', fontsize=9)
ax.set_xlim(-0.07,0.86)
ax.set_title('Somatic niche prioritization', fontsize=11, fontweight='bold', pad=4)
ax.grid(axis='x', alpha=0.25); clean(ax); label(ax,'b')
savefig(fig,'Fig3B_somatic_niche_prioritization')

# Panel C: LC/ST pathway bubbles
def pathway_ax(ax, ct, mapping, title, title_color):
    rows=list(mapping.keys())[::-1]
    xmap={g:i for i,g in enumerate(disease_order)}; ymap={r:i for i,r in enumerate(rows)}
    for disp, raw in mapping.items():
        sub=panelC_src[(panelC_src.panel_celltype==ct)&(panelC_src.pathway_readable==raw)]
        for _,r in sub.iterrows():
            if r.group not in xmap: continue
            size=16+min(float(r.neglog10_fdr),8)*18
            ax.scatter(xmap[r.group], ymap[disp], s=size, c=[r.NES], cmap='RdBu_r', vmin=-3, vmax=3, edgecolor='#555', linewidth=0.35)
    ax.set_xticks(list(xmap.values())); ax.set_xticklabels(disease_order, rotation=45, ha='right', fontsize=8)
    ax.set_yticks(np.arange(len(rows))); ax.set_yticklabels(rows, fontsize=8)
    ax.set_title(title, fontsize=11, fontweight='bold', color=title_color, pad=3)
    ax.grid(True, alpha=0.25); ax.set_axisbelow(True); clean(ax)

fig=plt.figure(figsize=(11.8,3.45))
gs=GridSpec(1,3,width_ratios=[1,1,0.06], wspace=0.33, figure=fig)
ax1=fig.add_subplot(gs[0,0]); ax2=fig.add_subplot(gs[0,1]); cax=fig.add_subplot(gs[0,2])
pathway_ax(ax1,'LC',lc_map,'LC pathways',cell_cols['Leydig cells'])
pathway_ax(ax2,'ST',st_map,'ST pathways',cell_cols['Sertoli cells'])
sm=mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=-3,vmax=3), cmap='RdBu_r')
cb=fig.colorbar(sm,cax=cax); cb.set_label('NES', fontsize=9)
for s,lab in zip([16+2*18,16+4*18,16+8*18],['2','4','8']):
    ax2.scatter([],[],s=s,color='white',edgecolor='black',label=lab)
ax2.legend(title='-log10(FDR)', frameon=False, fontsize=7, title_fontsize=7, loc='lower right', bbox_to_anchor=(1.05,-0.05))
ax1.text(-0.18,1.13,'c',transform=ax1.transAxes,fontsize=16,fontweight='bold')
fig.suptitle('LC and ST disease-associated pathways', y=1.02, fontsize=12, fontweight='bold')
savefig(fig,'Fig3C_LC_ST_pathway_enrichment')

# Panel D: long heatmap
fig=plt.figure(figsize=(12.4,3.6))
gs=GridSpec(2,2,height_ratios=[0.24,1],width_ratios=[1,0.025],hspace=0.08,wspace=0.03,figure=fig)
axg=fig.add_subplot(gs[0,0]); ax=fig.add_subplot(gs[1,0]); cax=fig.add_subplot(gs[1,1])
sample_groups=group_bar.set_index('Sample_ID').loc[ordered_samples,'Group']
rgb=np.array([mpl.colors.to_rgb(group_cols[g]) for g in sample_groups])[None,:,:]
axg.imshow(rgb,aspect='auto'); axg.set_xticks([]); axg.set_yticks([])
for sp in axg.spines.values(): sp.set_visible(False)
im=ax.imshow(score_core.values,cmap='RdBu_r',vmin=-2.3,vmax=2.3,aspect='auto')
ax.set_yticks(np.arange(score_core.shape[0])); ax.set_yticklabels(score_core.index.tolist(), fontsize=9)
ax.set_xticks(np.arange(len(ordered_samples))); ax.set_xticklabels(ordered_samples, rotation=45, ha='right', fontsize=7)
ax.set_xlabel('Donors ordered by niche-collapse axis', fontsize=9)
ax.hlines([3.5], -0.5, len(ordered_samples)-0.5, color='white', lw=1.5)
cb=fig.colorbar(im,cax=cax); cb.set_label('Z-score', fontsize=9)
axg.set_title('Integrated somatic niche scores', fontsize=12, fontweight='bold', pad=14)
axg.text(0,1.58,'preserved',color='#2b6cb0',fontsize=9,transform=axg.transAxes,ha='left')
axg.annotate('niche-collapse axis', xy=(0.82,1.58), xytext=(0.32,1.58), xycoords='axes fraction', textcoords='axes fraction', arrowprops=dict(arrowstyle='->', lw=1), fontsize=9, ha='center')
axg.text(1,1.58,'collapsed',color='#c53030',fontsize=9,transform=axg.transAxes,ha='right')
handles=[Line2D([0],[0],marker='s',linestyle='',color='none',markerfacecolor=group_cols[g],markersize=7,label=g) for g in group_order]
axg.legend(handles=handles,ncol=6,frameon=False,fontsize=7,loc='upper right',bbox_to_anchor=(1,1.45))
axg.text(-0.08,1.60,'d',transform=axg.transAxes,fontsize=16,fontweight='bold')
savefig(fig,'Fig3D_integrated_somatic_niche_scores')

# Panel E: association map
def draw_box(ax, xy, txt, color):
    x,y=xy; w,h=0.27,0.085
    patch=FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0.015,rounding_size=0.02',fc='white',ec=color,lw=1.2)
    ax.add_patch(patch); ax.text(x,y,txt,fontsize=7.7,ha='center',va='center')
def draw_edge(ax,p1,p2,rho):
    c='#b2182b' if rho>=0 else '#2166ac'; lw=1+3*min(abs(rho),1); ls='-' if abs(rho)>=0.6 else '--'
    arr=FancyArrowPatch(p1,p2,arrowstyle='-|>',mutation_scale=8,linewidth=lw,linestyle=ls,color=c,alpha=0.85,connectionstyle='arc3,rad=0.05')
    ax.add_patch(arr)
node_pos={'LC Mature':(0.16,0.70),'LC Immature':(0.16,0.42),'ST Maturity':(0.50,0.76),'BTB Integrity':(0.50,0.58),'Germ Maturity Index':(0.84,0.70),'Germ Latest Stage':(0.84,0.48)}
node_lab={'LC Mature':'LC mature','LC Immature':'LC immature','ST Maturity':'ST maturity','BTB Integrity':'BTB integrity','Germ Maturity Index':'Germ maturity','Germ Latest Stage':'Latest stage'}
fig,ax=plt.subplots(figsize=(4.8,3.6)); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
ax.set_title('Core LC-ST-germ association map', fontsize=11, fontweight='bold')
ax.text(0.16,0.92,'Leydig state',color=cell_cols['Leydig cells'],fontsize=8,fontweight='bold',ha='center')
ax.text(0.50,0.92,'Sertoli support',color=cell_cols['Sertoli cells'],fontsize=8,fontweight='bold',ha='center')
ax.text(0.84,0.92,'Germline output',color='#2f855a',fontsize=8,fontweight='bold',ha='center')
for _,r in edges.iterrows(): draw_edge(ax,node_pos[r.source],node_pos[r.target],r.spearman_rho)
for k,pos in node_pos.items():
    color=cell_cols['Leydig cells'] if k.startswith('LC') else cell_cols['Sertoli cells'] if k in ['ST Maturity','BTB Integrity'] else '#2f855a'
    draw_box(ax,pos,node_lab[k],color)
ax.plot([0.08,0.20],[0.12,0.12],color='#b2182b',lw=2.5); ax.text(0.22,0.12,'positive',fontsize=7,va='center')
ax.plot([0.08,0.20],[0.06,0.06],color='#2166ac',lw=2.5); ax.text(0.22,0.06,'negative',fontsize=7,va='center')
ax.text(0.56,0.08,'line width = |Spearman rho|',fontsize=7,color='#555')
label(ax,'e',x=-0.02,y=1.10)
savefig(fig,'Fig3E_core_LC_ST_germ_association_map')

# Panel F PCA
fig,ax=plt.subplots(figsize=(4.0,3.6))
for g in group_order:
    sub=pca[pca.Group==g]; ax.scatter(sub.PC1_v2,sub.PC2_v2,s=40,color=group_cols[g],edgecolor='none',label=g)
ax.axhline(0,color='#888',lw=0.7); ax.axvline(0,color='#888',lw=0.7)
ax.annotate('niche collapse',xy=(2.6,-1.15),xytext=(-0.5,-2.0),arrowprops=dict(arrowstyle='->',lw=1),fontsize=8)
ax.set_xlabel('PC1 (40.8%)',fontsize=9); ax.set_ylabel('PC2 (31.4%)',fontsize=9); ax.set_title('Donor-level niche axis',fontsize=11,fontweight='bold')
clean(ax); label(ax,'f',x=-0.17)
savefig(fig,'Fig3F_donor_level_niche_axis_PCA')

# Panel G
def scatter_fit(ax,df,ycol,summary,title,ylabel):
    for g in group_order:
        sub=df[df.Group==g]; ax.scatter(sub.x_niche_collapse_axis,sub[ycol],s=38,color=group_cols[g],edgecolor='none',label=g)
    x=df.x_niche_collapse_axis.to_numpy(float); y=df[ycol].to_numpy(float); m,b=np.polyfit(x,y,1); xx=np.linspace(x.min(),x.max(),100)
    ax.plot(xx,m*xx+b,color='black',lw=0.8)
    rho=float(summary.iloc[0].spearman_rho); p=float(summary.iloc[0].p_value)
    ax.text(0.05,0.94,f'Spearman rho = {rho:.2f}\nP = {p:.4g}',transform=ax.transAxes,fontsize=8,va='top')
    ax.set_title(title,fontsize=9); ax.set_xlabel('Niche-collapse axis (PC1v2)',fontsize=8); ax.set_ylabel(ylabel,fontsize=8); clean(ax)
fig=plt.figure(figsize=(6.5,3.6)); gs=GridSpec(1,2,wspace=0.34,figure=fig); ax1=fig.add_subplot(gs[0,0]); ax2=fig.add_subplot(gs[0,1])
scatter_fit(ax1,germ_mat,'y_germ_maturity_index',sum_mat,'Germ maturity','Germ maturity index')
scatter_fit(ax2,germ_latest,'y_germline_latest_stage_0to7',sum_latest,'Latest germ stage','Germline latest stage')
ax1.text(-0.22,1.14,'g',transform=ax1.transAxes,fontsize=16,fontweight='bold')
fig.suptitle('Niche collapse predicts germline outcome',fontsize=11,fontweight='bold',y=1.03)
savefig(fig,'Fig3G_niche_axis_predicts_germline_outcome')

# Composite preview via PIL, not intended as editable source
panel_files=['Fig3A_somatic_UMAP.png','Fig3B_somatic_niche_prioritization.png','Fig3C_LC_ST_pathway_enrichment.png','Fig3D_integrated_somatic_niche_scores.png','Fig3E_core_LC_ST_germ_association_map.png','Fig3F_donor_level_niche_axis_PCA.png','Fig3G_niche_axis_predicts_germline_outcome.png']
imgs={f:Image.open(PNG/f).convert('RGB') for f in panel_files}
# resize utility
def rw(img,w):
    scale=w/img.width; return img.resize((int(w),int(img.height*scale)))
W=2400; margin=60; gap=36
row1h=560; row2h=520; row3h=520; row4h=520
# target widths
A=rw(imgs['Fig3A_somatic_UMAP.png'],720); B=rw(imgs['Fig3B_somatic_niche_prioritization.png'],1500)
C=rw(imgs['Fig3C_LC_ST_pathway_enrichment.png'],W-2*margin)
D=rw(imgs['Fig3D_integrated_somatic_niche_scores.png'],W-2*margin)
E=rw(imgs['Fig3E_core_LC_ST_germ_association_map.png'],700); F=rw(imgs['Fig3F_donor_level_niche_axis_PCA.png'],560); G=rw(imgs['Fig3G_niche_axis_predicts_germline_outcome.png'],920)
H=90+A.height+gap+C.height+gap+D.height+gap+max(E.height,F.height,G.height)+margin
canvas=Image.new('RGB',(W,H),'white'); draw=ImageDraw.Draw(canvas)
try: font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf',36)
except: font=ImageFont.load_default()
draw.text((margin,20),'Figure 3',font=font,fill='black')
y=80; canvas.paste(A,(margin,y)); canvas.paste(B,(W-margin-B.width,y))
y+=max(A.height,B.height)+gap; canvas.paste(C,(margin,y))
y+=C.height+gap; canvas.paste(D,(margin,y))
y+=D.height+gap; x=margin; canvas.paste(E,(x,y)); x+=E.width+gap; canvas.paste(F,(x,y)); x+=F.width+gap; canvas.paste(G,(x,y))
canvas.save(PREV/'Figure3_realdata_layout_preview.png')
canvas.save(PNG/'Figure3_composite_layout_preview.png')
canvas.save(PDF/'Figure3_composite_layout_preview.pdf', 'PDF', resolution=300)

# Write code and README
shutil.copyfile(Path(__file__), CODE/'plot_figure3_realdata_AI_panels.py')
(DATA/'README.txt').write_text('Source data files for revised Figure 3 panels generated from uploaded real data. Panel D deliberately uses core ST/LC maturity, BTB and stress/SASP modules only; HH and metabolic candidate modules are saved as optional extended source data for later figures/supplement.\n',encoding='utf-8')
(OUT/'README.txt').write_text('Figure3_final_realdata_AI_panels\n\nContains individual PDF/SVG/PNG panels built from the uploaded real source data. Use panels_pdf or panels_svg for Adobe Illustrator assembly. The composite layout PDF/PNG is only a preview template.\n\nPanels:\nA Somatic UMAP\nB Somatic niche prioritization\nC LC/ST pathway enrichment\nD Integrated core somatic niche scores\nE Core LC-ST-germ association map\nF Donor-level niche axis PCA\nG Niche axis vs germline outcome\n',encoding='utf-8')
# Zip final
zip_path=BASE/'Figure3_final_realdata_AI_panels.zip'
if zip_path.exists(): zip_path.unlink()
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    for p in OUT.rglob('*'):
        if p.is_file(): z.write(p,p.relative_to(OUT.parent))
print('done', zip_path)
