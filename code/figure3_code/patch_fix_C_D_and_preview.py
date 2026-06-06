from pathlib import Path
import zipfile, shutil
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from PIL import Image, ImageDraw, ImageFont

BASE=Path('/mnt/data/Figure3_final_realdata_AI_panels')
DATA=BASE/'figure3_source_data'
PDF=BASE/'panels_pdf'; SVG=BASE/'panels_svg'; PNG=BASE/'panels_png'; PREV=BASE/'preview'
for d in [PDF,SVG,PNG,PREV]: d.mkdir(exist_ok=True)
# Formatting
mpl.rcParams['pdf.fonttype']=42; mpl.rcParams['ps.fonttype']=42; mpl.rcParams['svg.fonttype']='none'; mpl.rcParams['font.family']='DejaVu Sans'
cell_cols={'Leydig cells':'#ff7f0e','Sertoli cells':'#8c564b'}
group_cols={'Ctrl':'#1f77b4','OA':'#d62728','AZFc_Del':'#9467bd','iNOA_B':'#ff7f0e','iNOA_S':'#2ca02c','KS':'#8c564b'}
group_order=['Ctrl','OA','AZFc_Del','iNOA_B','iNOA_S','KS']
disease_order=['OA','AZFc_Del','iNOA_B','iNOA_S','KS']

def savefig(fig,name):
    fig.savefig(PDF/f'{name}.pdf', bbox_inches='tight')
    fig.savefig(SVG/f'{name}.svg', bbox_inches='tight')
    fig.savefig(PNG/f'{name}.png', dpi=500, bbox_inches='tight')
    plt.close(fig)

def clean(ax):
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# Panel C fixed
panelC_src=pd.read_csv(DATA/'Fig3C_LC_ST_pathway_enrichment_selected.csv')
lc_map={'Cholesterol\nhomeostasis':'Cholesterol Homeostasis','Fatty acid\nmetabolism':'Fatty Acid Metabolism','OXPHOS':'Oxidative Phosphorylation','TNF-alpha/\nNF-kB':'Tnfa Signaling Via Nfkb','ECM / adhesion':'Epithelial Mesenchymal Transition','Hypoxia':'Hypoxia','Apoptosis':'Apoptosis'}
st_map={'ECM / adhesion':'Epithelial Mesenchymal Transition','TNF-alpha/\nNF-kB':'Tnfa Signaling Via Nfkb','IFN-gamma\nresponse':'Interferon Gamma Response','Apoptosis':'Apoptosis','p53 stress':'P53 Pathway','Glycolysis':'Glycolysis','OXPHOS':'Oxidative Phosphorylation'}

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

fig=plt.figure(figsize=(12.8,3.6))
gs=GridSpec(1,4,width_ratios=[1.05,1.05,0.22,0.06],wspace=0.32,figure=fig)
ax1=fig.add_subplot(gs[0,0]); ax2=fig.add_subplot(gs[0,1]); axl=fig.add_subplot(gs[0,2]); cax=fig.add_subplot(gs[0,3])
pathway_ax(ax1,'LC',lc_map,'LC pathways',cell_cols['Leydig cells'])
pathway_ax(ax2,'ST',st_map,'ST pathways',cell_cols['Sertoli cells'])
# dedicated size legend axis
axl.axis('off')
axl.text(0.0,0.84,'-log10(FDR)',fontsize=8,ha='left',va='center')
for y,s,lab in zip([0.65,0.48,0.28],[16+2*18,16+4*18,16+8*18],['2','4','8']):
    axl.scatter(0.22,y,s=s,color='white',edgecolor='black')
    axl.text(0.52,y,lab,fontsize=8,va='center')
axl.set_xlim(0,1); axl.set_ylim(0,1)
sm=mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=-3,vmax=3),cmap='RdBu_r')
cb=fig.colorbar(sm,cax=cax); cb.set_label('NES',fontsize=9)
ax1.text(-0.18,1.15,'c',transform=ax1.transAxes,fontsize=16,fontweight='bold')
fig.suptitle('LC and ST disease-associated pathways',y=1.02,fontsize=12,fontweight='bold')
savefig(fig,'Fig3C_LC_ST_pathway_enrichment')

# Panel D fixed
score_core=pd.read_csv(DATA/'Fig3D_core_somatic_niche_scores_z_matrix.csv').set_index('score_display')
sample_order=pd.read_csv(DATA/'Fig3D_sample_order_by_PC1v2.csv')
ordered_samples=sample_order.sort_values('sample_rank_pc1v2')['Sample_ID'].tolist()
score_core=score_core[ordered_samples]
group_bar=pd.read_csv(DATA/'Fig3D_group_bar_annotation.csv')
fig=plt.figure(figsize=(12.8,3.8))
gs=GridSpec(3,2,height_ratios=[0.28,0.17,1.0],width_ratios=[1,0.025],hspace=0.08,wspace=0.03,figure=fig)
ax_title=fig.add_subplot(gs[0,0]); ax_title.axis('off')
ax_group=fig.add_subplot(gs[1,0]); ax=fig.add_subplot(gs[2,0]); cax=fig.add_subplot(gs[2,1])
ax_title.text(0.5,0.72,'Integrated somatic niche scores',ha='center',va='center',fontsize=12,fontweight='bold')
ax_title.text(0.00,0.28,'preserved',color='#2b6cb0',fontsize=9,ha='left',va='center')
ax_title.annotate('niche-collapse axis',xy=(0.83,0.28),xytext=(0.33,0.28),xycoords='axes fraction',textcoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=1),fontsize=9,ha='center',va='center')
ax_title.text(1.00,0.28,'collapsed',color='#c53030',fontsize=9,ha='right',va='center')
handles=[Line2D([0],[0],marker='s',linestyle='',color='none',markerfacecolor=group_cols[g],markersize=7,label=g) for g in group_order]
ax_title.legend(handles=handles,ncol=6,frameon=False,fontsize=7,loc='center right',bbox_to_anchor=(1.0,0.82))
ax_title.text(-0.08,0.98,'d',transform=ax_title.transAxes,fontsize=16,fontweight='bold',va='top')
# group bar
sample_groups=group_bar.set_index('Sample_ID').loc[ordered_samples,'Group']
rgb=np.array([mpl.colors.to_rgb(group_cols[g]) for g in sample_groups])[None,:,:]
ax_group.imshow(rgb,aspect='auto'); ax_group.set_xticks([]); ax_group.set_yticks([])
for sp in ax_group.spines.values(): sp.set_visible(False)
# heatmap
im=ax.imshow(score_core.values,cmap='RdBu_r',vmin=-2.3,vmax=2.3,aspect='auto')
ax.set_yticks(np.arange(score_core.shape[0])); ax.set_yticklabels(score_core.index.tolist(),fontsize=9)
ax.set_xticks(np.arange(len(ordered_samples))); ax.set_xticklabels(ordered_samples,rotation=45,ha='right',fontsize=7)
ax.set_xlabel('Donors ordered by niche-collapse axis',fontsize=9)
ax.hlines([3.5],-0.5,len(ordered_samples)-0.5,color='white',lw=1.5)
cb=fig.colorbar(im,cax=cax); cb.set_label('Z-score',fontsize=9)
savefig(fig,'Fig3D_integrated_somatic_niche_scores')

# Rebuild composite preview with updated panel C/D PNGs.
# Use existing panel pngs; if they are missing, package likely incomplete.
imgs={p.stem:Image.open(p).convert('RGB') for p in PNG.glob('Fig3*.png')}
def rw(img,w):
    sc=w/img.width
    return img.resize((int(w),int(img.height*sc)))
W=2400; margin=60; gap=36
A=rw(imgs['Fig3A_somatic_UMAP'],720); B=rw(imgs['Fig3B_somatic_niche_prioritization'],1500)
C=rw(imgs['Fig3C_LC_ST_pathway_enrichment'],W-2*margin)
D=rw(imgs['Fig3D_integrated_somatic_niche_scores'],W-2*margin)
E=rw(imgs['Fig3E_core_LC_ST_germ_association_map'],700); F=rw(imgs['Fig3F_donor_level_niche_axis_PCA'],560); G=rw(imgs['Fig3G_niche_axis_predicts_germline_outcome'],920)
H=90+max(A.height,B.height)+gap+C.height+gap+D.height+gap+max(E.height,F.height,G.height)+margin
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
canvas.save(PDF/'Figure3_composite_layout_preview.pdf','PDF',resolution=300)
# Rezip package
zip_path=BASE.parent/'Figure3_final_realdata_AI_panels.zip'
# actually BASE.parent is /mnt/data; good
if zip_path.exists(): zip_path.unlink()
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    for p in BASE.rglob('*'):
        if p.is_file() and '_tmp' not in p.parts:
            z.write(p,p.relative_to(BASE.parent))
print('patched and zipped',zip_path)
