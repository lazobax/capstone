import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

# ———————————————————————————————————————————
def parse_attrs(attr_str):
    return {k: v for k,v in (kv.split('=',1) for kv in attr_str.split(';') if '=' in kv)}

def get_transcript_and_exons(gff_path, gene_name):
    """Locate the gene feature by Name, then its first mRNA child, then its exons."""
    gene_id = None
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#'): continue
            feat, attrs = line.split('\t')[2], line.split('\t')[8]
            if feat=='gene' and parse_attrs(attrs).get('Name')==gene_name:
                gene_id = parse_attrs(attrs)['ID']
                break
    if gene_id is None:
        raise ValueError(f"No gene Name={gene_name}")

    transcript_id = tx_start = tx_end = None
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#'): continue
            cols = line.split('\t')
            attrs_d = parse_attrs(cols[8])
            if cols[2]=='mRNA' and attrs_d.get('Parent')==gene_id:
                transcript_id = attrs_d['ID']
                tx_start, tx_end = int(cols[3]), int(cols[4])
                break
    if transcript_id is None:
        raise ValueError(f"No mRNA under gene {gene_id}")

    exons = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#'): continue
            cols = line.split('\t')
            attrs_d = parse_attrs(cols[8])
            if cols[2]=='exon' and attrs_d.get('Parent')==transcript_id:
                exons.append((int(cols[3]), int(cols[4])))
    if not exons:
        raise ValueError(f"No exons for transcript {transcript_id}")

    exons_rel = sorted([(s - tx_start, e - tx_start) for s,e in exons], key=lambda x: x[0])
    return tx_start, tx_end, exons_rel

def group_sig(s):
    sl = str(s).lower()
    if 'pathogenic' in sl: return 'Pathogenic'
    if 'benign' in sl:     return 'Benign'
    return 'Unknown'

def collapse_introns(exons_rel, intron_shrink):
    segments, offset = [], 0
    for i,(s,e) in enumerate(exons_rel):
        segments.append({'type':'exon','orig_s':s,'orig_e':e,'new_s':offset,'new_e':offset+(e-s)})
        offset += (e-s)
        if i < len(exons_rel)-1:
            next_s = exons_rel[i+1][0]
            segments.append({'type':'intron','orig_s':e,'orig_e':next_s,
                             'new_s':offset,'new_e':offset+intron_shrink})
            offset += intron_shrink
    return segments, offset

def collapse_pos(p, segments, total_width):
    for seg in segments:
        if seg['orig_s'] <= p < seg['orig_e']:
            orig_len = seg['orig_e'] - seg['orig_s']
            frac = (p - seg['orig_s'])/orig_len if orig_len>0 else 0
            return seg['new_s'] + frac*(seg['new_e']-seg['new_s'])
    return total_width
# ———————————————————————————————————————————

# 0) load & pre-filter merged_df only once
merged_tsv = "merged.tsv"
df_all = pd.read_csv(merged_tsv, sep='\t', comment='#')
df_all = df_all[
    (df_all['Func.refGene']=='exonic') &
    (df_all['ExonicFunc.refGene']=='nonsynonymous SNV')
].dropna(subset=['POS','CONSEQ']).copy()
df_all['Pos'] = pd.to_numeric(df_all['POS'], errors='coerce').astype('Int64').astype(int)

# plotting parameters
colors      = {'Pathogenic':'#FF6961','Benign':'#8CD47E','Unknown':'grey'}
y_map       = {'Pathogenic':1.0,'Unknown':0.6,'Benign':0.2}
type_marker = {'SNP':'o','DEL':'s','INS':'^','Microsatellite':'D'}
intron_shrink = 100
alpha         = 0.15
marker_size   = 6
gff_path      = "data/GRCh37.gff"
output_dir    = "plots3"

genes = ["FANCA","FANCB","FANCC","BRCA2","FANCD2","FANCE",
         "FANCF","FANCG","FANCI","BRIP1","FANCL","FANCM",
         "PALB2","RAD51C","SLX4","ERCC4","RAD51","BRCA1",
         "UBE2T","XRCC2","MAD2L2","RFWD3"]

for gene in genes:
    tx_start, tx_end, exons_rel = get_transcript_and_exons(gff_path, gene)
    tx_len = tx_end - tx_start

    # 1) slice this gene’s mutations
    df = df_all[df_all['Gene.refGene']==gene].copy()  # adjust column name if needed
    df['rel_pos'] = df['Pos'] - tx_start
    df = df[(df['rel_pos']>=0)&(df['rel_pos']<=tx_len)]
    if df.empty:
        print(f"No exonic nonsyn SNVs for {gene}, skipping.")
        continue

    # 2) classify and collapse
    df['group'] = df['CONSEQ'].apply(group_sig)
    segments, total_width = collapse_introns(exons_rel, intron_shrink)
    df['x'] = df['rel_pos'].apply(lambda p: collapse_pos(p, segments, total_width))

    # 3) plot
    fig, ax = plt.subplots(figsize=(12,3))
    for seg in segments:
        if seg['type']=='exon':
            ax.add_patch(Rectangle((seg['new_s'],-0.05),
                                   seg['new_e']-seg['new_s'],0.1,
                                   color='black'))
    for _, row in df.iterrows():
        ax.plot(row['x'], y_map[row['group']],
                marker=type_marker.get(row.get('Type',''),'o'),
                color=colors[row['group']],
                markersize=marker_size, alpha=alpha)

    # x‐axis ticks
    ticks = [seg['new_s'] for seg in segments] + [total_width]
    labels = [str(int(seg['orig_s']+tx_start)) for seg in segments] + [str(tx_end)]
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, rotation=90, ha='right', fontsize=3)
    ax.set_xlim(0, total_width)
    ax.set_ylim(0, 1.05)
    ax.set_yticks([])
    ax.set_xlabel('Collapsed transcript coordinate')
    ax.set_title(f"Mutations on {gene}")

    # legend outside (colors only)
    color_legs = [Line2D([0],[0], color=c, lw=2, label=g)
                  for g,c in colors.items()]
    ax.legend(handles=color_legs,
              loc='upper left',
              bbox_to_anchor=(1.01,1),
              borderaxespad=0,
              title="Significance",
              fontsize=8)
    plt.tight_layout()
    plt.subplots_adjust(right=0.8)

    out = f"{output_dir}/{gene}_mutation_plot.png"
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved plot for {gene} to {out}")
