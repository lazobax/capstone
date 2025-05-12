from graphviz import Digraph

# Create a detailed flowchart with commented lines included
dot = Digraph(format='svg', engine='dot')
dot.attr(label='Comprehensive VCF Processing Pipeline', labelloc='t', fontsize='14', fontname='Helvetica')
dot.attr(rankdir='TB', splines='ortho', nodesep='0.3', ranksep='0.5')

# Cluster: Setup & Configuration
with dot.subgraph(name='cluster_config') as c:
    c.attr(label='Setup & Configuration', color='gray')
    c.node('Shebang', '#!/usr/bin/env bash\nset -euo pipefail', shape='note')
    c.node('Mapping', '# Default gene list & chromosome mapping', shape='note')
    c.node('Override', '# Allow override of gene list via args', shape='note')
    c.edge('Shebang', 'Mapping', style='invis')
    c.edge('Mapping', 'Override', style='invis')

# Cluster: Initial Downloads
with dot.subgraph(name='cluster_download') as d:
    d.attr(label='Initial Downloads', color='gray')
    d.node('Dirs', 'mkdir -p data ...', shape='box')
    d.node('ClinVarDL', '# Download & prepare ClinVar summary', shape='box')
    d.node('RefGRCh37', '# Download & prepare GRCh37 reference', shape='box')
    d.edge('Dirs', 'ClinVarDL')
    d.edge('ClinVarDL', 'RefGRCh37')

# Cluster: Gene Processing Loop
with dot.subgraph(name='cluster_genes') as g:
    g.attr(label='Per-Gene Processing Loop', color='gray')
    g.node('LoopStart', 'for gene in GENES', shape='box')
    g.node('ClinVarToVcf', '# ClinVar → VCF', shape='box')
    g.node('AddID', 'python3 scripts/addClinVarID.py', shape='box')
    g.node('LOVDGet', '# LOVD scraping and parsing', shape='box')
    g.edge('LoopStart', 'ClinVarToVcf')
    g.edge('ClinVarToVcf', 'AddID')
    g.edge('AddID', 'LOVDGet')

# Cluster: VCF Reformat
with dot.subgraph(name='cluster_vcf_reformat') as r:
    r.attr(label='VCF Reformat', color='gray')
    r.node('MakeDirs', '# Reformating into VCF format', shape='note')
    r.node('ToVcfCV', 'python3 scripts/tsvToVcf.py --indir ClinVar_vcf', shape='box')
    r.node('ToVcfLOVD', 'python3 scripts/tsvToVcf.py --indir LOVD_vcf', shape='box')
    r.edge('MakeDirs', 'ToVcfCV')
    r.edge('ToVcfCV', 'ToVcfLOVD')

# Cluster: Contig & Directories
with dot.subgraph(name='cluster_contigs') as ct:
    ct.attr(label='Contigs & Output Dirs', color='gray')
    ct.node('BuildContigs', '# Build contigs.txt from genome.fa.fai', shape='box')
    ct.node('PrepareVCFDirs', '# Prepare VCF/LOVD & VCF/CV dirs', shape='box')
    ct.edge('BuildContigs', 'PrepareVCFDirs')

# Cluster: VCF Cleanup & Stats Loop
with dot.subgraph(name='cluster_vcf_loop') as vl:
    vl.attr(label='VCF Clean & Stats Loop', color='gray')
    vl.node('InitLog', '# Initialize logfile lines.log', shape='note')
    vl.node('LoopVCF', 'for each sub in LOVD,CV', shape='box')
    vl.node('CountBefore', '# 1) count before', shape='box')
    vl.node('ReplaceCHR', '# a) replaceCHR → contig names', shape='box')
    vl.node('CleanAmbig', '# b) clean IUPAC & mismatches', shape='box')
    vl.node('BgzipClean', '# c) bgzip cleaned VCF', shape='box')
    vl.node('ReheaderContig', '# d) inject contig headers', shape='box')
    vl.node('SortCompress', '# e) sort & re-compress', shape='box')
    vl.node('IndexFinal', '# f) index (final)', shape='box')
    vl.node('LogStats', '# append results to log', shape='box')
    vl.edge('InitLog', 'LoopVCF')
    vl.edge('LoopVCF', 'CountBefore')
    vl.edge('CountBefore', 'ReplaceCHR')
    vl.edge('ReplaceCHR', 'CleanAmbig')
    vl.edge('CleanAmbig', 'BgzipClean')
    vl.edge('BgzipClean', 'ReheaderContig')
    vl.edge('ReheaderContig', 'SortCompress')
    vl.edge('SortCompress', 'IndexFinal')
    vl.edge('IndexFinal', 'LogStats')

# Cluster: Concatenation & Combined Processing
with dot.subgraph(name='cluster_concat') as cc:
    cc.attr(label='Concatenate & Combined Processing', color='gray')
    cc.node('Concat', '# 4) concatenate all sorted VCFs', shape='box')
    cc.node('AfterConcatEcho', 'echo "AFTER CONCAT"', shape='box')
    cc.node('CountCombined', 'count variants after concat', shape='box')
    cc.node('RemoveAmbigCombined', '# remove non-ATCGN bases', shape='box')
    cc.node('BgzipCombined', 'bgzip combined.clean.vcf', shape='box')
    cc.node('NormCombined', 'bcftools norm combined.clean.vcf.gz', shape='box')
    cc.node('IndexCombined', 'tabix combined.normalized.vcf.gz', shape='box')
    cc.edge('Concat', 'AfterConcatEcho')
    cc.edge('AfterConcatEcho', 'CountCombined')
    cc.edge('CountCombined', 'RemoveAmbigCombined')
    cc.edge('RemoveAmbigCombined', 'BgzipCombined')
    cc.edge('BgzipCombined', 'NormCombined')
    cc.edge('NormCombined', 'IndexCombined')

# Cluster: ANNOVAR Annotation
with dot.subgraph(name='cluster_annovar') as ca:
    ca.attr(label='ANNOVAR Annotation', color='gray')
    ca.node('AnnovarComment', '# Download ANNOVAR dbs (commented)', shape='note')
    ca.node('RenameHeader', '# Rename header & CHROM (awk)', shape='box')
    ca.node('TabixNum', 'bgzip & tabix num vcf', shape='box')
    ca.node('CountNumVCF', 'echo count after renaming', shape='box')
    ca.node('TableAnnovar', 'perl annovar/table_annovar.pl', shape='box')
    ca.node('MoveAnno', 'mv multianno.txt → finalAnnot.tsv', shape='box')
    ca.edge('AnnovarComment', 'RenameHeader')
    ca.edge('RenameHeader', 'TabixNum')
    ca.edge('TabixNum', 'CountNumVCF')
    ca.edge('CountNumVCF', 'TableAnnovar')
    ca.edge('TableAnnovar', 'MoveAnno')

# Cluster: Final Moves
with dot.subgraph(name='cluster_moves') as mv:
    mv.attr(label='Final File Moves & Prep', color='gray')
    mv.node('MakeCombinedDir', 'mkdir -p data/combined/', shape='box')
    mv.node('MoveAll', 'mv data/combined.* data/combined/', shape='box')
    mv.node('RenameFinalVCF', 'mv combined.normalized.vcf → finalVCF.vcf', shape='box')
    mv.edge('MakeCombinedDir', 'MoveAll')
    mv.edge('MoveAll', 'RenameFinalVCF')

# Cluster: LiftOver to hg38
with dot.subgraph(name='cluster_liftover') as lo:
    lo.attr(label='LiftOver to hg38', color='gray')
    lo.node('ChainDL', '# 1) download & unpack chain file', shape='box')
    lo.node('Hg38RefDL', '# 2) download & index hg38 reference', shape='box')
    lo.node('CrossMap', '# 4) liftover VCF via CrossMap.py', shape='box')
    lo.edge('ChainDL', 'Hg38RefDL')
    lo.edge('Hg38RefDL', 'CrossMap')

# Cluster: Post-LiftOver Processing
with dot.subgraph(name='cluster_post') as pp:
    pp.attr(label='Post-LiftOver Processing', color='gray')
    pp.node('Reheader38', 'bcftools reheader --fai hg38.fa.fai', shape='box')
    pp.node('AddChr', 'awk add "chr" prefix', shape='box')
    pp.node('Bgzip38', 'bgzip final.hg38.nums.vcf', shape='box')
    pp.node('Norm38', 'bcftools norm -f hg38.fa', shape='box')
    pp.edge('Reheader38', 'AddChr')
    pp.edge('AddChr', 'Bgzip38')
    pp.edge('Bgzip38', 'Norm38')

# Link clusters sequentially
dot.edge('Mapping', 'Dirs', style='dashed')
dot.edge('RefGRCh37', 'LoopStart', style='dashed')
dot.edge('LOVDGet', 'MakeDirs', style='dashed')
dot.edge('ToVcfLOVD', 'BuildContigs', style='dashed')
dot.edge('PrepareVCFDirs', 'InitLog', style='dashed')
dot.edge('LogStats', 'Concat', style='dashed')
dot.edge('IndexCombined', 'RenameHeader', style='dashed')
dot.edge('MoveAnno', 'MakeCombinedDir', style='dashed')
dot.edge('RenameFinalVCF', 'ChainDL', style='dashed')
dot.edge('Norm38', 'End', style='invis')

# Final node
dot.node('End', 'All Done', shape='oval')

# Render SVG
svg = dot.pipe().decode('utf-8')
with open('/mnt/data/comprehensive_pipeline.svg', 'w') as f:
    f.write(svg)

print("Comprehensive pipeline flowchart generated: [Download SVG](/mnt/data/comprehensive_pipeline.svg)")
