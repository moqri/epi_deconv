import pandas as pd
epic='https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC/EPIC.hg38.manifest.tsv.gz'
df=pd.read_table(epic,index_col='Probe_ID',usecols=['Probe_ID','CpG_chrm','CpG_beg','CpG_end'])
df=df.dropna().copy()
df[['CpG_beg','CpG_end']]=df[['CpG_beg','CpG_end']].astype(int)
dna='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110554/matrix/GSE110554_series_matrix.txt.gz'
dna=pd.read_table(dna,skiprows=97,index_col=0)
dna=dna.round(3).copy()
dg=df.merge(dna,left_index=True,right_index=True)
dg.sort_values(['CpG_chrm','CpG_beg']).to_csv('dna.csv',index=False,sep='\t')
