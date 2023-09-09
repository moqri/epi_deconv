import pandas as pd
atlas='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186458/matrix/GSE186458_series_matrix.txt.gz'
df=pd.read_table(atlas,skiprows=79,index_col=0)
df.T['!Sample_supplementary_file_3'].to_csv('ftps',index=False,header=False)
