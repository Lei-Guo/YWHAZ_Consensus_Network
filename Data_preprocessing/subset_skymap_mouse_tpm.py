import numpy as np
import pandas as pd

querySpecies = "Mus_musculus"
expression_metric = "tpm"

def loadDf(fname,mmap_mode='r'):
    with open(fname+'.index.txt') as f:
        myIndex=map(lambda s:s.replace("\n",""), f.readlines())
    with open(fname+'.columns.txt') as f:
        myColumns=map(lambda s:s.replace("\n",""), f.readlines())
    tmpMatrix=np.load(fname+".npy",mmap_mode=mmap_mode)
    tmpDf=pd.DataFrame(tmpMatrix,index=myIndex,columns=myColumns)
    tmpDf.columns.name='Run'
    return tmpDf

baseDir='/sc/orga/projects/zhangb03a/shared/skymap/' #Base directory
expression_metric='tpm' #offer Kallisto expression metric: ["tpm","est_counts"]
data_matrix_dir=baseDir+'/{specie}.gene_symbol.{expression_metric}'.format(specie=querySpecies, expression_metric=expression_metric)

rnaseqDf=loadDf(data_matrix_dir)

# load my curated meta file
my_meta_mouse = pd.read_csv('/sc/orga/projects/zhangb03a/lei_guo/YWHAZ/Datasets/meta_GSE66716_GSE60745.csv')

# filter the correct columns
filteredDf = rnaseqDf.loc[:,set(rnaseqDf.columns).intersection(set(my_meta_mouse['SRA_Run_ID']))]
print("Columns in meta file",len(my_meta_mouse['SRA_Run_ID']))
print("Columns in Skymap", len(filteredDf.columns))

# drop genes with 0 expression
print("Rows before",len(filteredDf))
cutoff = 0
#filteredDf = filteredDf[~(filteredDf<=cutoff).all(axis=1)]

# drop genes with 0 expression in any samples
filteredDf = filteredDf[(filteredDf != cutoff).all(1)]
print("Rows after",len(filteredDf))

# save output
filteredDf.to_csv("/sc/orga/projects/zhangb03a/lei_guo/YWHAZ/Datasets/YWHAZ_KO_KD_RNAseq_TPM_GSE66716_GSE60745.tsv", sep = "\t")