"""FIND filtering rules, re-implemented from the Methods for the revision analyses.

A variant is screened if ClinVar (as annotated in gnomAD) calls it pathogenic or likely
pathogenic, or if it is a pLoF (stop gained, frameshift, splice acceptor, splice donor) on the gene's
canonical transcript, whatever its ClinVar classification. This matches FIND's app.py, which
checks pLoF consequences only on the canonical transcript reported by gnomAD. Per-group allele counts below the minimum are set to 0;
a variant is returned if one group has AF > min_af and >= fold x every other group."""
import pandas as pd, numpy as np, sys
POPS=['African/African American','Admixed American','Ashkenazi Jewish','East Asian','European (Finnish)','Middle Eastern','European (non-Finnish)','Amish','South Asian']
LOF={'stop_gained','frameshift_variant','splice_acceptor_variant','splice_donor_variant'}
def load(path):
    d=pd.read_csv(path, low_memory=False)
    cv=d['ClinVar Germline Classification'].fillna('')
    plp=cv.str.contains('athogenic') & ~cv.str.contains('onflicting')
    lof=d['VEP Annotation'].fillna('').apply(lambda s: any(x in s for x in LOF))
    d['PLP']=plp; d['LOF']=lof
    # canonical transcript = the transcript gnomAD reports for most variants in the gene
    canon=d['Transcript'].value_counts().idxmax()
    lof=lof&(d['Transcript']==canon)
    return d[plp|lof].copy()
def apply(d,mincount=5,fold=10,minaf=8e-5):
    ac=d[[f'Allele Count {p}' for p in POPS]].to_numpy(float)
    an=d[[f'Allele Number {p}' for p in POPS]].to_numpy(float)
    acz=np.where(ac<mincount,0,ac)
    af=np.divide(acz,an,out=np.zeros_like(acz),where=an>0)
    top=af.argmax(1); topaf=af.max(1)
    other=np.sort(af,1)[:,-2]
    ok=(topaf>minaf)&(topaf>=fold*other)
    return ok, np.array(POPS)[top], topaf, other
if __name__=='__main__':
    d=load(sys.argv[1]); ok,pop,taf,oth=apply(d)
    print(len(d),'P/LP or pLoF variants')
    h=d[ok].copy(); h['pop']=pop[ok]; h['AF']=taf[ok]; h['other']=oth[ok]
    print(h[['gnomAD ID','HGVS Consequence','VEP Annotation','ClinVar Germline Classification','Filters - joint','pop','AF','other']].to_string())
