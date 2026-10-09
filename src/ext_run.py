import pandas as pd,numpy as np,json,warnings; warnings.filterwarnings('ignore')
from rdkit import Chem, RDLogger, DataStructs; RDLogger.DisableLog('rdApp.*')
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score,average_precision_score
d=pd.read_csv('/tmp/glue/MolGlueDB_full.csv').drop_duplicates('StdInChIKey').reset_index(drop=True)
fps=[];desc=[];grp=[];keep=[];bv=[]
for i,s in enumerate(d.SMILES):
    m=Chem.MolFromSmiles(s)
    if m is None: continue
    fp=AllChem.GetMorganFingerprintAsBitVect(m,2,2048); a=np.zeros(2048,dtype=np.uint8); DataStructs.ConvertToNumpyArray(fp,a)
    try: sc=MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    except Exception: sc=''
    grp.append(sc if sc else f'nosc{i}'); fps.append(a); bv.append(fp); keep.append(i); desc.append([f(m) for _,f in Descriptors.descList])
d=d.loc[keep].reset_index(drop=True); g_sc=np.array(grp)
X=np.hstack([np.array(fps),np.nan_to_num(np.array(desc,dtype=np.float64),posinf=0,neginf=0).clip(-1e6,1e6)])
phys=d[['MolWeight','tPSA','clogP']].apply(pd.to_numeric,errors='coerce').fillna(0).values
n=len(d); R={'n':n}; print(R,flush=True)
# Butina clusters
dist=[]
for i in range(1,n): dist.extend(1-np.array(DataStructs.BulkTanimotoSimilarity(bv[i],bv[:i])))
cl=Butina.ClusterData(dist,n,0.40,isDistData=True); g_cl=np.zeros(n,int)
for ci,c in enumerate(cl):
    for j in c: g_cl[j]=ci
g_src=d.SourceAddress_Website.astype(str).values
R['n_scaffold_groups']=int(len(set(g_sc))); R['n_cluster_groups']=int(len(set(g_cl))); R['n_source_groups']=int(len(set(g_src))); print(R,flush=True)
rf=lambda:RandomForestClassifier(300,class_weight='balanced',n_jobs=2,random_state=0)
def oof(F,y,mask,g):
    idx=np.where(mask)[0]; o=np.zeros(len(idx)); folds=[]
    for tr,te in GroupKFold(5).split(idx,groups=g[idx]):
        o[te]=rf().fit(F[idx[tr]],y[idx[tr]]).predict_proba(F[idx[te]])[:,1]; folds.append((idx[tr],idx[te]))
    return idx,o,folds
y1=d.IsActive.values.astype(int); allm=np.ones(n,bool)
rp=d.RecruitingProtein.astype(str).str.strip(); m2=(rp!='Unknown').values; y2=(rp!='CRBN').values.astype(int)
def maxsim(folds):
    mx=[]
    for tr,te in folds:
        trb=[bv[i] for i in tr]
        for j in te: mx.append(max(DataStructs.BulkTanimotoSimilarity(bv[j],trb)))
    return float(np.mean(np.array(mx)>=0.85))
# reference, scaffold-cold, same pipeline
i,o,f=oof(X,y1,allm,g_sc); R['ref_scaffold_D1']=dict(auroc=roc_auc_score(y1[i],o),share_ge085=maxsim(f)); print(R['ref_scaffold_D1'],flush=True)
i,o,_=oof(X,y2,m2,g_sc); R['ref_scaffold_D2']=dict(auprc=average_precision_score(y2[i],o),prev=float(y2[i].mean())); print(R['ref_scaffold_D2'],flush=True)
# D8
i,o,f=oof(X,y1,allm,g_cl); a=roc_auc_score(y1[i],o); R['D8']=dict(auroc=a,share_ge085=maxsim(f),gate_met=bool(a>=0.70)); print('D8',R['D8'],flush=True)
# D9
i,o,f=oof(X,y1,allm,g_src); a=roc_auc_score(y1[i],o); _,op,_=oof(phys,y1,allm,g_src)
R['D9']=dict(auroc=a,auroc_physchem=roc_auc_score(y1[i],op),prevalence=float(y1.mean()),share_ge085=maxsim(f),gate_met=bool(a>=0.65)); print('D9',R['D9'],flush=True)
# D10
i,o,f=oof(X,y2,m2,g_src); ap=average_precision_score(y2[i],o); pv=float(y2[i].mean())
R['D10']=dict(auprc=ap,prevalence=pv,auroc=roc_auc_score(y2[i],o),gate_met=bool(ap>=pv+0.10)); print('D10',R['D10'],flush=True)
json.dump(R,open('results/ext_results.json','w'),indent=1,default=float); print('done')
