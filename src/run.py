import pandas as pd,numpy as np,json,warnings; warnings.filterwarnings('ignore')
from rdkit import Chem, RDLogger; RDLogger.DisableLog('rdApp.*')
from rdkit.Chem import AllChem, Descriptors, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score,average_precision_score,f1_score
D='/tmp/glue/MolGlueDB_full.csv'
d=pd.read_csv(D).drop_duplicates('StdInChIKey').reset_index(drop=True)
fps=[];desc=[];grp=[];keep=[];bv=[]
names=[n for n,_ in Descriptors.descList]
for i,s in enumerate(d.SMILES):
    m=Chem.MolFromSmiles(s)
    if m is None: continue
    fp=AllChem.GetMorganFingerprintAsBitVect(m,2,2048); a=np.zeros(2048,dtype=np.uint8); DataStructs.ConvertToNumpyArray(fp,a)
    try: sc=MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    except Exception: sc=''
    grp.append(sc if sc else f'nosc{i}'); fps.append(a); bv.append(fp); keep.append(i)
    desc.append([f(m) for _,f in Descriptors.descList])
d=d.loc[keep].reset_index(drop=True); F=np.array(fps); Dsc=np.nan_to_num(np.array(desc,dtype=np.float64),posinf=0,neginf=0).clip(-1e6,1e6); g=np.array(grp)
X=np.hstack([F,Dsc]); R={'n':len(d),'n_scaffolds':int(len(set(g)))}; print(R,flush=True)
phys=d[['MolWeight','tPSA','clogP']].apply(pd.to_numeric,errors='coerce').fillna(0).values
meta=pd.get_dummies(d[['Pharmacophore','Core','ResearchStage']].astype(str)).values.astype(float)
rf=lambda s=0:RandomForestClassifier(300,class_weight='balanced',n_jobs=2,random_state=s,min_samples_leaf=1)
def oof(Xm,y,mask,seed=0,proba=True,groups=None):
    idx=np.where(mask)[0]; o=np.zeros((len(idx),len(set(y[idx])))) if not proba else np.zeros(len(idx)); gg=g[idx]
    for tr,te in GroupKFold(5).split(idx,groups=gg):
        m=rf(seed).fit(Xm[idx[tr]],y[idx[tr]])
        if proba: o[te]=m.predict_proba(Xm[idx[te]])[:,1]
        else: o[te]=m.predict_proba(Xm[idx[te]])
    return idx,o
# D1
y1=d.IsActive.values.astype(int); allm=np.ones(len(d),bool)
i1,o1=oof(X,y1,allm); a_fp=roc_auc_score(y1[i1],o1)
_,op=oof(phys,y1,allm); a_phys=roc_auc_score(y1[i1],op)
_,om=oof(meta,y1,allm); a_meta=roc_auc_score(y1[i1],om)
R['D1']=dict(auroc_fp=round(a_fp,4),auroc_phys=round(a_phys,4),auroc_meta=round(a_meta,4),prevalence=round(float(y1.mean()),3),gate_met=bool(a_fp>=0.70 and a_fp-a_phys>=0.03)); print('D1',R['D1'],flush=True)
# D2
rp=d.RecruitingProtein.astype(str).str.strip(); m2=(rp!='Unknown').values; y2=(rp!='CRBN').values.astype(int)
i2,o2=oof(X,y2,m2); ap=average_precision_score(y2[i2],o2); prev=float(y2[i2].mean())
_,o2m=oof(meta,y2,m2); apm=average_precision_score(y2[i2],o2m)
R['D2']=dict(auprc=round(ap,4),prevalence=round(prev,4),auprc_meta=round(apm,4),auroc=round(roc_auc_score(y2[i2],o2),4),gate_met=bool(ap>=prev+0.20)); print('D2',R['D2'],flush=True)
# D3
pt=d.PrimaryTarget.astype(str); top=pt.value_counts().index[:5].tolist(); m3=pt.isin(top).values; cls=sorted(top); y3=pt.map({c:i for i,c in enumerate(cls)}).fillna(-1).astype(int).values
idx=np.where(m3)[0]; pred=np.zeros(len(idx),int); gg=g[idx]
for tr,te in GroupKFold(5).split(idx,groups=gg):
    m=RandomForestClassifier(300,class_weight='balanced',n_jobs=2,random_state=0).fit(X[idx[tr]],y3[idx[tr]]); pred[te]=m.predict(X[idx[te]])
f1=f1_score(y3[idx],pred,average='macro'); maj=np.full(len(idx),np.bincount(y3[idx]).argmax()); f1m=f1_score(y3[idx],maj,average='macro')
R['D3']=dict(classes=top,n=int(len(idx)),macroF1=round(f1,4),majority_macroF1=round(f1m,4),gate_met=bool(f1>=f1m+0.15)); print('D3',R['D3'],flush=True)
# D4
sh=[]
for s in range(5):
    rng=np.random.RandomState(s); ys=y1.copy(); 
    # shuffle labels within each fold's train set only via permuting globally across groups (preserves prevalence)
    ys=rng.permutation(y1); _,os_=oof(X,ys,allm,seed=s); sh.append(roc_auc_score(ys[i1],os_))
# near-duplicate stats
mx=[];
for tr,te in GroupKFold(5).split(np.arange(len(d)),groups=g):
    trb=[bv[i] for i in tr]
    for j in te: mx.append(max(DataStructs.BulkTanimotoSimilarity(bv[j],trb)))
mx=np.array(mx)
R['D4']=dict(shuffle_auroc_mean=round(float(np.mean(sh)),4),shuffle_aurocs=[round(x,4) for x in sh],mean_max_tanimoto=round(float(mx.mean()),3),share_ge_085=round(float((mx>=0.85).mean()),3),gate_met=bool(0.45<=np.mean(sh)<=0.55)); print('D4',R['D4'],flush=True)
# D5
R['D5']=dict(delta_fp_minus_meta_D1=round(a_fp-a_meta,4),delta_fp_minus_meta_D2_auprc=round(ap-apm,4),gate_met=bool(a_fp-a_meta>=0.03))
# D6
ss=[a_fp]
for s in range(1,5):
    _,o=oof(X,y1,allm,seed=s); ss.append(roc_auc_score(y1[i1],o))
R['D6']=dict(aurocs=[round(x,4) for x in ss],std=round(float(np.std(ss)),4),gate_met=bool(np.std(ss)<=0.02)); print('D6',R['D6'],flush=True)
json.dump(R,open('results/results.json','w'),indent=1)
