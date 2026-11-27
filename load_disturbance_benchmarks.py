import json, numpy as np, sys

from safehx.loadsim import cl_with_load
from scipy.linalg import expm
from safehx.model import *
from safehx.benchplants import setup, static_gains, ideal_dynamic, pade
def iae_load(cfg, ch, H, dt):
    A,B,C,D = cl_with_load(cfg); n=A.shape[0]
    Mx=np.zeros((n+1,n+1)); Mx[:n,:n]=A*dt; Mx[:n,n]=B[:,ch]*dt; E=expm(Mx); Ad,bd=E[:n,:n],E[:n,n]
    x=np.zeros(n); s=0.0
    for k in range(int(H/dt)):
        y=C[:2]@x+D[:2,ch]; s+=np.sum(np.abs(y))*dt; x=Ad@x+bd
    return s
out={}
for name in ["Hylte","WoodBerry"]:
    base,P=setup(name); KF1,KF2=static_gains(name); d=json.load(open(f"results/bench2_{name}.json"))
    dd=1-KF1*KF2
    cf={"none":base,
        "static_inv":base.with_(K11=replace(base.K11,Kp=base.K11.Kp/dd),K22=replace(base.K22,Kp=base.K22.Kp/dd),D12_tf=(np.array([KF1]),np.array([1.0])),D21_tf=(np.array([KF2]),np.array([1.0]))),
        "pid":base.with_(D21=PIDCompact.from_vector(d["pid"][:6]),D12=PIDCompact.from_vector(d["pid"][6:]))}
    if d["ideal_dynamic_realizable"]:
        dyn={}
        for w in ["D12","D21"]:
            num,den,L=ideal_dynamic(name,w); pn,pd=pade(L); dyn[w+"_tf"]=(np.polymul(num,pn),np.polymul(den,pd))
        cf["ideal"]=base.with_(**dyn)
    r={k:(iae_load(c,2,P["H"],P["dt"])+iae_load(c,3,P["H"],P["dt"])) for k,c in cf.items()}
    out[name]={k:v/r["none"] for k,v in r.items()}
    print(name,{k:round(v,3) for k,v in out[name].items()})
json.dump(out,open("results/bench2_load.json","w"),indent=1)
