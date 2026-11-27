import numpy as np
from .model import *
from .model import _controller_ss, plant_ss, PLANT_TF
def cl_with_load(cfg):
    Ag,Bg,Cg,Dg=plant_ss(cfg.plant or PLANT_TF); Ak,Bk,Ck,Dk=_controller_ss(cfg)
    Bkr,Bky=Bk[:,:2],Bk[:,2:]; Dkr,Dky=Dk[:,:2],Dk[:,2:]; ng,nk=Ag.shape[0],Ak.shape[0]
    M=np.linalg.inv(np.eye(2)-Dg@Dky)
    # inputs [r1,r2,d1,d2]; plant input = u + d
    Cy=M@np.hstack([Cg,Dg@Ck]); Dyr=M@np.hstack([Dg@Dkr, Dg])
    Cu=np.hstack([np.zeros((2,ng)),Ck])+Dky@Cy; Dur=np.hstack([Dkr,np.zeros((2,2))])+Dky@Dyr
    A=np.block([[Ag,np.zeros((ng,nk))],[np.zeros((nk,ng)),Ak]])+np.vstack([Bg@Cu,Bky@Cy])
    B=np.vstack([Bg@(Dur+np.hstack([np.zeros((2,2)),np.eye(2)])), np.hstack([Bkr,np.zeros((nk,2))])+Bky@Dyr])
    return A,B,np.vstack([Cy,Cu]),np.vstack([Dyr,Dur])
