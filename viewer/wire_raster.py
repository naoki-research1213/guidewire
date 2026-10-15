"""Sparse, analytical ray intervals; only projected wire bounding boxes are visited."""
import math
import numpy as np


def capsule_intervals(source, rays, start, end, radius, tmax):
    """Finite cylinder plus endpoint spheres, clipped to source/detector."""
    axis = end-start
    length = np.linalg.norm(axis)
    axis = axis/length
    m = source-start
    de = rays@axis
    me = float(m@axis)
    dp = rays-de[:,None]*axis
    mp = m-me*axis
    qa = np.sum(dp*dp,axis=1)
    qb = 2*(dp@mp)
    qc = float(mp@mp)-radius**2
    disc = qb**2-4*qa*qc
    regular = qa>1e-14
    root = np.sqrt(np.maximum(disc,0))
    lo = np.where(regular,(-qb-root)/np.maximum(2*qa,1e-14),-np.inf)
    hi = np.where(regular,(-qb+root)/np.maximum(2*qa,1e-14),np.inf)
    valid = np.where(regular,disc>=0,qc<=0)
    moving = np.abs(de)>1e-14
    t0 = -me/np.where(moving,de,1.)
    t1 = (length-me)/np.where(moving,de,1.)
    lo = np.maximum(lo,np.where(moving,np.minimum(t0,t1),-np.inf))
    hi = np.minimum(hi,np.where(moving,np.maximum(t0,t1),np.inf))
    valid &= moving | ((me>=0)&(me<=length))
    yield np.maximum(lo,0), np.minimum(hi,tmax), valid
    for center in (start,end):
        m=source-center
        b=rays@m
        disc=b*b-(float(m@m)-radius**2)
        root=np.sqrt(np.maximum(disc,0))
        yield np.maximum(-b-root,0),np.minimum(-b+root,tmax),disc>=0


def wire_transmission(points,camera,matrix,diameter,mu,depth_mu):
    shape=(matrix,matrix,4)
    tau=np.zeros(shape,dtype=np.float32)
    if len(points)<2:
        return np.ones((matrix,matrix),dtype=np.float32)
    projected=camera.project_mm(points)
    offsets=np.array([[-.25,-.25],[.25,-.25],[-.25,.25],[.25,.25]])
    groups,starts,ends=[],[],[]
    radius=diameter*.5
    for i,(start,end) in enumerate(zip(points[:-1],points[1:])):
        if np.linalg.norm(end-start)<=1e-10:
            continue
        a,b=projected[i],projected[i+1]
        # Conservative sphere extent, including off-axis cone projection.
        depth=max(min(a[2],b[2])-radius,1.)
        r=radius*camera.source_detector_mm/depth
        rx=r*(1+max(abs(a[0]),abs(b[0]))/camera.source_detector_mm)
        ry=r*(1+max(abs(a[1]),abs(b[1]))/camera.source_detector_mm)
        x0=max(0,int(math.floor((min(a[0],b[0])-rx)/camera.detector_width_mm*matrix+matrix/2-1)))
        x1=min(matrix,int(math.ceil((max(a[0],b[0])+rx)/camera.detector_width_mm*matrix+matrix/2+1)))
        y0=max(0,int(math.floor(matrix/2-(max(a[1],b[1])+ry)/camera.detector_height_mm*matrix-1)))
        y1=min(matrix,int(math.ceil(matrix/2-(min(a[1],b[1])-ry)/camera.detector_height_mm*matrix+1)))
        if x0>=x1 or y0>=y1:
            continue
        px,py=np.meshgrid(np.arange(x0,x1),np.arange(y0,y1))
        u=((px[...,None]+.5+offsets[:,0])/matrix-.5)*camera.detector_width_mm
        v=(.5-(py[...,None]+.5+offsets[:,1])/matrix)*camera.detector_height_mm
        rays=(camera.source_detector_mm*camera.forward + u[...,None]*camera.horizontal
              + v[...,None]*camera.vertical).reshape(-1,3)
        rays/=np.linalg.norm(rays,axis=1)[:,None]
        tmax=camera.source_detector_mm/(rays@camera.forward)
        ids=(((py*matrix+px)[...,None]*4)+np.arange(4)).ravel()
        for lo,hi,valid in capsule_intervals(camera.source,rays,start,end,radius,tmax):
            valid &= hi>lo
            if valid.any():
                groups.append(ids[valid]);starts.append(lo[valid]);ends.append(hi[valid])
    if not groups:
        return np.ones((matrix,matrix),dtype=np.float32)
    ids=np.concatenate(groups);lo=np.concatenate(starts);hi=np.concatenate(ends)
    order=np.lexsort((lo,ids));ids,lo,hi=ids[order],lo[order],hi[order]
    unique,group=np.unique(ids,return_inverse=True)
    # Positive offsets reset a vectorized running maximum between sample groups.
    span=float(hi.max()+1.)
    shifted=hi+group*span
    running=np.maximum.accumulate(shifted)
    prior=np.r_[-np.inf,running[:-1]]-group*span
    fresh=np.r_[True,ids[1:]!=ids[:-1]]
    prior[fresh]=-np.inf
    lengths=np.bincount(group,weights=np.maximum(hi-np.maximum(lo,prior),0),minlength=len(unique))
    sample_tau=mu*lengths
    if depth_mu>0:
        sample_tau+=depth_mu*lo[fresh]
    tau.ravel()[unique]=sample_tau
    return np.exp(-tau).mean(axis=2)
