import numpy as np
from scipy import ndimage as ndi
from skimage.morphology import ball
import sample2_reconstruction as S
HU = np.load("grid1.npy"); sp = np.array([0.9765625,0.9765625,1.0]); vox = float(np.prod(sp))
d = np.load("/mnt/user-data/uploads/sample1.npz"); origin = np.array([d["x"].min(), d["y"].min(), d["z"].min()], float)
seg = np.load("seg1.npy"); present = HU >= 500
ceram = ndi.binary_dilation(seg > 0, structure=ball(2))
rest = present & ~ceram
print("stored volume %.1f cm3; ceramic pieces (dilated) %.1f cm3; the rest %.1f cm3" % (present.sum()*vox/1000, (present & ceram).sum()*vox/1000, rest.sum()*vox/1000))
opened = ndi.binary_opening(rest, structure=ball(3))
lab, n = ndi.label(opened, structure=np.ones((3,3,3))); sizes = np.bincount(lab.ravel())[1:]
order = np.argsort(sizes)[::-1]
print("components of the rest after opening: %s cm3" % [round(sizes[i]*vox/1000,1) for i in order[:5]])
core = lab == (order[0]+1)
rec = core.copy()
for _ in range(3): rec = ndi.binary_dilation(rec, structure=ball(1)) & rest
print("brown piece: %.1f cm3 after recovering its edge" % (rec.sum()*vox/1000))
# remove protruding debris: opening with a larger ball, keeping the largest component
clean = ndi.binary_opening(rec, structure=ball(5))
lab2, n2 = ndi.label(clean, structure=np.ones((3,3,3))); s2 = np.bincount(lab2.ravel())[1:]
clean = lab2 == (np.argmax(s2)+1)
clean = ndi.binary_fill_holes(clean)
print("cleaned (lumps under about 10 mm removed): %.1f cm3; debris removed %.1f cm3" % (clean.sum()*vox/1000, (rec.sum()-clean.sum())*vox/1000))
np.save("brown_raw.npy", rec); np.save("brown_clean.npy", clean)
X = np.argwhere(clean)*sp + origin; c = X.mean(0)
U, Sg, Vt = np.linalg.svd(X-c, full_matrices=False); sd = Sg/np.sqrt(len(X))
print("cleaned piece PCA std %s mm (two equal values + one distinct = solid of revolution); extents along PCs %s mm" % (np.round(sd,1), np.round(((X-c)@Vt.T).max(0)-((X-c)@Vt.T).min(0),0)))
print("PC directions:", np.round(Vt,3).tolist())
mean_hu = HU[clean & present].astype(float).mean() if False else None
