import ismrmrd as mrd
import numpy as np
import matplotlib.pyplot as plt

# load dataset
mrdName = r'D:\Willmering_Lab\Projects\General_Testing\20251003_Florida_GasRemoval_FLORETGX\XMRI_019-Test\XMRI_019\5)Xenon_3D_radial_Dixon\5)Xenon_3D_radial_Dixon_dixon.h5'
dset = mrd.Dataset(mrdName, "dataset", create_if_needed=False)
header = mrd.xsd.CreateFromDocument(dset.read_xml_header())

# get trajectories
traj_read = np.array([])
for acqnum in range(dset.number_of_acquisitions()-1):
    acq_temp = dset.read_acquisition(acqnum)
    if acqnum == 0:
        traj_read = acq_temp.traj[:]
        traj_read = traj_read[np.newaxis, :,:]
    else:
        traj_read = np.append(traj_read,acq_temp.traj[np.newaxis,:,:],axis=0)

plt.figure()
ax = plt.axes(projection='3d')
N_proj = np.size(traj_read,0)
N_visual = 100  # Number of projections you want to show, inefficient for large n
color = iter(plt.cm.viridis(np.linspace(0, 1, N_visual)))
for i in np.linspace(0, N_proj-1, N_visual):
    i = int(i)
    c = next(color)
    ax.scatter(traj_read[i, :, 0], traj_read[i, :, 1],
               traj_read[i, :, 2], color=c, s=0.5, marker='.')
ax.set_zlabel('$k_z$')
ax.set_ylabel('$k_y$')
ax.set_xlabel('$k_x$')
plt.title('Trajectory coordinates')


dset.close()

