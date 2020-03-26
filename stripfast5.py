#!/usr/bin/python
import h5py
import os
import shutil
import sys

def runalong(f,f2):

    for i,j in  f.attrs.items():
        f2.attrs[i]=j

    keys=f.keys()

    for i in keys:
        try:
            if isinstance(f[i],h5py.Group):
  #              if i not in ['Hairpin_Split_000','Alignment_000','Basecall_1D_000','Basecall_2D_000','Calibration_Strand_000']:
                    f2.create_group(i)
                    runalong(f[i],f2[i])
            elif isinstance(f[i],h5py.Dataset):
                name=f[i].name.split('/')[-1]
                shape=f[i].shape
                dtype=f[i].dtype

                if i not in ['Events','Signal']:
                  f2.create_dataset(i,shape,dtype,data=f[i].value)
        except: pass


f=h5py.File(sys.argv[1])
f2=h5py.File(sys.argv[1]+".strip",'w')


runalong(f,f2)
f2.close()

shutil.move(sys.argv[1]+".strip",sys.argv[1])
