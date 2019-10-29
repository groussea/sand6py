
#The output files HDF5 are saved with the hdf5_Write function and may be read with the hdf5_Read function

import h5py


def hdf5_Write(filename,coordinates,variables,chunks=True,compression="gzip", compression_opts=4):
    with h5py.File(filename, 'w') as f:
        g1=f.create_group('coordinates')
        g1.attrs['description']='Coordinates'
        for c in coordinates:             
            g1.create_dataset(c[0], data=c[1],chunks=chunks,compression=compression, compression_opts=compression_opts)        
        g2=f.create_group('variables')
        g2.attrs['description']='matricies'
        for v in variables:             
            g2.create_dataset(v[0], data=v[1],chunks=chunks,compression=compression, compression_opts=compression_opts)
        print('[log] Data saved in '+filename)

def hdf5_Read(filename):
    with h5py.File(filename, 'r') as f:
        coordinates=[]
        for vk in f['coordinates'].keys():
            coordinates.append([vk,f['coordinates'][vk][:]])
        variables=[]
        for vk in f['variables'].keys():
            variables.append([vk,f['variables'][vk][:]])
       
    return coordinates,variables


        
def hdf5_WriteUnstructured2DTimeserie(filename,Time,PointsSelected,VelSelected,chunks=True,compression="gzip", compression_opts=4):

    with h5py.File(filename, 'w') as f:
        f.create_dataset('time',data=Time,chunks=chunks,compression=compression, compression_opts=compression_opts)
        g=f.create_group('velocities') 
        for t,v in zip(Time,VelSelected):                     
            g.create_dataset(str(t), data=v[:,0:2],chunks=chunks,compression=compression, compression_opts=compression_opts)       
        g=f.create_group('coordinates')  
        for t,p in zip(Time,PointsSelected):                     
            g.create_dataset(str(t), data=p[:,0:2],chunks=chunks,compression=compression, compression_opts=compression_opts)            
    print('[log] Data saved in '+filename)        
 
def hdf5_ReadUnstructured2DTimeserie(filename):
    with h5py.File(filename, 'r') as f:
        Time=f['time'][:]
        grp=f['velocities']
        VelSelected=[i[:] for i in grp.values()]
        grp=f['coordinates']
        PointsSelected=[i[:] for i in grp.values()]

    return Time,PointsSelected,VelSelected




