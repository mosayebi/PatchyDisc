# PatchyDisk
This is an implementation of the [Gaussian patchy disk](https://arxiv.org/abs/1607.06626) model. The code was built based on the [VMMC library](https://github.com/lohedges/vmmc) of Lester Hedges.

## Requirements
- cmake >= 3.11
- [nlohmann json library](https://github.com/nlohmann/json) 

## Compiling Patchy Disk
```
git clone git@github.com:mosayebi/PatchyDisk.git
cd PatchyDisk
mkdir build
cd build
cmake ..
make
```
After a successfull compilation, the executable can be found in `build/PatchyDisk`.

## Running examples
```
cd examples
../build/PatchyDisk input_mixture.json
```
This will run a simulation using the parameters configureed in the `input_mixtures.json`.

## Visualisation
The trajectory is saved in XYZ format. You can vissualise it in the visualaisation package of your choice. However, I recommend using [ovito](https://www.ovito.org/).

When importing the trajectory in ovito, please choose the follwing order for the column mappings: Particle type, Position X, Position Y, Orientation Z, Orientation W. You can then load particle shape STL files saved in the `UTILS/particle_shapes/` if you are simulating a standard patchy disk, or make your own particle shape using the python script provided.

To load a particle shape in ovito under the "Data source" panel click on "Particle Types" and choose a particle type. The "Load shape..."  bottom can then be found in the bottom of the panel.


![Binary mixture of patchy disks](mixture.gif)
