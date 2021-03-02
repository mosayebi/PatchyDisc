# PatchyDisc
This is an implementation of the [Gaussian patchy disc](https://arxiv.org/abs/1607.06626) model. The code was built based on the [VMMC library](https://github.com/lohedges/vmmc) of Lester Hedges.

## Requirements
- cmake >= 3.11
- [nlohmann json library](https://github.com/nlohmann/json) 

## Compiling PatchyDisc
```
git clone git@github.com:mosayebi/PatchyDisc.git
cd PatchyDisc
mkdir build
cd build
cmake ..
make
```
A successfull compilation will create the executable file `build/PatchyDisc`.

## Running examples
```
cd examples
../build/PatchyDisc input_mixture.json
```
This will run a simulation using the parameters configureed in the `input_mixture.json`.

## Input file

Below is the contents of a sample input file for simulating a mixture of tri-valant and bi-valant patchy discs.

```json
{
    "initialisation" :
    {
        "mode": "from_init_conf",
        "from_random_conf":
        {
            "types": 2,
            "particle_numbers" : [{"type": 0, "N":300}, {"type": 1, "N":100}],
            "box": [40.0, 40.0]
        },
        "from_init_conf":
        {
            "init_conf": "last_conf.xyz",
            "restart_step_counter": false
        }
    },


    "topology" :
    {
        "patches": [
            {"type": 0, "nPatches": 3, "angles" : [0.0, 120, 240.0]},
            {"type": 1, "nPatches": 2, "angles" : [0.0, 120.0]}
        ],
        "pair_coeff" : [
            {"type1":0, "type2":0, "epsilon":15.0, "delta":1.0, "sigma":1.0, "sigma_p": 0.3, "rcut": 1.5},
            {"type1":0, "type2":1, "epsilon":5.0, "delta":1.0, "sigma":1.0, "sigma_p": 0.3, "rcut": 1.5},
            {"type1":1, "type2":0, "epsilon":5.0, "delta":1.0, "sigma":1.0, "sigma_p": 0.3, "rcut": 1.5},
            {"type1":1, "type2":1, "epsilon":15.0, "delta":1.0, "sigma":1.0, "sigma_p": 0.3, "rcut": 1.5}
        ]
    },

    "simulation parameters":
    {
        "seed": 42,
        "interaction": "GaussianPatchyDisc",
        "trajectory": "trajectory.xyz",
        "last_conf": "last_conf.xyz",
        "log_file" : "log.dat",
        "output_every": 1000,
        "sweeps": 1e5
    }
}
```
Currently the following interactions are implemented

- `GaussianPatchyDisc` -- It is described in the [reference](https://arxiv.org/abs/1607.06626).
- `GaussianPatchyDiscHR` -- Same as the `GaussianPatchyDisc`, but with a hard core repulsion when $r_{ij} \lt \sigma_{ij}$.
- `GaussianPatchyDiscHRSW` -- Same as the `GaussianPatchyDiscHRSW`, but with a Square Well poential when  $r_{ij} \le {r_{cut}}_{ij}$.

## Visualisation
The trajectory is saved in XYZ format. You can vissualise it in the visualaisation package of your choice. However, I recommend using [ovito](https://www.ovito.org/).

When importing the trajectory in ovito, please choose the follwing order for the column mappings: Particle type, Position X, Position Y, Orientation Z, Orientation W. You can then load particle shape STL files saved in the `UTILS/particle_shapes/` if you are simulating a standard patchy disc, or make your own particle shape using the python script provided.

To load a particle shape in ovito under the "Data source" panel click on "Particle Types" and choose a particle type. The "Load shape..."  bottom can then be found in the bottom of the panel.


![Binary mixture of patchy discs](mixture.gif)
