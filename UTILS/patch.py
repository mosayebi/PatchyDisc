import numpy as np
import os
import sys
import json

class Reader:
    def __enter__(self):
        return self

    def __exit__(self,  exc_type, exc_val, exc_tb):
        self.__del__()

    def __del__(self):
        try:
            if self._conf:
                self._conf.close()
        except:
            print("ERROR: The reader could not load the provided configuration")
            sys.exit(1)

    def __init__(self, configuration):
        if not os.path.isfile(configuration):
            print("Configuration file '{}' is not readable".format(configuration))
            sys.exit(1)
        self._conf = open(configuration, "r")
        self._time = None
        self._box = np.zeros(2)
        self._len = 0

    """
    Special reader that handles the first line and allocates the memory for storing the system
    """

    def read(self, n_skip=0):
        for i in range(n_skip):
            self._line = self._conf.readline().split()
            if  len(self._line) == 0:
                return False
            else:
                self._N = int(self._line[0])
                for j in range(self._N+1):
                    self._line = self._conf.readline().split()

        self._line = self._conf.readline().split()
        if len(self._line) == 0:
            return False
        self._N = int(self._line[0])
        self._line = self._conf.readline().split()
        self._time, self._box = np.float(self._line[0]), np.array(self._line[1:3], dtype=float)

        self._positions = np.zeros((self._N, 2), dtype=float)
        self._orientations = np.zeros((self._N, 2), dtype=float)
        self._type = np.zeros(self._N, dtype=int)

        for i in range(self._N):
            self._line = self._conf.readline().split()
            if (len(self._line) == 0):
                print("ERROR: Reader encountered a partial configuration with only {} lines ({} expected)".format(i, self._N))
            for j in range(2):
                self._positions[i, j] = float( self._line[1+j] )
                self._orientations[i, j] = float( self._line[3+j] )
        self._configuration = {'types': self._type,
                               'positions': self._positions,
                               'orientations':self._orientations,
                               'box': self._box,
                               'time': self._time
        }
        return self._configuration

def read_top(input_json):
    with open(input_json) as f:
        inp = json.load(f)

    nTypes = inp['types']
    Ni = inp['particle_numbers']

    pair_coeff = {}
    items = ['epsilon', 'delta', 'sigma', 'sigma_p', 'rcut']
    for item in items:
        pair_coeff[item] = np.zeros((nTypes, nTypes), dtype=float)
    patchAngles, ts, counter = [], [], -1
    for i in range(nTypes):
         ts.append(inp['patches'][i]['type'])
         patchAngles.append(inp['patches'][i]['angles'])
         for j in range(nTypes):
             counter += 1
             t1 = inp['pair_coeff'][counter]['type1']
             t2 = inp['pair_coeff'][counter]['type2']
             for item in items:
                 pair_coeff[item][t1, t2] = inp['pair_coeff'][counter][item]
    patchAngles = [x for y, x in sorted(zip(ts, patchAngles))] #sort based on type
    return patchAngles, pair_coeff

def _geenrate_stl(filename, sigma=1, angles=[0,120,240]):
    import cadquery as cq
    r0 = (cq.Workplane("front")
        .moveTo(sigma/7, sigma/10)
        .lineTo(sigma/2 - sigma/9, sigma/10)
        .lineTo(sigma/2 - sigma/10, sigma/9)
        .threePointArc((sigma/2, 0), (sigma/2 - sigma/10, -sigma/9))
        .lineTo(sigma/2 - sigma/9, -sigma/10)
        .lineTo(sigma/7, -sigma/10)
        .threePointArc((-sigma/7, 0), (sigma/7, sigma/10))
        .close()
        .extrude(sigma/15).edges().fillet(0.01)
    )
    for i, a in enumerate(angles):
        r0_copy = r0.rotate((0,0,0),(0,0,1), a)
        if i==0 :
            result = r0_copy
        else:
            result = result.union(r0_copy)
    cq.exporters.export(result, filename)


def generate_monomer_stl_file(input_json):
    patchAngles, pair_coeff = read_top(input_json)
    for t, angles in enumerate(patchAngles):
        sigma = pair_coeff['sigma'][t][t]
        _geenrate_stl(f'monomer_{t}.stl', sigma=sigma, angles=angles)


conf = '../examples/trajectory.xyz'
system = Reader(conf)
conf = system.read()
top = read_top('../examples/input.json')
print(top[0])
# while conf:
#     print(system.read()['time'])
#     conf = system.read()