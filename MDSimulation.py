import sys, os
import getopt
import math

# Objects of type Atom store their coordinates, velocities and forces as well as their connectivity information
class Atom:
    def __init__(self, x, y, z, vx, vy, vz, bonded_atoms=None):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.bonded_atoms = bonded_atoms
        self.non_bonded_atoms = []
        # Atoms automatically initialize with zero forces acting on them
        self.fx = 0.0
        self.fy = 0.0
        self.fz = 0.0

    def distance(self, other):
        dx = self.x - other.x
        dy = self.y - other.y
        dz = self.z - other.z
        return math.sqrt(dx**2 + dy**2 + dz**2)
   
    
def toRVCString(coord):
    s = str(round(coord,4))
    if '-' in s:
        return "  "+s
    else: return "   "+s


def simulateMD(inFilename, outFilename, kB=40000.0, kN=400.0, nbCutoff=0.50, dt=0.0010, mass=12.0):
    fIn = open(inFilename)
    atoms = [Atom(0, 0, 0, 0, 0, 0)]  # So indexes on list line up with atom numbers
    # Parse the input file and store all the atoms as Atom objects in a list
    headerfound = 0
    for line in fIn:
        if (line[0] == '#') & (headerfound == 0):
            headerfound = 1
        if line.strip()=='' or line[0]=='#': continue
        data = line.split()
        x, y, z = map(float, data[1:4])
        vx, vy, vz = map(float, data[4:7])
        bonded_atoms = [[i,0] for i in list(map(int, data[7:]))]
        atom = Atom(x, y, z, vx, vy, vz, bonded_atoms)
        atoms.append(atom)
    fIn.close()
    # Calculate and store reference bond and nonbond distances
    for i in range(len(atoms)):
        if i == 0: continue
        for j in atoms[i].bonded_atoms:
            j[1] = atoms[i].distance(atoms[j[0]])
        for j in range(i+1,len(atoms)):
            # If the atoms are bonded there are no nonbond interactions
            isBonded = False
            for k in atoms[i].bonded_atoms:
                if k[0] == j:
                    isBonded = True
            if (not isBonded) & (atoms[i].distance(atoms[j]) < nbCutoff):
                atoms[i].non_bonded_atoms.append([j , atoms[i].distance(atoms[j])])
                atoms[j].non_bonded_atoms.append([i , atoms[i].distance(atoms[j])])
    # Run the simulation
    fOutEnergies = open(outFilename,'w')
    fOutRVC = open(outFilename.replace(".out", ".rvc"),'w')
    fOutEnergies.write("# step E_k E_b E_nb E_tot\n")
    for i in range(1,1001):
        eK = 0
        eB = 0
        eNb = 0
        # Update each atom using the first two steps of Verlet Integration
        for j in range(1,len(atoms)):
            atom = atoms[j]
            atom.vx = atom.vx+(atom.fx/mass)*(dt/2)
            atom.vy = atom.vy+(atom.fy/mass)*(dt/2)
            atom.vz = atom.vz+(atom.fz/mass)*(dt/2)
            atom.x = atom.x + atom.vx*dt
            atom.y = atom.y + atom.vy*dt
            atom.z = atom.z + atom.vz*dt
        # Calculate forces and store the potential energies
        for j in range(1,len(atoms)):
            atom = atoms[j]
            atom.fx = 0
            atom.fy = 0
            atom.fz = 0
            for other in atom.bonded_atoms:
                d = atom.distance(atoms[other[0]])
                d0 = other[1] # This is the reference distance 
                atom.fx += kB*(d-d0)*((atoms[other[0]].x-atom.x)/d)
                atom.fy += kB*(d-d0)*((atoms[other[0]].y-atom.y)/d)
                atom.fz += kB*(d-d0)*((atoms[other[0]].z-atom.z)/d)
                if other[0]>j: # Wont sum the energy if the bond has already been counted
                    eB += (kB*(d-d0)*(d-d0)) / 2
            for other in atom.non_bonded_atoms:
                d = atom.distance(atoms[other[0]])
                d0 = other[1] # This is the reference distance 
                atom.fx += kN*(d-d0)*((atoms[other[0]].x-atom.x)/d)
                atom.fy += kN*(d-d0)*((atoms[other[0]].y-atom.y)/d)
                atom.fz += kN*(d-d0)*((atoms[other[0]].z-atom.z)/d)
                if other[0]>j: # Wont sum the energy if the interaction has already been counted
                    eNb += (kN*(d-d0)*(d-d0)) / 2
        # Update Velocities again and use those values to compute kinetic energy
        for j in range(1,len(atoms)):
            atom = atoms[j]
            atom.vx = atom.vx+(atom.fx/mass)*(dt/2)
            atom.vy = atom.vy+(atom.fy/mass)*(dt/2)
            atom.vz = atom.vz+(atom.fz/mass)*(dt/2)
            eK += mass * (atom.vx * atom.vx + atom.vy * atom.vy + atom.vz * atom.vz)/2
        # Write to .out file and to .rvc file if appropriate
        eTotal = eK + eB + eNb
        fOutEnergies.write(str(i)+" "+str(round(eK,1))+" "+str(round(eB,1))+" "+str(round(eNb,1))+" "+str(round(eTotal,1))+"\n")
        if i % 10 == 0:
            fOutRVC.write("# At step "+str(i)+", energy = "+str(eTotal)+"kJ\n")
            for j in range(1,len(atoms)):
                atom = atoms[j]
                fOutRVC.write("    "+str(j)+toRVCString(atom.x)+toRVCString(atom.y)+toRVCString(atom.z)+toRVCString(atom.vx)+toRVCString(atom.vy)+toRVCString(atom.vz))
                for other in atom.bonded_atoms:
                    fOutRVC.write("     "+str(other[0]))
                fOutRVC.write("\n")
    fOutEnergies.close()
    fOutRVC.close()


def main():
    inFile = sys.argv[1]
    outFile = sys.argv[2]
    kB = float(sys.argv[3])
    kN = float(sys.argv[4])
    nbCutoff = float(sys.argv[5])
    dt = float(sys.argv[6])
    mass = float(sys.argv[7])
    simulateMD(inFile,outFile,kB,kN,nbCutoff,dt,mass)
    

if __name__=='__main__':
    main()