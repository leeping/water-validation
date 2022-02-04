#!/usr/bin/env python

import os
import shutil
from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *
from forcebalance.nifty import *
from forcebalance.nifty import _exec
# from forcebalance.output import getLogger, RawStreamHandler, INFO
import warnings
warnings.simplefilter('ignore')
from forcebalance.molecule import *
from forcebalance.gmxio import *
import argparse
from collections import OrderedDict

parser = argparse.ArgumentParser()
parser.add_argument('--nve', action='store_true', help='Launch NVE simulations.')
parser.add_argument('-p', type=int, default=8241, help='Specify Work Queue port.')
parser.add_argument('-r', action='store_true', help='Read NVE results instead of launching.')

args, sys.argv = parser.parse_known_args(sys.argv)
calc = (not args.r)

if calc:
    wq_port = args.p
    work_queue.set_debug_flag('all')
    wq = work_queue.WorkQueue(port=wq_port, exclusive=False, shutdown=False)
    wq.specify_name('mdrun')
    wq.specify_keepalive_interval(8640000)
    print('Work Queue listening on %d' % (wq.port))

toptmp = """#include "{ffitp}"

[ system ]
{nmol} water molecules

[ molecules ]
SOL {nmol}
"""

sizes = [216, 343, 512, 1000, 1600, 2500, 4000] #, 16000]
# sizes = [216]
# temps = [260, 273, 285, 298.15, 310, 323]
temps = [298.15]
press = [1.0]

# Spacing of 10 means 1000 simulations
# Spacing of 100 means 100 simulations
spac = {216:10, 343:10, 512:20, 1000:50, 1600:50, 2500: 100, 4000:100, 16000:200}
for ffd in ["TIP3P", "SPC-E", "TIP4P", "TIP4P-Ew", "TIP4P2005", "TIP3P-FB", "TIP4P-FB"]:
    ffd = os.path.join(ffd, "Gmx-Diffuse")
    if len([i for i in os.listdir(ffd) if i.endswith('.xml')]) != 1:
        raise Exception('Please ensure only one .xml file in this folder')
    ffxml = [i for i in os.listdir(ffd) if i.endswith('.xml')][0]
    ffitp = ffxml[:-4]+'.itp'
    if not os.path.exists(os.path.join(ffd,ffitp)):
        raise Exception('Please ensure %s is here' % ffitp)

    for sz in sizes:
        if not calc: continue
        szd = os.path.join(ffd, "%i" % sz)
        printcool("Working in " + szd)
        # Create the GROMACS coordinate file if it does not exist.
        if not os.path.exists(os.path.join(szd, "conf.gro")):
            # Create an OpenMM PDB object.
            pdb = PDBFile(os.path.join(szd, "water.pdb"))
            # Use OpenMM functions to add virtual particles to water molecules.
            m = Modeller(pdb.topology, pdb.positions)
            ff = ForceField(os.path.join(ffd, ffxml))
            m.addExtraParticles(ff)
            system = ff.createSystem(m.topology, nonbondedMethod=PME, nonbondedCutoff=0.9*nanometer)
            integ = VerletIntegrator(0.002*picosecond)
            print "OpenMM minimizing energy.."
            simul = Simulation(pdb.topology, system, integ)
            simul.context.setPositions(m.positions)
            simul.minimizeEnergy()
            print "done"
            # Write the OpenMM output to a temporary PDB file.
            with open(os.path.join(szd, ".temp.pdb"),'w') as f: PDBFile.writeModel(m.topology, simul.context.getState(getPositions=True).getPositions(), f)
            # Load both the new and the original PDB files using the Molecule class.
            M = Molecule(os.path.join(szd, ".temp.pdb"))
            M0 = Molecule(os.path.join(szd, "water.pdb"))
            # Copy over the periodic box, and correct the residue and atom names.
            M.boxes = M0.boxes
            M.resname = ['SOL' for i in range(M.na)]
            # This is the bit that needs to be changed for 3pt vs 4pt
            if 'TIP3' in ffd:
                M.atomname = [['OW','HW1','HW2'][i%3] for i in range(M.na)]
            elif 'SPC' in ffd:
                M.atomname = [['OW','HW1','HW2'][i%3] for i in range(M.na)]
            elif 'TIP4' in ffd:
                M.atomname = [['OW','HW1','HW2','MW'][i%4] for i in range(M.na)]
            else:
                raise Exception('Spoo!')
            # Write the GROMACS coordinate file.
            M.write(os.path.join(szd, "conf.gro"))
    
        # Create the GROMACS topology file if it does not exist.
        if not os.path.exists(os.path.join(szd, "topol.top")):
            # Write the GROMACS topology file.
            fout = open(os.path.join(szd, "topol.top"), "w")
            print >> fout, toptmp.format(nmol=sz, ffitp=ffitp)
            fout.close()
    
        # Copy the force field and the mdp files.
        if not os.path.exists(os.path.join(szd, ffitp)) : os.symlink("../%s" % ffitp, os.path.join(szd, ffitp))
    
        # Call GROMACS to minimize the energy.
        if not os.path.exists(os.path.join(szd, "min.gro")) : 
            _exec("grompp -maxwarn 1 -f min.mdp -c conf.gro -o min.tpr", print_to_screen=True, expand_cr=True, cwd=szd)
            _exec("mdrun -v -deffnm min", print_to_screen=True, expand_cr=True, cwd=szd)
    
        for temp in temps:
            for pres in press:
                ptd = os.path.join(szd, "%.2fK-%.1fatm" % (temp, pres))
                if not os.path.exists(ptd): os.makedirs(ptd)
                eqdict = OrderedDict([("ref_t", "%.2f" % temp), ("gen_temp", "%.2f" % temp), ("ref_p", "%.1f" % pres)]) 
                mddict = OrderedDict([("ref_t", "%.2f" % temp), ("ref_p", "%.1f" % pres)])
                if args.nve:
                    if not os.path.exists(os.path.join(ptd, 'md.edr')):
                        if not os.path.exists(os.path.join(ptd, 'result.tar.bz2')):
                            print "result.tar.bz2 does not exist in", ptd
                            continue
                        else:
                            print "extracting from result.tar.bz2 in", ptd
                            os.system("tar xjf %s --directory %s 2> /dev/null" % (os.path.join(ptd,"result.tar.bz2"),ptd))
                            # os.remove(os.path.join(ptd,"result.tar.bz2"))
                    else:
                        print "results already extracted in", ptd
                    nums = [i+spac[sz] for i in range(0, 10000, spac[sz])]
                    do_extract = True
                    if os.path.exists(os.path.join(ptd, "NVE")):
                        if len(nums) == int(_exec("find %s -name conf.gro | wc -l" % os.path.join(ptd, "NVE"))[0].strip()):
                            do_extract = False
                    if do_extract:
                        _exec("trjconv -f %s/md.trr -s %s/md.tpr -o %s/out.gro -sep -b %i -skip %i -nzero 5"% (ptd, ptd, ptd, spac[sz], spac[sz]/10), stdin="System")
                        for t0, conf in zip(nums, sorted([i for i in os.listdir('%s' % ptd) if (i.startswith('out') and i.endswith('.gro'))])):
                            nved = os.path.join(ptd, "NVE", "%05i" % t0)
                            if not os.path.exists(nved): os.makedirs(nved)
                            shutil.move(os.path.join(ptd, conf), os.path.join(nved, "conf.gro"))
                    for t0 in nums:
                        nved = os.path.join(ptd, "NVE", "%05i" % t0)
                        if not os.path.exists(os.path.join(nved, 'msd.xvg')):
                            input_files = [(os.path.join(szd, "%s" % i), "%s" % i) for i in "topol.top", ffitp, "nve.mdp"] + [("gmxnve.sh", "gmxnve.sh")]
                            input_files.append((os.path.join(nved, "conf.gro"), "conf.gro"))
                            output_files = [(os.path.join(nved, "%s" % i), "%s" % i) for i in "nve.log", "msd.xvg", "energy.xvg"]
                            queue_up_src_dest(wq, "sh gmxnve.sh &> nve.log", input_files=input_files, output_files=output_files)
                else:
                    write_mdp(os.path.join(ptd, "eq.mdp"), eqdict, fin="%s/eq.mdp" % szd)
                    write_mdp(os.path.join(ptd, "md.mdp"), mddict, fin="%s/md.mdp" % szd)
                    if not os.path.exists(os.path.join(ptd, "min.gro")) : os.symlink("../min.gro", os.path.join(ptd, "min.gro"))
                    if not os.path.exists(os.path.join(ptd, "topol.top")) : os.symlink("../topol.top", os.path.join(ptd, "topol.top"))
                    if not os.path.exists(os.path.join(ptd, ffitp)) : os.symlink("../../%s" % ffitp, os.path.join(ptd, ffitp))
                    if not os.path.exists(os.path.join(ptd, "md.trr")) :
                        input_files = [(os.path.join(ptd, "%s" % i), "%s" % i) for i in "min.gro", "topol.top", ffitp, "eq.mdp", "md.mdp"] + [("rungmx.sh", "rungmx.sh")]
                        output_files = [(os.path.join(ptd, "%s" % i), "%s" % i) for i in "result.tar.bz2", "rungmx.log"]
                        queue_up_src_dest(wq, "sh rungmx.sh &> rungmx.log", input_files=input_files, output_files=output_files)
        if calc: wq_wait1(wq, wait_time=10)

    for temp in temps:
        if not args.r: continue
        for pres in press:
            print "%.2fK-%.1fatm" % (temp, pres)
            for sz in sizes:
                szd = os.path.join(ffd, "%i" % sz)
                ptd = os.path.join(szd, "%.2fK-%.1fatm" % (temp, pres))
                diffuse = []
                for t0 in [i+spac[sz] for i in range(0, 10000, spac[sz])]:
                    nved = os.path.join(ptd, "NVE", "%05i" % t0)
                    if os.path.exists(os.path.join(nved, 'msd.xvg')):
                        diffuse.append(float(os.popen("grep 'D\[' %s" % os.path.join(nved, 'msd.xvg')).readlines()[0].split()[4]))
                diffuse = np.array(diffuse)
                L = (sz / AVOGADRO_CONSTANT_NA * (18.015 * gram/mole) / (995 * (kilogram / meter**3))).value_in_unit(angstrom**3)**(1./3)
                print "Mol = %5i" % sz, "N = %5i" % len(diffuse), "1/L = %.3f" % (1.0/L), "D = % .3f" % np.mean(diffuse), "+-", "%.3f" % (np.std(diffuse)/np.sqrt(len(diffuse)))
   
if calc: wq_wait(wq)
