import subprocess
import os
from sys import platform as _platform


def relax(rosetta_path, pdb_path, outfolder='./rosetta_out'):
    if _platform == "linux" or _platform == "linux2":
        # linux
        binary_name = "relax.linuxgccrelease"
    elif _platform == "darwin":
        # MAC OS X
        binary_name = "relax.macosclangrelease"
    elif _platform == "win32":
        # Windows
        print "Warning: Add binary name"
        pass
    binary_path = os.path.join(rosetta_path, "main", "source", "bin", binary_name)
    database_path = os.path.join(rosetta_path, "main", "database")

    command = "{} -database {} -s {} -ignore_unrecognized_res -overwrite -nstruct 1 -relax:constrain_relax_to_start_coords -out:path:pdb {} -out:path:score {}".format(
        '\''+binary_path+'\'',
        '\''+database_path+'\'',
        '\''+pdb_path+'\'',
        '\''+outfolder+'\'',
        '\''+outfolder+'\'')
    print command
    subprocess.check_output(command,
                            stderr=subprocess.STDOUT,
                            shell=True)


def remodel(rosetta_path, pdb_path, blueprint_path, outfolder='./rosetta_out'):
    if _platform == "linux" or _platform == "linux2":
        # linux
        binary_name = "remodel.linuxgccrelease"
    elif _platform == "darwin":
        # MAC OS X
        binary_name = "remodel.macosclangrelease"
    elif _platform == "win32":
        # Windows
        print "Warning: Add binary name"
        pass
    binary_path = os.path.join(rosetta_path, "main", "source", "bin", binary_name)
    database_path = os.path.join(rosetta_path, "main", "database")

    command = "{} -database {} -s {} -remodel:blueprint {} -run:chain A -remodel:dr_cycles 3 -nstruct 1 -overwrite -ex1 -ex2 -ignore_zero_occupancy false -ignore_unrecognized_res true".format(
        '\''+binary_path+'\'',
        '\''+database_path+'\'',
        '\''+pdb_path+'\'',
        '\''+blueprint_path+'\'',
        '\''+outfolder+'\'',
        '\''+outfolder+'\'')
    print command
    #remodel:dr_cycles 3
    #ex1
    #ex2
    #nstruct 1
    #ignore_zero_occupancy false
    #-remodel:quick_and_dirty
    #command = "/Users/Naofumi/Downloads/rosetta/main/source/bin/remodel.macosclangrelease -database /Users/Naofumi/Downloads/rosetta/main/database/ -s 1c20.pdb -remodel:blueprint 1c20.remodel -run:chain A  -remodel:quick_and_dirty"
    subprocess.check_output(command,
                            stderr=subprocess.STDOUT,
                            shell=True)


if __name__ == '__main__':
    # config file just contain a absolute path to the root of rosetta
    rosetta_path = open("rosetta_path.config").read().strip('\n')
    pdb_path = os.path.join("pdb", "1c20.pdb")
    blueprint_path = os.path.join("blueprint", "1c20.remodel")
    remodel(rosetta_path, pdb_path, blueprint_path)
