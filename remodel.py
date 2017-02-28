import subprocess
import os
from sys import platform as _platform


def remodel(rosetta_path, pdb_path, blueprint_path):
    if _platform == "linux" or _platform == "linux2":
        # linux
        print "Warning: Add binary name"
        pass
    elif _platform == "darwin":
        # MAC OS X
        binary_name = "remodel.macosclangrelease"
    elif _platform == "win32":
        # Windows
        print "Warning: Add binary name"
        pass
    binary_path = os.path.join(rosetta_path, "main", "source", "bin", binary_name)
    database_path = os.path.join(rosetta_path, "main", "database")

    command = "{} -database {} -s {} -remodel:blueprint {} -run:chain A  -remodel:quick_and_dirty".format(
        binary_path,
        database_path,
        pdb_path,
        blueprint_path)
    #command = "/Users/Naofumi/Downloads/rosetta/main/source/bin/remodel.macosclangrelease -database /Users/Naofumi/Downloads/rosetta/main/database/ -s 1c20.pdb -remodel:blueprint 1c20.remodel -run:chain A  -remodel:quick_and_dirty"
    subprocess.check_output(command,
                            stderr=subprocess.STDOUT,
                            shell=True)


if __name__ == '__main__':
    # config file just contain a absolute path to the root of rosetta
    rosetta_path = open("rosetta_path.config").read()
    pdb_path = os.path.join("pdb", "1c20.pdb")
    blueprint_path = os.path.join("blueprint", "1c20.remodel")
    remodel(rosetta_path, pdb_path, blueprint_path)