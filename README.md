# MD-molprep
## Short preparation for placement of small molecules in MD simulation

This will allow you to quickly build a "peptide".pdb out of small molecules that can be used in CharrmGUI to add to a membrame.

```1_rotater.py``` takes as input any .pdb that is found within ```inp_rotater``` and will make copies of it. The file should have the tag `N_LIGA.pdb`, where *N* is the number of copies you want, and *LIGA* is the ligand tag that it should have.\
If you run the script, it will take all the valid .pdb files, translate the center of mass (COM) to 0, randomly rotate the ligand, and store the multiples in `inp_placer`.

The next step is to input the files generated in the folder `upper`or `lower`, depending on in which leaflet you want the ligands to be placed. Then set `var_placer/params.txt` to your desired values (box size and z value for upper and lower leaflet).

At last, run `2_placer.py`. It takes all the .pdbs in `inp_placer`, as well as the parameters from the `params.txt` and combines them into one big .pdb file, aptly named `BigPDB.pdb`, found in `inp_placer`.\
It calculates the distribution of the COMs according to the Thomson problem in a box with repulsive walls, the results of which can be seen in `var_placer/coords/`. If coordinates were calculated once, they are stored and reused for further identical problems.


Install requirements with `pip install -r requirements.txt`.

Good Luck.

-V
