
##Generate the topology
gmx pdb2gmx -f pknot.pdb -o pknot_processed.gro -water spce -ignh
## 1: AMBER03 protein, nucleic AMBER94 (Duan et al., J. Comp. Chem. 24, 1999-2012, 2003)

gmx editconf -f pknot_processed.gro -o pknot_newbox.gro -c -d 2 -bt triclinic
gmx solvate -cp pknot_newbox.gro -cs spc216.gro -o pknot_solv.gro -p topol.top
gmx grompp -f ions.mdp -c pknot_solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o pknot_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
##Chose Grup  3:SOL
##Group     3 (            SOL) has 383157 elements

#gmx grompp -f minim.mdp -c pknot_solv_ions.gro -p topol.top -o em.tpr

#gmx mdrun -v -deffnm em

#gmx energy -f em.edr -o potential.xvg

##change in the mdp file 
##RNA Water_and_ions
#gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr

#gmx mdrun -deffnm nvt

#gmx energy -f nvt.edr -o temperature.xvg
## 18 0

#gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
#gmx mdrun -deffnm npt
#gmx energy -f npt.edr -o pressure.xvg
#gmx energy -f npt.edr -o density.xvg
## 24 0

#The Simulation, finally
#gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
#gmx mdrun -deffnm md_0_1

##Fix the trajectory 
#gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center

#####Analysis build-in calculations

#gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns
#gmx rms -s em.tpr -f md_0_1_noPBC.xtc -o rmsd_xtal.xvg -tu ns
#gmx gyrate -s md_0_1.tpr -f md_0_1_noPBC.xtc -o gyrate.xvg

#Select 1 and 1 for "RNA" RMSD for both.
