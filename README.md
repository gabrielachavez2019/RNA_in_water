# RNA_in_water
RNA Molecular Dynamics Simulation is a tutorial to simulate an RNA pseudoknot in water and how changes on temperature and ions affect its structure

![](https://img.shields.io/twitter/url/http/@gabvet.svg?label=%40gabvet&style=social)

This folder contains all the files you need to start MD on a viral pseudoknot using the structure of: pknot.pdb 

This example will guide newbies through the process of setting up a simulation system containing an RNA molecule (pseudoknot) in a box of water, with ions. Each step will contain an explanation of input and output, using typical settings for general use.

This tutorial assumes you are using a [GROMACS](http://www.gromacs.org/) version 2018.7

## Requeriments

- Install ThinLinc Client
- Connect to earth.auburn.edu
- Use your ID Auburn credentials
- Open a terminal
- Make a work directory ie `mkdir pknot_MD` 
- Fisrt time you might need to source gromacs  `source /automnt/gromacs/default/bin/GMXRC`
- Get the [pknot.pdb](pknot.pdb) file (ie using `scp` )
- Get all the mdp files cloning this directory using [git](https://github.com/joaks1/au-bootcamp-git-intro) or downloading it as a zip file.

## Create topology

Verified that all the necessary atoms are present. Always check your pdb file for MISSING entries, as these entries indicate either atoms or whole residues that are not present in the crystal structure. Terminal regions may be absent, and may not present a problem for dynamics. Incomplete internal sequences or nucleic residues that have missing atoms will cause pdb2gmx to fail. These missing atoms/residues must be modeled in using other software packages: rosetta, robetta or molprobity. 

Please note that pdb2gmx is not magic. It cannot generate topologies for arbitrary molecules, just the residues defined by the force field (in the *.rtp files - generally proteins, nucleic acids, and a very finite amount of cofactors, like NAD(H) and ATP).
 
The purpose of pdb2gmx is to generate three files:

- The topology for the molecule.
- A position restraint file.
- A post-processed structure file.

The topology (topol.top by default) contains all the information necessary to define the molecule within a simulation. This information includes nonbonded parameters (atom types and charges) as well as bonded parameters (bonds, angles, and dihedrals). We will take a more detailed look at the topology once it has been generated.

Execute pdb2gmx by issuing the following command

`
gmx pdb2gmx -f pknot.pdb -o pknot_processed.gro -water spce -ignh
`

The structure will be processed by pdb2gmx, and you will be prompted to choose a force field:

```
GROMACS:      gmx pdb2gmx, version 2018.1
Executable:   /usr/bin/gmx
Data prefix:  /usr
Working dir:  /home/chris/pknot_MD
Command line:
  gmx pdb2gmx -f pknot.pdb -o pknot_processed.gro -water spce -ignh
Select the Force Field:
From '/usr/share/gromacs/top':
 1: AMBER03 protein, nucleic AMBER94 (Duan et al., J. Comp. Chem. 24, 1999-2012, 2003)
 2: AMBER94 force field (Cornell et al., JACS 117, 5179-5197, 1995)
 3: AMBER96 protein, nucleic AMBER94 (Kollman et al., Acc. Chem. Res. 29, 461-469, 1996)
 4: AMBER99 protein, nucleic AMBER94 (Wang et al., J. Comp. Chem. 21, 1049-1074, 2000)
 5: AMBER99SB protein, nucleic AMBER94 (Hornak et al., Proteins 65, 712-725, 2006)
 6: AMBER99SB-ILDN protein, nucleic AMBER94 (Lindorff-Larsen et al., Proteins 78, 1950-58, 2010)
 7: AMBERGS force field (Garcia & Sanbonmatsu, PNAS 99, 2782-2787, 2002)
 8: CHARMM27 all-atom force field (CHARM22 plus CMAP for proteins)
 9: GROMOS96 43a1 force field
10: GROMOS96 43a2 force field (improved alkane dihedrals)
11: GROMOS96 45a3 force field (Schuler JCC 2001 22 1205)
12: GROMOS96 53a5 force field (JCC 2004 vol 25 pag 1656)
13: GROMOS96 53a6 force field (JCC 2004 vol 25 pag 1656)
14: GROMOS96 54a7 force field (Eur. Biophys. J. (2011), 40,, 843-856, DOI: 10.1007/s00249-011-0700-9)
15: OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)
```
The force field will contain the information that will be written to the topology. This the most important choice! You should always read thoroughly about each force field and decide which is most applicable to your situation. For this tutorial, we will use AMBER03 for nucleic acid, so type `1` at the command prompt, followed by 'Enter'.

Please note that we are using the option:

-ignh: Ignore H atoms in the PDB file; especially useful for NMR structures. Otherwise, if H atoms are present, they must be in the named exactly how the force fields in GROMACS expect them to be. Different conventions exist, so dealing with H atoms can occasionally be a headache! If you need to preserve the initial H coordinates, but renaming is required, then the Linux sed command is your friend.

You have now generated three new files: pknot_processed.gro, topol.top, and posre.itp. pknot_processed.gro is a GROMACS-formatted structure file that contains all the atoms defined within the force field. The topol.top file is the system topology (more on this in a minute). The posre.itp file contains information used to restrain the positions of heavy atoms (more on this later).

Let's look at what is in the output topology (topol.top). Inspect its contents using vi or nano, don't change anythig!. After several comment lines (preceded by ;), you will find the following:

`
#include "amber03.ff/forcefield.itp"
`
This line calls the parameters within the AMBER03 force field. It is at the beginning of the file, indicating that all subsequent parameters are derived from this force field. The next important line is [ moleculetype ], below which you will find

```
[ moleculetype ]
; Name            nrexcl
RNA_chain_A         3

```

The next section defines the [ atoms ] in the RNA. The information is presented as columns:

```
[ atoms ]
;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB
; residue   1 G   rtp RG5  q -0.3
     1         OH      1      G    O5'      1    -0.6223         16   ; qtot -0.6223
     2         HO      1      G    H5T      2     0.4295      1.008   ; qtot -0.1928
     3         CT      1      G    C5'      3     0.0558      12.01   ; qtot -0.137
     4         H1      1      G   H5'1      4     0.0679      1.008   ; qtot -0.0691

```

The interpretation of this information is as follows:
nr: Atom number
type: Atom type
resnr: nucleotide residue number
residue: The  nucleotide residue name
atom: Atom name
cgnr: Charge group number
Charge groups define units of integer charge; they aid in speeding up calculations
charge: Self-explanatory
The "qtot" descriptor is a running total of the charge on the molecule
mass: Also self-explanatory
typeB, chargeB, massB: Used for free energy perturbation (not discussed here)
Subsequent sections include [ bonds ], [ pairs ], [ angles ], and [ dihedrals ]. Some of these sections are self-explanatory (bonds, angles, and dihedrals). The parameters and function types associated with these sections are elaborated on in Chapter 5 of the GROMACS manual. Special 1-4 interactions are included under "pairs" (section 5.3.4 of the GROMACS manual).

The remainder of the file involves defining a few other useful/necessary topologies, starting with position restraints. The "posre.itp" file was generated by pdb2gmx; it defines a force constant used to keep atoms in place during equilibration (more on this later).

```
; Include Position restraint file
#ifdef POSRES
#include "posre.itp"
#endif

; Include water topology
#include "amber03.ff/spce.itp"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif

; Include topology for ions
#include "amber03.ff/ions.itp"

[ system ]
; Name
pknot_RNA

[ molecules ]
; Compound        #mols
RNA_chain_A         1

```
This ends the moleculetype definition. The remainder of the topology file is dedicated to defining other molecules and providing system-level descriptions. The next moleculetype (by default) is the solvent, in this case SPC/E water. Other typical choices for water include SPC, TIP3P, and TIP4P. We chose this by passing "-water spce" to pdb2gmx. For an excellent summary of the many different water models, click here, but be aware that not all of these models are present within GROMACS.

The only parameter we are going to change at this stage is the [system] name. Gromacs is more commonly used for proteins so the defauld name is Protein, should be substitute for Zika_RNA for example.

A few key notes about the [ molecules ] directive:

The order of the listed molecules must exactly match the order of the molecules in the coordinate (in this case, .gro) file.
The names listed must match the [ moleculetype ] name for each species, not residue names or anything else.
If you fail to satisfy these concrete requirements at any time, you will get fatal errors from grompp (discussed later) about mismatched names, molecules not being found, or a number of others.

Now that we have examined the contents of a topology file, we can continue building our system.

## Generate a box and add water (solvatation)

Now that you are familiar with the contents of the GROMACS topology, it is time to continue building our system. In this example, we are going to be simulating a simple aqueous system. It is possible to simulate proteins and other molecules in different solvents, provided that good parameters are available for all species involved.

There are two steps to defining the box and filling it with solvent:
Define the box dimensions using the editconf module.
Fill the box with water using the solvate module (formerly called genbox).
You are now presented with a choice as to how to treat the unit cell. For the purpose of this tutorial, we will use a simple triclinic box as the unit cell. As you become more comfortable with periodic boundary conditions and box types, I highly recommend the rhombic dodecahedron, as its volume is ~71% of the cubic box of the same periodic distance, thus saving on the number of water molecules that need to be added to solvate the RNA molecule.

Let's define the box using editconf:

```
gmx editconf -f pknot_processed.gro -o pknot_newbox.gro -c -d 2.5 -bt triclinic
```

The above command centers the RNA fragment in the box (-c), and places it at least 2.5 nm from the box edge (-d 2.5). The box type is defined as a cube (-bt triclinic). The distance to the edge of the box is an important parameter. Since we will be using periodic boundary conditions, we must satisfy the minimum image convention. That is, an RNA molecule should never see its periodic image, otherwise the forces calculated will be spurious. Specifying a solute-box distance of 2.5 nm will mean that there are at least 2.0 nm between any two periodic images of the molecule. This distance will be sufficient for just about any cutoff scheme commonly used in simulations.

![](Zika_box02.png)

Now that we have defined a box, we can fill it with solvent (water). Solvation is accomplished using solvate:

```
gmx solvate -cp pknot_newbox.gro -cs spc216.gro -o pknot_solv.gro -p topol.top
```

The configuration of the RNA molecule (-cp) is contained in the output of the previous editconf step, and the configuration of the solvent (-cs) is part of the standard GROMACS installation. We are using spc216.gro, which is a generic equilibrated 3-point solvent model. You can use spc216.gro as the solvent configuration for SPC, SPC/E, or TIP3P water, since they are all three-point water models. The output is called pknot_solv.gro, and we tell solvate the name of the topology file (topol.top) so it can be modified. Note the changes to the [ molecules ] directive of topol.top:

```
[ molecules ]
; Compound        #mols
RNA_chain_A         1
SOL             127719
```

What solvate has done is keep track of how many water molecules it has added, which it then writes to your topology to reflect the changes that have been made. Note that if you use any other (non-water) solvent, solvate will not make these changes to your topology! Its compatibility with updating water molecules is hard-coded.

## Add ions

Since life does not exist at a net charge, we must add ions to our system.
We have a solvated system that contains a charged RNA. The output of pdb2gmx told us that the RNA has a net charge of +8e (based on its atomic composition). If you missed this information in the pdb2gmx output, look at the last line of your [ atoms ] directive in topol.top; it should read (in part) "qtot 8." 

The tool for adding ions within GROMACS is called genion. What genion does is read through the topology and replace water molecules with the ions that the user specifies. The input is called a run input file, which has an extension of .tpr; this file is produced by the GROMACS grompp module (GROMACS pre-processor), which will also be used later when we run our first simulation. What grompp does is process the coordinate file and topology (which describes the molecules) to generate an atomic-level input (.tpr). The .tpr file contains all the parameters for all of the atoms in the system.

To produce a .tpr file with grompp, we will need an additional input file, with the extension .mdp (molecular dynamics parameter file); grompp will assemble the parameters specified in the .mdp file with the coordinates and topology information to generate a .tpr file.

An .mdp file is normally used to run energy minimization or an MD simulation, but in this case is simply used to generate an atomic description of the system.

In reality, the .mdp file used at this step can contain any legitimate combination of parameters. I typically use an energy-minimization script, because they are very basic and do not involve any complicated parameter combinations. Please note that the files provided with this tutorial are intended only for use with the AMBER03 force field. Settings, particularly nonbonded interaction settings, will be different for other force fields.

Assemble your .tpr file with the following:

`gmx grompp -f ions.mdp -c pknot_solv.gro -p topol.top -o ions.tpr`


Now we have an atomic-level description of our system in the binary file ions.tpr. We will pass this file to genion:

`gmx genion -s ions.tpr -o pknot_solv_ions.gro -p topol.top -pname NA -nname CL -neutral`

Choose group 3 "SOL" for embedding ions. You do not want to replace parts of your RNA with ions.

In the genion command, we provide the structure/state file (-s) as input, generate a .gro file as output (-o), process the topology (-p) to reflect the removal of water molecules and addition of ions, define positive and negative ion names (-pname and -nname, respectively), and tell genion to add only the ions necessary to neutralize the net charge on the molecule by adding the correct number of negative ions (-neutral, which in this case will add 8 Cl- ions to offset the +8 charge on the molecule). You can also use genion to add a specified concentration of ions in addition to simply neutralizing the system by specifying the -neutral and -conc options in conjunction. Refer to the genion man page for information on how to use these options.

The names of the ions specified with -pname and -nname are always the elemental symbol in all capital letters, which is the [ moleculetype ] name that is then written to the topology. Residue or atom names may or may not append the sign of the charge (+/-), depending on the force field. Do not use atom or residue names in the genion command, or you will encounter errors in subsequent steps.

Visualize your molecule, water and ions with VMD and you should be able to see this:

![](Zika_water_ions.png)

## Energy Minimization

The solvated, electroneutral system is now assembled. Before we can begin dynamics, we must ensure that the system has no steric clashes or inappropriate geometry. The structure is relaxed through a process called energy minimization (EM).

The process for EM is much like the addition of ions. We are once again going to use grompp to assemble the structure, topology, and simulation parameters into a binary input file (.tpr), but this time, instead of passing the .tpr to genion, we will run the energy minimization through the GROMACS MD engine, mdrun.

Assemble the binary input using grompp using **em.tpr** input parameter file with nsteps = 5000:

`gmx grompp -f minim.mdp -c pknot_solv_ions.gro -p topol.top -o em.tpr`

Make sure you have been updating your topol.top file when running genbox and genion, or else you will get lots of nasty error messages ("number of coordinates in coordinate file does not match topology," etc).

We are now ready to invoke mdrun to carry out the EM:

`gmx mdrun -v -deffnm em`

The -v flag is for the impatient: it makes mdrun verbose, such that it prints its progress to the screen at every step. The -deffnm flag will define the file names of the input and output. So, if you did not name your grompp output "em.tpr," you will have to explicitly specify its name with the mdrun -s flag. In our case, we will get the following files:

- em.log: ASCII-text log file of the EM process
- em.edr: Binary energy file
- em.trr: Binary full-precision trajectory
- em.gro: Energy-minimized structure

There are two very important factors to evaluate to determine if EM was successful. The first is the potential energy (printed at the end of the EM process, even without -v). Epot should be negative, and (for a simple protein or RNA in water) on the order of 105-106, depending on the system size and number of water molecules. The second important feature is the maximum force, Fmax, the target for which was set in minim.mdp - "emtol = 1000.0" - indicating a target Fmax of no greater than 1000 kJ mol-1 nm-1. It is possible to arrive at a reasonable Epot with Fmax > emtol. If this happens, your system may not be stable enough for simulation. Evaluate why it may be happening, and perhaps change your minimization parameters (integrator, emstep, etc).

Let's do a bit of analysis. The em.edr file contains all of the energy terms that GROMACS collects during EM. You can analyze any .edr file using the GROMACS energy module:

`gmx energy -f em.edr -o potential.xvg`

At the prompt, type "10 0" to select Potential (10); zero (0) terminates input. You will be shown the average of Epot, and a file called "potential.xvg" will be written. To plot this data, you will need the Grace plotting tool. The resulting plot should look something like this, demonstrating the nice, steady convergence of E_pot:

`xmgrace potential.xvg`

![](Energy_Minimization_Zika.png)

Now that our system is at an energy minimum, we can begin real dynamics.

## Equilibration in Temperature

EM ensured that we have a reasonable starting structure, in terms of geometry and solvent orientation. To begin real dynamics, we must equilibrate the solvent and ions around the RNA molecule. If we were to attempt unrestrained dynamics at this point, the system may collapse. The reason is that the solvent is mostly optimized within itself, and not necessarily with the solute. It needs to be brought to the temperature we wish to simulate and establish the proper orientation about the solute (the RNA). After we arrive at the correct temperature (based on kinetic energies), we will apply pressure to the system until it reaches the proper density.

Remember that **posre.itp** file that pdb2gmx generated a long time ago? We're going to use it now! The purpose of **posre.itp** is to apply a position restraining force on the heavy atoms of the molecule (anything that is not a hydrogen). Movement is permitted, but only after overcoming a substantial energy penalty. The utility of position restraints is that they allow us to equilibrate our solvent around our RNA molecule, without the added variable of structural changes in the RNA. The origin of the position restraints (the coordinates at which the restraint potential is zero) is provided via a coordinate file passed to the -r option of grompp.

Equilibration is often conducted in two phases. The first phase is conducted under an NVT ensemble (constant Number of particles, Volume, and Temperature). This ensemble is also referred to as "isothermal-isochoric" or "canonical." The timeframe for such a procedure is dependent upon the contents of the system, but in NVT, the temperature of the system should reach a plateau at the desired value. If the temperature has not yet stabilized, additional time will be required. Typically, 50-100 ps should suffice, and we will conduct a 100-ps NVT equilibration for this exercise. Depending on your machine, this may take a while (just under an hour if run in parallel on 16 cores or so).

We will call grompp and mdrun just as we did at the EM step:

```
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt
```

A full explanation of the parameters used can be found in the [GROMACS](http://www.gromacs.org) manual, in addition to the comments provided. Take note of a few parameters in the .mdp file:

- gen_vel = yes: Initiates velocity generation. Using different random seeds (gen_seed) gives different initial velocities, and thus multiple (different) simulations can be conducted from the same starting structure.
- tcoupl = V-rescale: The velocity rescaling thermostat is an improvement upon the Berendsen weak coupling method, which did not reproduce a correct kinetic ensemble.
- pcoupl = no: Pressure coupling is not applied.

Let's analyze the temperature progression, again using energy:

`gmx energy -f nvt.edr -o temperature.xvg`

Type "16 0" at the prompt to select the temperature of the system and exit. The resulting plot should look something like the following after `xmgrace temperature.xvg`

![](Temperature_Zika.png)

From the plot, it is clear that the temperature of the system quickly reaches the target value (300 K), and remains stable over the remainder of the equilibration. For this system, a shorter equilibration period (on the order of 50 ps) may have been adequate.

## Equilibration in Pressure

The previous step, *NVT* equilibration, stabilized the temperature of the system. Prior to data collection, we must also stabilize the pressure (and thus also the density) of the system. Equilibration of pressure is conducted under an NPT ensemble, wherein the Number of particles, Pressure, and Temperature are all constant. The ensemble is also called the "isothermal-isobaric" ensemble, and most closely resembles experimental conditions.

The npt.mdp it is not drastically different from the parameter file used for NVT equilibration. Note the addition of the pressure coupling section, using the Parrinello-Rahman barostat.

A few other changes:

continuation = yes: We are continuing the simulation from the NVT equilibration phase
gen_vel = no: Velocities are read from the trajectory (see below)
We will call grompp and mdrun just as we did for NVT equilibration. Note that we are now including the -t flag to include the checkpoint file from the NVT equilibration; this file contains all the necessary state variables to continue our simulation. To conserve the velocities produced during NVT, we must include this file. The coordinate file (-c) is the final output of the NVT simulation.

```
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt
```

Let's analyze the pressure progression, again using energy:

`gmx energy -f npt.edr -o pressure.xvg`

Type "18 0" at the prompt to select the pressure of the system and exit. The resulting plot should look something like the following:

![](Pressure_Zika.png)

The pressure value fluctuates widely over the course of the 100-ps equilibration phase, but this behavior is not unexpected. The running average of these data are plotted as the red line in the plot. Over the course of the equilibration, the average value of the pressure is 7.5 ± 160.5 bar. Note that the reference pressure was set to 1 bar, so is this outcome acceptable? Pressure is a quantity that fluctuates widely over the course of an MD simulation, as is clear from the large root-mean-square fluctuation (160.5 bar), so statistically speaking, one cannot distinguish a difference between the obtained average (7.5 ± 160.5 bar) and the target/reference value (1 bar).

Let's take a look at density as well, this time using energy and entering "24 0" at the prompt.

`gmx energy -f npt.edr -o density.xvg`

As with the pressure, the running average of the density is also plotted in red. The average value over the course of 100 ps is 1019 ± 3 kg m-3, close to the experimental value of 1000 kg m-3 and the expected density of the SPC/E model of 1008 kg m-3. The parameters for the SPC/E water model closely replicate experimental values for water. The density values are very stable over time, indicating that the system is well-equilibrated now with respect to pressure and density.

if your density values do not match. Pressure-related terms are slow to converge, and thus you may have to run NPT equilibration slightly longer.

## Molecular Dynamics

Upon completion of the two equilibration phases, the system is now well-equilibrated at the desired temperature and pressure. We are now ready to release the position restraints and run production MD for data collection. The process is just like we have seen before, as we will make use of the checkpoint file (which in this case now contains preserve pressure coupling information) to grompp. We will run a 1-ns MD simulation, the script for which can be found here.

`gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr`

grompp will print an estimate for PME load, which will dictate how many processors should be dedicated to the PME calculation, and how many for the PP calculations. Refer to the GROMACS 4 publication and the manual for details.

Estimate for the relative computational load of the PME mesh part: 0.22
For a cubic box, the optimal setup will have a PME load of 0.25 (3:1 PP:PME - we're very close to optimal!); for a dodecahedral box, the optimal PME load is 0.33 (2:1 PP:PME). When executing mdrun, the program should automatically determine the best number of processors to assign for the PP and PME calculations. Thus, make sure you indicate an appropriate number of threads/cores for your calculation (the value of -nt X), so that you can get the best performance.

Now, execute mdrun:

`gmx mdrun -deffnm md_0_1`

In GROMACS 2018, the PME calculations can be offloaded to graphical processing units (GPU), which speeds up the simulation substantially.

## Analysis

Now that we have simulated our RNA, we should run some analysis on the system. What types of data are important? This is an important question to ask before running the simulation, so you should have some ideas about the types of data you will want to collect in your own systems. For this tutorial, a few basic tools will be introduced.

The first is trjconv, which is used as a post-processing tool to strip out coordinates, correct for periodicity, or manually alter the trajectory (time units, frame frequency, etc). For this exercise, we will use trjconv to account for any periodicity in the system. The RNA molecule will diffuse through the unit cell, and may appear "broken" or may "jump" across to the other side of the box. To account for such actions, issue the following:

`gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center`

Select 1 ("RNA") as the group to be centered and 0 ("System") for output. We will conduct all our analyses on this "corrected" trajectory. Let's look at structural stability first. 

### RMSD calculations
GROMACS has a built-in utility for RMSD calculations called rms. To use rms, issue this command:

`gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns`

Choose 1 ("RNA") for both the least-squares fit and the group for RMSD calculation. The -tu flag will output the results in terms of ns, even though the trajectory was written in ps. This is done for clarity of the output (especially if you have a long simulation - 1e+05 ps does not look as nice as 100 ns). The output plot will show the RMSD relative to the structure present in the minimized, equilibrated system:

If we wish to calculate RMSD relative to the crystal structure, we could issue the following:

`gmx rms -s em.tpr -f md_0_1_noPBC.xtc -o rmsd_xtal.xvg -tu ns`

Results look something like:

![](RMDS_100ns.png)

### Radius of gyration
The radius of gyration of an RNA molecule is a measure of its compactness. If the RNA is stably folded, it will likely maintain a relatively steady value of Rg. If the RNA unfolds, its Rg will change over time. Let's analyze the radius of gyration for Zika in our simulation:

`gmx gyrate -s md_0_1.tpr -f md_0_1_noPBC.xtc -o gyrate.xvg`

Choose group 1 (RNA) for analysis.
![](Gyration.png)

We can see from the reasonably invariant Rg values that the RNA remains very unestable, does not remain folded in the same way ovet the time form over the course of 1 ns at 300 K. This result is not unexpected, but illustrates an advanced capacity of GROMACS analysis that comes built-in.

Happy simulating!

![](Zika_Toomer.png)

# Acknowledgments

## Material
This exercise was inspired by, and borrowed heavily from, the Lysozyme in Water [Gromacs tutorial](http://www.mdtutorials.com/gmx/lysozyme/index.html) can be found at <http://www.mdtutorials.com/gmx/lysozyme/index.html>.

## Support
This work was made possible by funding provided to [Joanna Sztuba-Solinska](https://jsztuba.wixsite.com/sztubasolinska)


# License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/deed.en_US"><img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/deed.en_US">Creative Commons Attribution 4.0 International License</a>.


# PS: For adding GROMACS to your PATH

```
cd /automnt/gromacs/default/bin/
echo "" >> ~/.bashrc 
echo "export PATH=\"\$PATH:$(pwd)\"" >> ~/.bashrc 
source ~/.bashrc 
cd
```


