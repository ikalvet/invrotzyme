# InvrotZyme

Script for building inverse rotamer assemblies out of a Rosetta matcher/enzdes constraint file.

This script will place sidechains according to the constraint file definitions, sample backbone positions, and optionally grow out extended backbone stubs (idealized helix or strand).
This script will perform an exhaustive analysis of all allowed rotamers and CST samplings.

You can also provide a motif PDB that will serve as a host for a particular constrained catalytic residue. That residue must exist in the PDB file, and only the rotamer will then be used for that residue.

The purpose of this tool is to find combinations of inverse rotamers that can be placed (on small extended backbones) without clashes. The outputs of this script can subsequently be used as inputs for RFdiffusion All-Atom to create protein backbones that host these active sites.



## Examples

A few usage examples are provided in `examples/`

**Kemp eliminase example:**
Places three catalytic residues around a benzisoxazole substrate. A HIS-GLU/ASP dyad on one side, and a SER/THR/TYR/GLN/ASN H-bond donor on the other side.
`cd examples/Kemp_eliminase ; python ../../invrotzyme.py --cstfile inputs/BIO_His_ED_oxy_nosample.cst --params inputs/BIO.params --dunbrack_prob 0.6 --frac_random_rotamers_per_cst 0.5 0.5 0.5 0.5 --secstruct_per_cst H H E --prefix outputs/ --suffix HHE`


**P450 example:**
Places a custom Heme ligand in complex with a substrate against a CYS-containing motif from a cytochrome P450 enzyme.
`cd examples/P450 ; python ../../invrotzyme.py --cstfile inputs/HBA_CYS_P450_nosample.cst --params inputs/HBA_unique.params --motif_for_cst 1:3:inputs/P450_motif.pdb --frac_random_rotamers 0.1 --prefix outputs/`


## Usage

First prepare a matcher/enzdes Constraint file according to the standard format outlined in Rosetta documentation:<br>
https://docs.rosettacommons.org/docs/latest/rosetta_basics/file_types/match-cstfile-format

This script requires all six degrees of freedom to be defined, so you msut provide distance, 2 angles, and 3 torsions for each interaction.

You can then run the script using many of the options below, perhaps taking inspiration from the provided examples.

```
options:
  -h, --help            show this help message and exit
  --cstfile CSTFILE     CST file used for matching. Keep sampling to minimum to avoid combinatorial explosion.
  --params PARAMS [PARAMS ...]
                        params files used by ligands and residues
  --keep_his_tautomer KEEP_HIS_TAUTOMER
                        Per cst, should a specific HIS tautomer (HIS or HIS_D) be used. Keeps only one the requested HIS tautomers. Format: 'cst_no:HIS/HIS_D,..'
  --dunbrack_prob DUNBRACK_PROB
                        Cumulative Dunbrack probability of used rotamers for any residue. As used by the -packing:dunbrack_prob_... flag in Rosetta.
  --dunbrack_prob_per_cst DUNBRACK_PROB_PER_CST [DUNBRACK_PROB_PER_CST ...]
                        Cumulative Dunbrack probability of used rotamers for each CST residue.
  --N_len N_LEN         Number of residues added to the stub N-term
  --C_len C_LEN         Number of residues added to the stub C-term
  --N_len_per_cst N_LEN_PER_CST [N_LEN_PER_CST ...]
                        Number of residues added to the stub N-term, per CST
  --C_len_per_cst C_LEN_PER_CST [C_LEN_PER_CST ...]
                        Number of residues added to the stub C-term, per CST
  --prune_ligand_rotamers PRUNE_LIGAND_ROTAMERS
                        Prunes the set of used ligand rotamers based on clashcheck, AND rmsd similarity cutoff.
  --max_random_rotamers MAX_RANDOM_ROTAMERS
                        Number of random rotamers picked for each residue for the sampling. Reasonable number would be below 20 for quick sampling.
  --max_random_rotamers_per_cst MAX_RANDOM_ROTAMERS_PER_CST [MAX_RANDOM_ROTAMERS_PER_CST ...]
                        Number of random rotamers picked for each CST block for the sampling. First value is for the ligand.
  --frac_random_rotamers FRAC_RANDOM_ROTAMERS
                        Fraction of rotamers that are randomly picked for each residue for the sampling.
  --frac_random_rotamers_per_cst FRAC_RANDOM_ROTAMERS_PER_CST [FRAC_RANDOM_ROTAMERS_PER_CST ...]
                        Fraction of rotamers that are randomly picked for each CST block for the sampling. First value is for the ligand.
  --secstruct SECSTRUCT
                        What secondary structure stub should be generated for each residue.
  --secstruct_per_cst SECSTRUCT_PER_CST [SECSTRUCT_PER_CST ...]
                        Per CST, what secondary structure stub should be generated for reaach residue.
  --motif_for_cst MOTIF_FOR_CST [MOTIF_FOR_CST ...]
                        Per CST, an external motif that should be used, instead of inverse rotamers. Only works for the first CST right now.
                        Format: cst_no:resno_in_motif:filepath ...
  --use_best_rotamer_cstids USE_BEST_ROTAMER_CSTIDS [USE_BEST_ROTAMER_CSTIDS ...]
                        CST ID's that should only use the best rotamer from each secondary structure bin. Numbering starts from 1.
  --extra_chi EXTRA_CHI
                        Enables extra CHI sampling on a given level for all CST's. Input format: chino:level,chino2:level2
  --extra_chi_per_cst EXTRA_CHI_PER_CST [EXTRA_CHI_PER_CST ...]
                        Enables extra CHI sampling on a given level for specific CST's. Input format: CSTNO1-chino:level,chino2:level2 CSTNO2-chino:level,chino2:level2
                        Sampling levels:
                          0 Default original dihedral only; same as using no flag at all
                          1 +/- one standard deviation (sd); 3 samples 
                          2 +/- 0.5 sd; 3 samples 
                          3 +/- 1 & 2 sd; 5 samples 
                          4 +/- 0.5 & 1 sd; 5 samples
                          5 +/- 0.5, 1, 1.5 & 2 sd; 9 samples 
                          6 +/- 0.33, 0.67, 1 sd; 7 samples 
                          7 +/- 0.25, 0.5, 0.75, 1, 1.25 & 1.5 sd; 13 samples.
  --suffix SUFFIX       Suffix to be added to the end of output PDB files
  --prefix PREFIX       Prefix to be added to the beginning of output PDB files
  --tip_atom            Inverse rotamers will be pre-selected based on whether the tip atoms are placed geometrically differently. Rotamer diversity is ignored.
  --debug               Debug mode. Printing more stuff out and running single-threaded
```

The script runs by default on multiple CPU cores using python multiprocessing. When submitted as a Slurm job, it will adjust the number of cores based on the environment variable `SLURM_CPUS_ON_NODE`.


### Best practices

Keep conformational sampling levels in the CST file to a minimum to avoid combinatorial explosion. Only sample torsions that are expectes to lead different valid assemblies.<br>

It's possible to limit the sampling by randomly picking rotamers for each residue, and limiting how the sidechain placements are sampled in the CST file.<br>
It's possible to control the length of the generated idealized backbone stub (from zero to ...).<br>
It's possible control most of the parameters separately for each constraint block.<br>
With using the `--tip_atom` argument it is possible to skip the inverse rotamer clash analysis, and only output assemblies based on their unique placement of catalytic atoms.

The output PDB files of this script will also contain the `REMARK 666 ...` lines which are required by the Rosetta enzdes constraint parser. As such, the outputs are suitable for building more complex enzyme design pipelines.<br>
For example, the published all-atom diffusion pipeline (https://github.com/ikalvet/heme_binder_diffusion) is directly compatible with the outputs of this script.


### Requirements

Python packages that are required:
```
pyrosetta
numpy
pandas
scipy
```
