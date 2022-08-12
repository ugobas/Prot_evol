# Prot_evol
Author: Ugo Bastolla Centro de Biologia Molecular Severo Ochoa (CSIC-UAM) ubastolla@cbm.csic.es

Citations:
1) Jimenez MJ, Arenas M, Bastolla U (2018) Substitution rates predicted by stability-constrained models of protein evolution are not consistent with empirical data. Mol Biol Evol. 35: 743-755
2) Arenas M, Weber CC, Liberles DA, Bastolla U (2017) ProtASR: An Evolutionary Framework for Ancestral Protein Reconstruction with Selection on Folding Stability. Syst Biol. 66:1054-1064. 
3) Arenas M, Sanchez-Cobos A, Bastolla U (2015) Maximum-Likelihood Phylogenetic Inference with Selection on Protein Folding Stability. Mol Biol Evol. 32:2195-207.
4) Minning J, Porto M, Bastolla U (2013) Detecting selection for negative design in proteins through an improved model of the misfolded state. Proteins 81:1102-12.

==================================
OVERVIEW
==================================

The program Prot_evol computes Structure and Stability Constrained Protein Evoluion (SSCPE) models in the form of independent site-specific amino acid substitution models (stationary frequencies and Halpern-Bruno exchangeability matrices) for each protein site. These substitution processes model selection for structural conservation (Str-CPE, structure constrained protein evolution) and thermodynamic stability of the native state against unfolding and misfolding (Stab-CPE).

The program implements eight selection models:
(1) The Stab-CPE mean-field (MF) in which stability is estimated in the mean-field produced by the other sites;
(2) The Stab-CPE wild-type (WT) model that models fitness as the change in folding free energy (DDG) of all possible mutations of the wild-type sequence.
(3,4) The Str-CPE (RMSD) and (DE) models that consider the predicted effect of each amino acid mutation of the wild-type sequence either on the mutant structure (RMSD) or on the free energy barrier between the mutant and the native structure (DE).
These predictions are produced by the program TNM (https://github.com/ugobas/tnm) that models each specific mutation as a perturbing force and predicts the deformation that it produces on the native structure as the linear response of the Torsional Network Model (TNM). These predictions are output in the files <>_mut_RMSD.dat and <>_mut_DE.dat, which must be input to Prot_evol (either both, just one or none).
(5-8) The structure and stability constrained model that combines the structure constrained DE or RMSD model with the stability constrained MF and WT model (DEMF, DEWT, RMSDMF or RMSDWT).

==================================
INPUT FILES
==================================

- Configuration file Input_Prot_evol.in included in the package that contains the input parameters that can be modified and their explanation.
- Protein data bank (PDB) file
- Multiple sequence alignment that includes the PDB sequence (optional)
- Structural effects of mutations predicted by the program TNM (optional)
- File "structures.in" included in the package with representative set of contact matrices used for pre-computing the misfolding model

The configuration file includes six sections:
A) Input files (PDB file, chain identifier)
B) Thermodynamic model. Modifying the configurational entropy with respect to
   unfolding SU1 and with respect to misfolding SC1 one type of stability will
   prevail on the other one.
C) Selection model
D) Evolution simulations. To activate them, set TMAX to a non-zero number.
E) Mutation model to build the background distribution and the 
   exchangeability matrices
F) Output control.
IMPORTANT: 
 - Modify the line PDB=... with the local path to your target PDB file
 - Modify the line FILE_STR=... with the local path to the file 
   structures.in, which is included in the package.

==================================
OUTPUT FILES
==================================

A) Site-specific substitution matrices are output for three models, selected through maximum likelihood of the input MSA: Stab-CPE (usually WT), Str-CPE model (usually RMSD) and SSCPE (usually RMSDWT, which is the one with best performances).
    A1) <PDB>_SSCPE_<model>_AA_profiles.txt 
    Site-specific amino acid frequencies and entropy (one site per line)

    A2) <PDB>_SSCPE_<model>_exchangeability_sites_<EXCH_MODEL>.txt
    Site-specific exchangeability matrices (one for each site)
    IMPORTANT: Due to their large size, site-specific exchangeability matrices
    are only printed if PRINT_E=1 in configuration file

    A3) <PDB>_SSCPE_<model>_AA_profile_global.txt
    Average amino acid frequencies across all sites
    
    A4) <PDB>_SSCPE_<model>_exchangeability_glob_<EXCH_MODEL>.txt
    Exchangeability matrix averaged over all sites.

    A5) <PDB>_SSCPE_<model>_rate_profile.dat 
    For each protein site we report:
    (1) the site identifier in the PDB,
    (2) the amino acid at the site in the PDB,
    (3) the secondary structure,
    (4) the number of native contacts,
    (5) the entropy of the predicted amino acid distribution,
    (6) the predicted average hydrophobicity.
    (7-9) the predicted substitution rate for the mut, exc and flux
    exchangeability models.
    Pearson correlation coefficients between these variables are 
    reported at the end of the file.
    
    A6) <PDB>_summary.dat
    Performances of all models and general information:
    PDB file, length, predicted folding free energy DG, hydrophobicity,
    thermodynamic parameters, input MSA (seq, length, entropy, mean identity
    with PDB), method for optimizing the selection parameter Lambda and for
    applied regularization.
    
    For each selection model, the following measures are given (for the SSCPE
    model, the global frequencies are determined twice and various
    regularization parameters are applied): (KL=Kullback-Leibler):
    - Model
    - Lambda0 Optimized selection parameter
    - Lambda1 Second selection parameter for combined model 
    - lik(MSA) Likelihood of site-specific MSA frequencies with respect
      to site-specific model (phylogenetically unaware)
    - KL(mod,reg.obs) KL divergence from model to regularized MSA frequencies
    - KL(reg.obs,mod) KL div. from regularized MSA frequencies to model
    - KL(mod,mut) KL div. from model to global backround frequencies
    - entropy(mod) Site-specific entropy of the model
    - DG Folding free energy of the model
    - Tf Freezing temperature of the misfolding model
    - h  Average hydrophobicity of the model

    A7) <PDB>_SSCPE_Cv.dat
    The same data as above are reported for every value of the regularization
    parameter reg (first column) for the optimal (usually RMSDWT) model,
    including the Cv measure dlik(MSA)/dreg used for selecting the
    regularization parameter and consequently the selection parameter Lambda

    A8) <PDB>_entropy.dat
    Site-specific entropy of the input MSA

B) Folding free energy:

    B1) <PDB>_DeltaG.dat
    Predicted folding free energy of the wild-type sequence

    B2) <ALI>.DeltaG
    Predicted folding free energy of each sequence in the MSA

    B2) <PDB>_Threading.dat
    Statistics of contacts used for computing the misfolding model

C) Simulations of evolution:

   Additionally, the program can also simulate stability constrained protein
   evolution, without applying the approximation that protein sites evolve
   independently of each other adopted to get the simulated matrices.

   C1) MSA of simulated sequences in the file <>_msa.fasta

   C2) Statistics of accepted mutations and substitutions that increase
   stability in <>_sel.dat

   C3) <PDB>_stab.dat
   For each amino acid substitution, it reports
   (1) the mutated site, (2) the mutated amino acid,
   (3) the native energy, (4) DG, (5) fitness, 
   (6) number of synonymous substitutions since the previous aa substitution
   (7) number of attempted mutations.
   This file can be used to reconstruct the protein sequences generated
   by the simulated evolution.

   C4) <PDB>_ave.dat
   Every it_print (internal parameter) substitutions, it prints average
   and standard error of: (1,2) fitness, (3,4) native energy, 
   (5,6) DeltaG, (7) difference between amino acid entropy of mutation and 
   selection, (8) non_synonymos/synonymous subst. rate, (9) acceptance rate.

   C4) <PDB>_final.dat
   Same information as in stab file, but printed at the end of the simulation.

==================================
DETAILED DESCRIPTION
==================================

1) SUBSTITUTION MATRICES
----------------------------
This computation is performed if MEANFIELD=1 in the configuration file.

Site-specific amino acid frequencies are modelled at each site i as

P^i(a,Lambda)=P^mut(a)*exp(Lambda*phi^i(a)/Z_i (1)

- phi^i(a) is a site-specific factor that represents log fitness (for instance, it is the predicted stability change DDG produced by the mutation to amino acid a at site i starting from the wild-type sequence in the WT model, or is the predicted RMSD produced by the same mutation in the RMSD model, or it is the self consistent effect of amino acid a at site i in the mean-field generated by the other sites in the MF model). These factors are predicted without any free parameter.
- Z_i=sum_a P^mut(a)*exp(Lambda*phi^i(a) is a normalization factor.
- P^mut(a) is a global amino acid frequency, determined in two steps in order to maximize the log-likelihood of the MSA at column i, MSA^i(a)=Number of sequences with amino acid a at column i
- Lambda is a global selection parameter that is optimized by one of two methods:

(1) Maximizing the log-likelihood of the MSA plus a regularization term analogous to ridge regression: LL(MSA^i)-REG*Lambda^2
(2) Maximizing the lok-likelihood of the MSA regularized with site-independent pseudo-counts and regularization parameter REG.
The regularization is necessary if MSA^i(a) is zero for some site and amino acid. The parameter REG is crucial since it determines the value of the selection parameter Lambda, which is the only free parameter of the SSCPE model.

There are three strategies to obtain REG and consequently Lambda.
1) REG is input through the configuration file Input_Prot_evol.in through OPT_REG=0 REG=something. Only 0<REG<=1 is allowed.
2) REG is optimized with OPT_REG=1 SCORE_CV=1 (reccomended). The chosen value of REG maximizes the "specific heat" (minus derivative of the likelihood of the observed site-specific amino acid distributions with respect to the regularization, see https://pubmed.ncbi.nlm.nih.gov/28555214/)
3) REG is optimized with OPT_REG=1 SCORE_CV=0. The chosen value of REG is such that the Kullback-Leibler divergence between the predicted distributions and the regularized observed distribution is the same in both directions.

   Global Frequencies

   The global distribution P^mut(a) can be computed in three different ways.

   1a) As the amino acid distribution observed across all protein 
   sites (command GET_FREQ=2, 19 degrees of freedom, Reccomended).
   If one of the amino acids is not present, this distribution is combined
   with the one at point (1b).
     
   1b) Fitted from a codon based mutation model, whose parameters 
   are the nucleotide frequencies pi_A, pi_T, pi_G, pi_C, the transition/
   tranversion ratio TT_RATIO and the enhancement of the mutation rate at
   CpG dinucleotides kCpG. These parameters are fitted imposing that the
   background distribution is as similar as possible to the observed one (1a)
   (GET_FREQ=1, 5 degrees of freedom). 
     
    1c) The parameters of the mutation model can be input manually 
    from the input file (GET_FREQ=0, 0 degrees of freedom).

    Global EXCHANGEABILITY

    The site-specific substitution process is given by
    (Q^i)_ab=(E^i)_ab*(P^i)_b

    where the EXCHANGEABILITY matrices(E^i)_ab are computed through the Halpern-Bruno model that relates the fixation probability with the stationary distributions, using the global exchangeability matrix (E^mut)_ab

    Four kinds of global exchangeability models are implemented.
    Three are based on an empirical substitution model that can be chosen with
    the command
    MATRIX=... (available options are WAG and JTT).

        2a) EXCHANGE=MUT The global exchangeability matrix is computed from
    the same mutation process that generates the background distribution.
        
        2b) EXCHANGE=EXCH The global exchangeability matrix is equal to the 
    empirical one (JTT or WAG, specified with MATRIX=WAG). The two former
    choices yield poor results.
    
        2c) EXCHANGE=RATE The exchangeability matrix is computed 
    imposing that the average amino acid substitution rate across 
    all sites is equal to the one generated by the empirical 
    rates, which is Q_ab=(E^emp)_ab*(f^emp)_b
        
        2d) EXCHANGE=FLUX (default) The exchangeability matrix is 
    computed imposing that the average flux of amino acids across all 
    sites is equal to the one generated by the empirical rates, which is 
    F_ab=(E^emp)_ab*(f^emp)_a*(f^emp)_b.

B) Simulations of protein evolution with selection on folding stability.

    Prot_evol simulates protein evolution subject to the constraint of selection
    on the folding stability of the native state against both unfolding and
    misfolding. 
    If only the mean-field distributions are of interest, this computation
    can be switched off with the command TMAX=0.


    At each time step, one mutated sequence is generated, its folding 
stability against unfolding and misfolding DG_mut is estimated with the 
model of Minning, Porto and Bastolla (Proteins 2013, 81:1102-112), and 
the mutation is fixed or rejected with one of three possible SELECTION 
models:

    3a) Neutral: fixation if DG_mut < DG_thr, rejected otherwise.
    
    DG_thr= 0.95*DG_PDB
    Commands NEUTRAL=1 MEANFIELD=0

    3b) Fixation probability of the Moran model with very low mutation 
    rate, P_fix=(1-f_wt/f_mut)/[1-(f_wt/f_mut)^N], where f_wt and f_mut 
    are the wild type and mutant fitness, respectively, 
    f=exp(-DG/T)/[1+exp(-DG/T)], and N is the effective population size,
    which is input with the command NPOP=... (default is 100).
    
    Commands NEUTRAL=0 MEANFIELD=0


==================================
COMPILE AND RUN:
==================================

To compile and run the program execute the following commands:

    unzip Prot_evol.zip
    make
    ./Prot_evol # Without parameter file, it will print an help screen
    ./Prot_evol Input_Prot_evol.in

==================================
INPUT FILE
==================================
    
    A detailled documentation of the input parameters needed will be found on 
Input_Prot_evol.in file included in this package.
It includes three sections:
A) Input files
   PDB= (PDB file)
   CHAIN= (chain identifier)
   ALI= File with MSA that includes the PDB sequence, optional
   STR_MUT= File with structural effects of mutations (mut_DE or mut_RMSD)
   computed by the program TNM, optional.
B) Substitution models
C) Amino acid frequencies
Mutation model to build the background distribution and the 
   exchangeability matrices
D) Exchangeability matrix
C) Selection model (recommended: do not modify)
E) Thermodynamic model
Modifying the configurational entropy with respect to
unfolding SU1 and with respect to misfolding SC1 one type of stability will
prevail on the other one.
F) Output control.
G) Evolution simulations. To activate them, set TMAX to a large number.






