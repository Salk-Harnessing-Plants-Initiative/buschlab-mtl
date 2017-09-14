#!/usr/bin/python
import sys
import os
import getpass
import shutil


# general options
species = "arabidopsis"  # one of "arabidopsis", "lotus", "medicago"

# mtl options
mtl_snpsdb = "/lustre/scratch/datasets/genotypes_for_pygwas/1.0.0/1/all_chromosomes_binary.hdf5"
mtl_traitfile2 = "/lustre/scratch/projects/busch_gwas/christian.goeschl/mtl-strigo-20170425/GWASinput_2016_control_means-na.csv"
mtl_filecolumns1 = "[63]" #"range(3,97)" #"range(1, 217)" #[1,2]
mtl_traitprefix2 = "c"
mtl_traitfile1 = "/lustre/scratch/projects/busch_gwas/christian.goeschl/mtl-strigo/GWASinput_2016_strigolactone_means-na.csv"
mtl_filecolumns2 = "[63]" #"range(3,97)" #"range(1, 217)" #[1,2]
mtl_traitprefix1 = "s" #+rnr
mtl_mac_threshold = 5
mtl_fileseparator = ','

# genehunter options
gh_upstream_dist = 4000
gh_downstream_dist = 4000
gh_pval_threshold = 1.0e-6
gh_mac_threshold = 10
gh_fdr_alpha = 0.05

###              no changes below this line                ###

# staging in
# this is only 2 files, so we do it on the login node
# since it is on cluster scratch, we assume we can savely remove any existing file with the same name
mtlworkdir = os.path.join("/lustre/scratch/projects/busch_gwas", getpass.getuser(), "_vs_".join([mtl_traitprefix1, mtl_traitprefix2]))
mtldestdir = os.path.join(mtlworkdir, "mtl_results")
if os.path.exists(mtlworkdir):
    shutil.rmtree(mtlworkdir)

os.makedirs(mtldestdir, mode=0o660)
shutil.copy2(mtl_traitfile1, mtlworkdir)
if mtl_traitfile1 != mtl_traitfile2:
    shutil.copy2(mtl_traitfile2, mtlworkdir)

# create mtl job
walltime = "12:00:00"
mem = "48gb"
mtl_opts_dict = {'MTL_FILE1': os.path.join(mtlworkdir, os.path.basename(mtl_traitfile1)),
                 'MTL_COLS1': mtl_filecolumns1,
                 'MTL_FILE2': os.path.join(mtlworkdir, os.path.basename(mtl_traitfile2)),
                 'MTL_COLS2': mtl_filecolumns2,
                 'MTL_PREFIX1': mtl_traitprefix1,
                 'MTL_PREFIX2': mtl_traitprefix2,
                 'MTL_SNPS': mtl_snpsdb,
                 'MTL_MAC': mtl_mac_threshold,
                 'MTL_OUTDIR': mtldestdir,
                 'MTL_FILE_SEPARATOR': mtl_fileseparator}

njobs = len(eval(mtl_filecolumns1))
mtlcmd = '/net/gmi.oeaw.ac.at/software/shared/busch_gwas/mtl/0.2.0/bin/mtl_process.sh'
mtlopts = ",".join(["=".join([key, "\"{}\"".format(val)]) for key, val in mtl_opts_dict.items()])

import subprocess

mtlqsubcmd = ['qsub']
if njobs > 1:
    mtlqsubcmd.append('-J0-{}'.format(njobs - 1))
mtlqsubcmd.append('-v{}'.format(mtlopts))
mtlqsubcmd.append('-P busch_gwas')
mtlqsubcmd.append('-l walltime={}'.format(walltime))
mtlqsubcmd.append('-l mem={}'.format(mem))
mtlqsubcmd.append('{}'.format(mtlcmd))

mtlproc = subprocess.Popen(mtlqsubcmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
mtlcomm = mtlproc.communicate()

#terminate unsuccessful if we were not able to create the mtl processing job.
if mtlcomm[0] == "":
    exit(2)


# create genehunter job

def select_db(sp):
    return {
        'arabidopsis:': os.environ['GENEHUNTER_AT_DB'],
        'lotus': os.environ['GENEHUNTER_LJ_DB'],
        'medicago': os.environ['GENEHUNTER_MT_DB']
    }[sp]

gh_opts_dict = {'GHUNTER_DATABASE': select_db(species),
                'GHUNTER_WORKDIR': mtldestdir,
                'GHUNTER_UDISTANCE': gh_upstream_dist,
                'GHUNTER_DDISTANCE': gh_downstream_dist,
                'GHUNTER_PVALTHRES': gh_pval_threshold,
                'GHUNTER_MACTHRES': gh_mac_threshold,
                'GHUNTER_FDR': gh_fdr_alpha,
                'GHUNTER_OUTPUTPREFIX': "_vs_".join([mtl_traitprefix1, mtl_traitprefix2])}

mtlcmd = '/net/gmi.oeaw.ac.at/software/shared/busch_gwas/mtl/0.2.0/bin/mtl_process.sh'
mtlopts = ",".join(["=".join([key, "\"{}\"".format(val)]) for key, val in mtl_opts_dict.items()])

ghqsubcmd = ['qsub']
ghqsubcmd.append('-Wdepend=afterok:{}'.format(mtlcomm[0]))

