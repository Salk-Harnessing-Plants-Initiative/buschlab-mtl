"""
    mtl (multi trait limix)
    author: christian goeschl
    date:   2016-09-16
"""

import sys
import time

import h5py
import limix.qtl as qtl
import numpy as np
import os
import pandas as pd
import scipy as sp

import pygwas_modules.kinship as kinship
import pygwas_modules.plotting as gplt
import pygwas_modules.result as res


class MTL:
    def __init__(self, mac_thres=0):
        self.mac_thres = mac_thres
        self.phenotypes = None
        # self.snps = None
        self.iid = None
        self.ibs = None
        self.ts_norm = None
        # self.bnorm_K = None
        self.used_snp_pos = None
        self.macs = None
        self.mafs = None
        self.chromosomes = None
        self.chr_names = None

    def read_phenotype_col(self, phenotype_filepath, colnr, colprefix="", sep='\t'):
        sys.stdout.write("reading phenotypes: {}, col: {}\n".format(phenotype_filepath,colnr))
        with open(phenotype_filepath, 'r') as phenofile:
            hcols = phenofile.readline().strip().split(sep)
            p = []

            for l in phenofile:
                dcols = l.strip().split(sep)
                try:
                    p.append([dcols[0], np.float64(dcols[colnr])])
                except ValueError:
                    sys.stdout.write("excluding accession {} because of trait value {}\n".format(dcols[0], dcols[colnr]))
                    continue

            data = pd.DataFrame(p)

            ids = data[0].values
            data.index = ids
            data = data[list(range(1, data.shape[1]))]
            data.columns = ["_".join([colprefix, hcols[colnr]])]

            if self.phenotypes is None:
                self.phenotypes = data
            else:
                pheno_acc_ids = list(set(self.phenotypes.index) & set(ids))
                self.phenotypes = pd.concat([self.phenotypes.loc[pheno_acc_ids], data.loc[pheno_acc_ids]], axis=1)
                sys.stdout.write("phenotype intersection is {} accessions.\n".format(len(pheno_acc_ids)))

        self.phenotypes.sort_index(axis=0, inplace=True)
        return

    def write_phenotypes(self, path):
        self.phenotypes.to_csv(path, sep=',', index=True, index_label="Acc_ID")

    def read_genotypes(self, genotype_filepath):
        sys.stdout.write("reading genotypes ... ")
        sys.stdout.flush()
        with h5py.File(genotype_filepath, 'r') as genofile:
            geno_acc_ids = list(genofile["/accessions"].value)
            pheno_geno_acc_intersect = list(set(geno_acc_ids) & set(self.phenotypes.index))
            geno_acc_idx = np.in1d(genofile['/accessions'], pheno_geno_acc_intersect)
            snps = genofile["/snps"][:, geno_acc_idx]
            geno_acc_ids = np.array(geno_acc_ids)[geno_acc_idx]

            chr_names = genofile['positions'].attrs.get('chrs')
            chr_regions = np.array(genofile['positions'].attrs.get('chr_regions'))
            geno_chroms = []
            for i, reg in enumerate(chr_regions):
                geno_chroms.extend(np.repeat(chr_names[i], reg[1]-reg[0]))
            pos = genofile['/positions'].value
        sys.stdout.write("ok.\n")

        macs = np.array(snps.sum(axis=1)).astype(int)
        macs_th = (macs >= self.mac_thres) & (macs <= snps.shape[0]-self.mac_thres)
        snps = snps[macs_th, :]
        sys.stdout.write("removed {:d} snps because of MAC threshold {:d}. (Remaining snps: {:d}.)\n"
                         .format(pos.shape[0]-snps.shape[0], self.mac_thres, snps.shape[0]))
        pos = pos[macs_th]
        geno_chroms = np.array(geno_chroms)[macs_th]

        snps = pd.DataFrame(snps, index=pd.MultiIndex.from_arrays([geno_chroms, pos]), columns=geno_acc_ids)
        snps = snps.reindex_axis(sorted(snps.columns), axis=1)

        accs_no_geno_info = np.array(self.phenotypes.index)[np.invert(np.in1d(self.phenotypes.index, pheno_geno_acc_intersect))]
        self.phenotypes.drop(accs_no_geno_info, inplace=True)
        sys.stdout.write("no genotype information for accessions: {}. Removed them from list of phenotypes.".format(accs_no_geno_info))

        self.ibs = np.array(kinship.calc_ibs_kinship(snps.values))

        # snps.index = pd.MultiIndex.from_arrays([geno_chroms, pos])
        #geno_acc_ids
        for ix, reg in enumerate(chr_regions):
            self.chromosomes[reg[0]:reg[1]] = self.chr_names[ix]

        self.iid = sorted(list(set(geno_acc_ids) & set(self.phenotypes.index)))
        del geno_acc_ids

        sys.stdout.write("genotype-phenotype intersection is {} accessions.\n".format(len(self.iid)))

        snps = np.array(snps.loc[self.iid])
        snpsshape = snps.shape


        ts = snps[:, macs_th]
        sys.stdout.write("creating kinship matrix ... ")
        sys.stdout.flush()
        start = time.time()
        self.ibs = kinship.calc_ibs_kinship(ts.T)

        # self.bnorm_K = kinship.scale_k(ibs).astype(np.float64)
        elapsed = time.time() - start
        sys.stdout.write("ok. ({} s)\n".format(elapsed))

        self.used_snp_pos = pos[macs_th]
        self.macs = macs[macs_th]
        self.mafs = self.macs/float(snpsshape[0])
        self.chromosomes = self.chromosomes[macs_th]
        # ts=sub_snps[:,(sumts<sub_snps.shape[0]*0.99)&(sumts>(sub_snps.shape[0]*0.01))]
        self.ts_norm = (ts - ts.mean(axis=0)) / ts.std(axis=0)
        return

    # def create_kinship(self):
    #     sys.stdout.write("creating kinship matrix ... ")
    #     sys.stdout.flush()
    #     sub_snps = self.snps.loc[self.iid]
    #     sub_snps = np.array(sub_snps)
    #     sub_snps = sub_snps.astype(np.float64)
    #
    #     # self.snps = self.snps.loc[self.iid]
    #     # sumts = self.snps.sum(axis=0)
    #     # ts = self.snps[:, (sumts != 0) & (sumts != self.snps.shape[0])]
    #     # self.used_snp_pos = self.snps.columns[(sumts != 0) & (sumts != self.snps.shape[0])].astype(np.float64)
    #     mac = sub_snps.sum(axis=0)
    #     maf = float(mac)/sub_snps.shape[0]
    #     ts = sub_snps[:, (mac != 0) & (mac != sub_snps.shape[0])]
    #     self.used_snp_pos = self.snps.columns[(mac != 0) & (mac != sub_snps.shape[0])].astype(np.float64)
    #     # ts=sub_snps[:,(sumts<sub_snps.shape[0]*0.99)&(sumts>(sub_snps.shape[0]*0.01))]
    #     self.ts_norm = (ts - ts.mean(axis=0)) / ts.std(axis=0)
    #
    #     start = time.time()
    #     self.ibs = kinship.calc_ibd_kinship(ts.T)
    #     elapsed = time.time() - start
    #     sys.stdout.write("ok. ({} s)\n".format(elapsed))

    def box_cox_transform(self, values, lambda_range=(-2.0, 2.0), lambda_increment=0.1, verbose=False, method='standard'):
        """
        Performs the Box-Cox transformation, over different ranges, picking the optimal one w. respect to normality.
        """
        from scipy import stats
        a = sp.array(values)
        if method == 'standard':
            vals = (a - min(a)) + 0.1 * sp.std(a)
        else:
            vals = a
        sw_pvals = []
        lambdas = sp.arange(lambda_range[0], lambda_range[1] + lambda_increment, lambda_increment)
        for l in lambdas:
            if l == 0:
                vs = sp.log(vals)
            else:
                vs = ((vals ** l) - 1) / l
            r = stats.shapiro(vs)
            if sp.isfinite(r[0]):
                pval = r[1]
            else:
                pval = 0.0
            sw_pvals.append(pval)
        # log.info(sw_pvals)
        i = sp.argmax(sw_pvals)
        l = lambdas[i]
        if l == 0:
            vs = sp.log(vals)
        else:
            vs = ((vals ** l) - 1) / l
        # self._perform_transform(vals,"box-cox")
        sys.stdout.write('optimal lambda was %0.1f\n' % l)
        return vals

    def do_qtl(self, outputdir, rnr=None):
        pheno_norm = self.phenotypes.loc[self.iid].values

        # p1 = np.array(self.phenotypes[[0]].loc[self.iid]).astype(np.float64)
        # p1 = (p1 - p1.min()) / (p1.max() - p1.min())
        # p2 = np.array(self.phenotypes[[1]].loc[self.iid]).astype(np.float64)
        # p2 = (p2 - p2.min()) / (p2.max() - p2.min())
        # pheno_norm = np.concatenate((p1, p2), axis=1)


        #exp transform (does not converge)

        #sqrt transform (does not converge)
        # p1 = np.array(self.phenotypes[[0]].loc[self.iid]).astype(np.float64)
        # p1 = np.sqrt((p1 - p1.min()) + 0.1 * np.std(p1))
        # p2 = np.array(self.phenotypes[[1]].loc[self.iid]).astype(np.float64)
        # p2 = np.sqrt((p2 - p2.min()) + 0.1 * np.std(p2))
        # pheno_norm = np.concatenate((p1, p2), axis=1)

        # pheno = np.array(self.phenotypes[[0, 1]].loc[self.iid])
        # pheno_norm = (pheno - pheno.mean(axis=0)) / pheno.std(axis=0)

        # box - cox - transform (does not converge)
        # p1 = np.array(self.phenotypes[[0]].loc[self.iid]).astype(np.float64)
        # p1 = self.box_cox_transform(p1)
        # p2 = np.array(self.phenotypes[[1]].loc[self.iid]).astype(np.float64)
        # p2 = self.box_cox_transform(p2)
        # pheno_norm = np.concatenate((p1, p2), axis=1)

        # ascombe transform (does not converge)
        # p1 = np.array(self.phenotypes[[0]].loc[self.iid]).astype(np.float64)
        # p1 = 2.0 * sp.sqrt(p1 + 3.0 / 8.0)
        # p2 = np.array(self.phenotypes[[1]].loc[self.iid]).astype(np.float64)
        # p2 = 2.0 * sp.sqrt(p2 + 3.0 / 8.0)
        # pheno_norm = np.concatenate((p1, p2), axis=1)



        # QTL
        n_pheno = pheno_norm.shape[1]  # number of traits
        # N = len(self.iid)  # number of accessions
        covs = None
        Acovs = None
        K1r = self.ibs
        covar_type = 'freeform'

        # Testing for GxE effect
        Asnps0 = sp.ones((1, n_pheno))  # common effects: degree of freedom is 1
        Asnps1 = sp.zeros((2, n_pheno))
        Asnps1[0, :] = 1.0
        Asnps1[1, 0] = 1.0

        sys.stdout.write("calulating qtl ... \n")
        sys.stdout.flush()
        start = time.time()
        pvalues_inter = qtl.qtl_test_interaction_lmm_kronecker(snps=self.ts_norm, phenos=pheno_norm, covs=covs, Acovs=Acovs,
                                                           Asnps0=Asnps0,
                                                           Asnps1=Asnps1, K1r=K1r)
        elapsed = time.time() - start
        sys.stdout.write("qtl finished. ({} s)\n".format(elapsed))

        if not os.path.isdir(outputdir):
            sys.stdout.write("creating output directory: {} ... ".format(outputdir))
            sys.stdout.flush()
            os.makedirs(outputdir)
            sys.stdout.write("ok.\n")
            sys.stdout.flush()

        sys.stdout.write("plotting and writing results ... \n")
        sys.stdout.flush()
        pvalues_inter = np.array(pvalues_inter)
        pvalues_inter = pvalues_inter[:, 0, :]

        if rnr is not None:
            fileprefix = "{}-mac{}-run{}".format("-x-".join(self.phenotypes.columns), self.mac_thres, rnr)
        else:
            fileprefix = "{}-mac{}".format("-x-".join(self.phenotypes.columns), self.mac_thres)

        # specific (G x E)
        sys.stdout.write("... writing specific interaction results ... ")
        start = time.time()
        gwas_result = res.GWASResult(self.chr_names, self.chromosomes,
                                     self.used_snp_pos, pvalues_inter[0],
                                     dict(mafs=self.mafs, macs=self.macs),
                                     additional_columns={})
        gwas_result.save_as_csv(os.path.join(outputdir, "{}_specific_pvals.csv".format(fileprefix)))
        gplt.plot_gwas_result(gwas_result, os.path.join(outputdir, "{}_specific_manhattan.png".format(fileprefix)),
                              mac=self.mac_thres)
        gplt.plot_qq(gwas_result, os.path.join(outputdir, "{}_specific_qq.png".format(fileprefix)))
        sys.stdout.write("ok ({:f} s)\n".format(time.time() - start))

        # common
        sys.stdout.write("... writing common interaction results ... ")
        start = time.time()
        gwas_result = res.GWASResult(self.chr_names, self.chromosomes,
                                     self.used_snp_pos, pvalues_inter[1],
                                     dict(mafs=self.mafs, macs=self.macs),
                                     additional_columns={})
        gwas_result.save_as_csv(os.path.join(outputdir, "{}_common_pvals.csv".format(fileprefix)))
        gplt.plot_gwas_result(gwas_result, os.path.join(outputdir, "{}_common_manhattan.png".format(fileprefix)),
                              mac=self.mac_thres)
        gplt.plot_qq(gwas_result, os.path.join(outputdir, "{}_common_qq.png".format(fileprefix)))
        sys.stdout.write("ok ({:f} s)\n".format(time.time() - start))

        # any
        sys.stdout.write("... writing any interaction results ... ")
        start = time.time()
        gwas_result = res.GWASResult(self.chr_names, self.chromosomes,
                                     self.used_snp_pos, pvalues_inter[2],
                                     dict(mafs=self.mafs, macs=self.macs),
                                     additional_columns={})
        gwas_result.save_as_csv(os.path.join(outputdir, "{}_any_pvals.csv".format(fileprefix)))
        gplt.plot_gwas_result(gwas_result, os.path.join(outputdir, "{}_any_manhattan.png".format(fileprefix)),
                              mac=self.mac_thres)
        gplt.plot_qq(gwas_result, os.path.join(outputdir, "{}_any_qq.png".format(fileprefix)))
        sys.stdout.write("ok ({:f} s)\n".format(time.time() - start))
        sys.stdout.write("ok.\n")


def run_by_environment_vars():
    tfile1 = os.environ['MTL_FILE1']
    tfile2 = os.environ['MTL_FILE2']
    tcols1str = os.environ['MTL_COLS1']
    tcols2str = os.environ['MTL_COLS2']
    tprefix1 = os.environ['MTL_PREFIX1']
    tprefix2 = os.environ['MTL_PREFIX2']
    snpsdb = os.environ['MTL_SNPS']
    macthres = os.environ['MTL_MAC']
    outputdir = os.environ['MTL_OUTDIR']
    jobid = int(os.getenv('PBS_ARRAY_INDEX', '0'))
    filesep = os.getenv('MTL_FILE_SEPARATOR', '\t')

    tcols1 = eval(tcols1str)
    tcols2 = eval(tcols2str)
    # tcols1 = [int(x) for x in tcols1str.lstrip('[').rstrip(']').split(',')]
    # tcols2 = [int(x) for x in tcols2str.lstrip('[').rstrip(']').split(',')]

    mt = MTL(int(macthres))
    mt.read_phenotype_col(tfile1, tcols1[jobid], tprefix1, sep=filesep)
    mt.read_phenotype_col(tfile2, tcols2[jobid], tprefix2, sep=filesep)
    mt.read_genotypes(snpsdb)
    mt.do_qtl(outputdir)


if __name__ == "__main__":
    # run_by_environment_vars()

    workdir = "/home/GMI/christian.goeschl/devel/pycharm/mtmm-limix-take"
    genotypedir = "/data/gwas/genotypes_for_pygwas/1.0.0/regmap_horton_et_al_2012"

    limtmm = MTL(mac_thres=0)
    i = 1
    j = 1
    # limtmm.read_phenotype_col(os.path.join(workdir, "bao_Std.txt"), i, colprefix="ctrl{:d}".format(i), sep="\t")
    # limtmm.read_phenotype_col(os.path.join(workdir, "bao_Cd+.txt"), j, colprefix="cd+{:d}".format(j), sep="\t")
    limtmm.read_phenotype_col(os.path.join(workdir, "LRD_IAAvsMock.csv"), i, colprefix="IAA{:d}".format(i), sep=",")
    limtmm.read_phenotype_col(os.path.join(workdir, "LRD_IBAvsMock.csv"), j, colprefix="IBA{:d}".format(j), sep=",")
    # limtmm.write_phenotypes(os.path.join(workdir, "used_phenotypes_dbg_{}-{}.csv".format(i, j)))
    limtmm.read_genotypes(os.path.join(genotypedir, "all_chromosomes_binary.hdf5"))
    limtmm.do_qtl(os.path.join(workdir, "IBA-vs-IAA-mtl-run"))

    # for i in range(3, 23):
    #     limtmm = MtmmLimix(mac_thres=5)
    #     limtmm.read_phenotype_col(os.path.join(workdir, "20170301_Zn_MS_RLd2.csv"), i, colprefix="s{:d}".format(i), sep=",")
    #     limtmm.read_phenotype_col(os.path.join(workdir, "20170301_Zn_MS_RLd2.csv"), 8, colprefix="c{:d}".format(i), sep=",")
    #     limtmm.read_genotypes(os.path.join(workdir, "all_chromosomes_binary.hdf5"))
    #     limtmm.do_qtl(os.path.join(workdir, "debug-strigo-vs-ctrl"))
    #
    #     limtmm = MtmmLimix(mac_thres=5)
    #     limtmm.read_phenotype_col(os.path.join(workdir, "GWASinput_2016_strigolactone_means-na.csv"), i, colprefix="s{:d}".format(i), sep=",")
    #     limtmm.read_phenotype_col(os.path.join(workdir, "GWASinput_2016_control_means-na.csv"), 8, colprefix="c{:d}".format(i), sep=",")
    #     limtmm.read_genotypes(os.path.join(workdir, "all_chromosomes_binary.hdf5"))
    #     limtmm.do_qtl(os.path.join(workdir, "debug-strigo-vs-ctrl"))
