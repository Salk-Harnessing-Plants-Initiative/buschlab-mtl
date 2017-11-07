import logging
import os
import pandas as pd
import scipy as sp
import h5py as h5
import pygwas_modules.mtcorr as mtcorr
import pygwas_modules.statistics as stats

from core.data_transform import DataTransform

log = logging.getLogger('__name__')


class Phenotype(object):
    def __init__(self):
        self.data = pd.DataFrame()
        self.transform = ""

    @property
    def data(self):
        return self.__data

    @data.setter
    def data(self, data):
        self.__data = data

    @property
    def transform(self):
        return self.__transform

    @transform.setter
    def transform(self, transform):
        self.__transform = transform

    def read_csv_col(self, phenotype_filepath, colnr, colprefix="", sep='\t'):
        log.info("reading phenotypes: {}, colnr: {}".format(phenotype_filepath, colnr))

        with open(phenotype_filepath, 'r') as phenofile:
            hcols = phenofile.readline().strip().split(sep)

            ids = []
            vals = []
            for l in phenofile:
                dcols = l.strip().split(sep)
                try:
                    vals.append(sp.float64(dcols[colnr]))
                    ids.append(dcols[0])
                except ValueError:
                    log.warn("excluding accession {} because of trait value: {}".format(dcols[0], dcols[colnr]))
                    continue

            tmp_data = pd.DataFrame(vals)
            tmp_data.index = ids
            tmp_data.columns = ["_".join([colprefix, hcols[colnr]])]

            if self.data.size == 0:
                self.data = tmp_data
            else:
                self.data = pd.concat([self.data, tmp_data], join='inner', axis=1)
                log.info("phenotype intersection is {} accessions.".format(self.data.shape[0]))
        # self.data.sort_index(inplace=True, key=float);
        return


class GwasData(object):
    def __init__(self):
        self.__data_h5 = None
        # self.data = pd.DataFrame()

    @property
    def data_h5(self):
        return self.__data_h5

    @data_h5.setter
    def data_h5(self, data_h5):
        self.__data_h5 = data_h5

    def read_csv(self, path):
        self.__data_h5 = h5.File('memfile.h5', 'w', driver='core', backing_store=False)
        self.__data_h5.create_group('/pvalues')
        self.__data_h5.create_group('/quantiles')
        csvdata = pd.read_csv(path, header=0)

        #convert chromosome numbers to strings if not already done.
        if csvdata['chromosomes'].dtype == sp.dtype('int'):
            csvdata['chromosomes'] = ["chr{:d}".format(x) for x in csvdata['chromosomes']]
        chr_names = sorted(set(csvdata['chromosomes'].values))

        csvdata['scores'] = -sp.log10(csvdata['pvals'].values)

        for chr_name in chr_names:
            subdata = csvdata[csvdata['chromosomes'] == chr_name]
            h5_chr_group = self.__data_h5.create_group('/pvalues/{}'.format(chr_name))
            for col in subdata.columns:
                if col == 'chromosomes' or col == 'pvals':
                    continue
                h5_chr_group.create_dataset(col, data=subdata[col].values)

        # create attributes
        h5_pvalues_group = self.__data_h5['/pvalues']
        h5_pvalues_group.attrs['analysis_method'] = 'N/A'
        h5_pvalues_group.attrs['bh_thres'] = -sp.log10(mtcorr.get_bhy_thres(csvdata['pvals'].values)['thes_pval'])
        h5_pvalues_group.attrs['bonferroni_threshold'] = -sp.log10(0.05 / csvdata.shape[0])
        ks_stats = stats.calc_ks_stats(csvdata['pvals'].values)
        h5_pvalues_group.attrs['ks_pval'] = ks_stats['p_val']
        h5_pvalues_group.attrs['ks_stat'] = ks_stats['D']
        h5_pvalues_group.attrs['max_score'] = csvdata['scores'].values.max()
        h5_pvalues_group.attrs['med_pval'] = sp.median(csvdata['pvals'].values)
        h5_pvalues_group.attrs['numberOfSNPs'] = csvdata.shape[0]
        h5_pvalues_group.attrs['transformation'] = 'raw'

    def write_hdf5(self, path):
        with h5.File(path, 'w') as real_h5:
            for item in self.__data_h5:
                real_h5.copy(self.__data_h5[item], item, expand_refs=True)
            for attr_name in self.__data_h5.attrs:
                real_h5.attrs[attr_name] = self.__data_h5.attrs[attr_name]


if __name__ == "__main__":
    # log = logging.getLogger()
    # log.basicConfig(level=log.INFO, format='%(asctime)s %(levelname)s: %(message)s')
    # workdir = "/data/christian.goeschl/mtmm"
    # pheno = Phenotype()
    # pheno.read_csv_col(os.path.join(workdir, "bao_Std.txt"), colnr=1, colprefix="ctrl-{:d}".format(1), sep="\t")
    # pheno.read_csv_col(os.path.join(workdir, "bao_Cd+.txt"), colnr=1, colprefix="ctrl-{:d}".format(1), sep="\t")
    # DataTransform.transform(pheno.data.values, 'most-normal')
    gwas_data = GwasData()
    gwas_data.read_csv('/net/gmi.oeaw.ac.at/busch/lab_new/Christian/mtl-tempstress/10LT-mtl-250k-results/LT_median_Total_length_day001-x-Cli_Bio01_Annu-mac1_any_pvals.csv')
    gwas_data.write_hdf5('/home/GMI/christian.goeschl/LT_median_Total_length_day001-x-Cli_Bio01_Annu-mac1_any_pvals.hdf5')
    pass
    # pheno.sqr_transform()
