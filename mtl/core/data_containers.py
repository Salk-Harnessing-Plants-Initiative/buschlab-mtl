import logging
import os
import pandas as pd
import scipy as sp

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


class GWAS_Data(object):
    def __init__(self):
        self.data = pd.DataFrame()

    @property
    def data(self):
        return self.__data

    @data.setter
    def data(self, data):
        self.__data = data

if __name__ == "__main__":
    # log = logging.getLogger()
    # log.basicConfig(level=log.INFO, format='%(asctime)s %(levelname)s: %(message)s')
    workdir = "/data/christian.goeschl/mtmm"
    pheno = Phenotype()
    pheno.read_csv_col(os.path.join(workdir, "bao_Std.txt"), colnr=1, colprefix="ctrl-{:d}".format(1), sep="\t")
    pheno.read_csv_col(os.path.join(workdir, "bao_Cd+.txt"), colnr=1, colprefix="ctrl-{:d}".format(1), sep="\t")
    # DataTransform.transform(pheno.data.values, 'most-normal')
    pass
    # pheno.sqr_transform()
