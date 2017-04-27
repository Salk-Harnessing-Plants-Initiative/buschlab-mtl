import logging
import os
import pandas as pd
import scipy as sp

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

    def log_transform(self):
        for i in range(self.data.shape[1]):
            p = self.data.iloc[:, i].values
            p = sp.log((p - p.min()) + 0.1 * sp.std(p))
            self.data.iloc[:, i] = p
        self.transform = "log({})".format(self.transform)

    def sqrt_transform(self):
        for i in range(self.data.shape[1]):
            p = self.data.iloc[:, i].values
            p = sp.sqrt((p - p.min()) + 0.1 * sp.std(p))
            self.data.iloc[:, i] = p
        self.transform = "sqrt({})".format(self.transform)

    def ascombe_transform(self):
        for i in range(self.data.shape[1]):
            p = self.data.iloc[:, i].values
            p = 2.0 * sp.sqrt(p + 3.0 / 8.0)
            self.data.iloc[:, i] = p
        self.transform = "ascombe({})".format(self.transform)

    def sqr_transform(self):
        for i in range(self.data.shape[1]):
            p = self.data.iloc[:, i].values
            p = ((p - p.min())+0.1*p.std())
            p = p*p
            self.data.iloc[:, i] = p
        self.transform = "sqr({})".format(self.transform)

    def exp_transform(self):
        for i in range(self.data.shape[1]):
            p = self.data.iloc[:, i].values
            p = sp.exp((p - p.min()) + 0.1 * p.std())
            self.data.iloc[:, i] = p
        self.transform = "exp({})".format(self.transform)

    def box_cox_transform(self, lambda_range=(-2.0, 2.0), lambda_increment=0.1):
        from scipy import stats
        for i in range(self.data.shape[1]):
            p = self.data.iloc[:, i].values
            p = ((p - p.min())+0.1*p.std())

        sw_pvals = []
        lambdas = sp.arange(lambda_range[0], lambda_range[1] + lambda_increment, lambda_increment)
        for l in lambdas:
            if l == 0:
                vs = sp.log(p)
            else:
                vs = ((p ** l) - 1) / l
            r = stats.shapiro(vs)
            if sp.isfinite(r[0]):
                pval = r[1]
            else:
                pval = 0.0
            sw_pvals.append(pval)
        log.info(sw_pvals)
        i = sp.argmax(sw_pvals)
        l = lambdas[i]
        if l == 0:
            vs = sp.log(p)
        else:
            vs = ((p ** l) - 1) / l
        self.transform = "box-cox({})".format(self.transform)
        log.debug('optimal lambda was %0.1f' % l)

        #
        # p1 = np.array(self.phenotypes[[0]].loc[self.iid]).astype(np.float64)
        # p1 = np.log((p1 - p1.min()) + 0.1 * np.std(p1))
        # p2 = np.array(self.phenotypes[[1]].loc[self.iid]).astype(np.float64)
        # p2 = np.log((p2 - p2.min()) + 0.1 * np.std(p2))
        # pheno_norm = np.concatenate((p1, p2), axis=1)


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


if __name__ == "__main__":
    # log = logging.getLogger()
    # log.basicConfig(level=log.INFO, format='%(asctime)s %(levelname)s: %(message)s')

    workdir = "/data/christian.goeschl/mtmm"
    pheno = Phenotype()
    pheno.read_csv_col(os.path.join(workdir, "bao_Std.txt"), colnr=1, colprefix="ctrl-{:d}".format(1), sep="\t")
    pheno.read_csv_col(os.path.join(workdir, "bao_Cd+.txt"), colnr=1, colprefix="ctrl-{:d}".format(1), sep="\t")
    pheno.sqr_transform()
