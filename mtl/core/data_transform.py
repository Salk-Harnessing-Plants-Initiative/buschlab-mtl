import scipy as sp
import numpy as np




class DataTransform(object):

    @staticmethod
    def transform(data, method):
        trafo_mapping = {'log':             DataTransform.__log_transform,
                         'exp':             DataTransform.__exp_transform,
                         'sqrt':            DataTransform.__sqrt_transform,
                         'sqr':             DataTransform.__sqr_transform,
                         'ascombe':         DataTransform.__ascombe_transform,
                         'arcsin-sqrt':     DataTransform.__arcsin_sqrt_transform,
                         'box-cox':         DataTransform.__box_cox_transform,
                         'feature-scale':   DataTransform.__feature_scale,
                         'zscore-scale':    DataTransform.__zscore_scale}

        if method != 'most-normal':
            func1d = trafo_mapping[method]
            data[:] = np.apply_along_axis(func1d, axis=0, arr=data)[:]
        else:
            data[:] = DataTransform.__most_normal_transform(data)[:]
        return

    @staticmethod
    def __feature_scale(a):
        return (a - a.min()) / (a.max() - a.min())

    @staticmethod
    def __zscore_scale(a):
        return (a - a.mean()) / a.std()

    @staticmethod
    def __log_transform(a):
        return sp.log((a - a.min()) + 0.1 * a.std())

    @staticmethod
    def __sqrt_transform(a):
        return sp.sqrt((a - a.min()) + 0.1 * a.std())

    @staticmethod
    def __ascombe_transform(a):
        return 2.0 * sp.sqrt(a + 3.0 / 8.0)

    @staticmethod
    def __sqr_transform(a):
        p = ((a - a.min()) + 0.1 * a.std())
        return p * p

    @staticmethod
    def __exp_transform(a):
        return sp.exp((a - a.min()) + 0.1 * a.std())

    @staticmethod
    def __arcsin_sqrt_transform(a):
        if a.min() < 0 or a.max() > 1:
            # log.debug('Some values are outside of range [0,1], using feature scaling!')
            a = DataTransform.__feature_scale(a)
        return sp.arcsin(sp.sqrt(a))

    @staticmethod
    def __box_cox_transform(a):
        from scipy import stats
        boxcox, maxlog = stats.boxcox(a, lmbda=None, alpha=None)
        return boxcox

    @staticmethod
    def __most_normal_transform(data, trafo_types=None):
        from scipy import stats
        if trafo_types is None:
            trafo_types = ['none', 'sqrt', 'log', 'sqr', 'exp', 'arcsin-sqrt']

        shapiro_pvals = []
        for trafo_type in trafo_types:
            dcpy = data.copy()
            if trafo_type != 'none':
                DataTransform.transform(dcpy, trafo_type)
            pvals = []
            for col in range(dcpy.shape[1]):
                r = stats.shapiro(dcpy[:, col])
                if sp.isfinite(r[0]):
                    pvals.append(r[1])
                else:
                    pvals.append(0.0)
            shapiro_pvals.append(pvals)

        sa = sp.array(shapiro_pvals)
        return
