import pyarrow.parquet as pq
import argparse
import statsmodels.discrete.count_model as smdc

class DataSet(dict):
    def __init__(self, path):
        print('initializing dataset')
        self.filepath = path
        self.parquet = pq.ParquetFile(self.filepath)
    
    def __getitem__(self, key):
        try:
            return self.parquet.read([key]).to_pandas()[key]
        except:
            raise KeyError

    def __reduce__(self):
        #return self.parquet.read().to_pandas().__reduce__()
        return (self.__class__, (self.filepath, ))


def poisson_model(counts_parq, patsy_formula):
    counts_model_poisson = smdc.GeneralizedPoisson.from_formula(formula = patsy_formula, data = counts_parq)

    return counts_model_poisson

def zi_poisson_model(counts_parq, patsy_formula):
    counts_model_zi_poisson = smdc.ZeroInflatedPoisson.from_formula(formula = patsy_formula, data = counts_parq)

    return counts_model_zi_poisson

def negative_binomial_model(counts_parq, patsy_formula):
    counts_model_negative_binomial = smdc.NegativeBinomialP.from_formula(formula = patsy_formula, data = counts_parq)

    return counts_model_negative_binomial

def zi_negative_binomial_model(counts_parq, patsy_formula):
    counts_model_zi_negative_binomial = smdc.ZeroInflatedNegativeBinomialP.from_formula(formula = patsy_formula, data = counts_parq)

    return counts_model_zi_negative_binomial

def get_stats(fit_model):
    print('getting stats')
    # available stats come from here https://www.statsmodels.org/dev/generated/statsmodels.discrete.discrete_model.CountResults.html

    aic = fit_model.aic
    bic = fit_model.bic
    loglike = fit_model.llf
    llr_chi2 = fit_model.llr
    llr_chip = fit_model.llr_pvalue
    pseudr = fit_model.prsquared
    resid = fit_model.resid
    converged = fit_model.converged

    return [converged, aic, bic, loglike, llr_chi2, llr_chip, pseudr, resid]

def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scmpra_counts_file', type=str)
    parser.add_argument('--model_choice', type=str)
    parser.add_argument('--formula', type=str)
    parser.add_argument('--maxiter', type=int, default=50)
    parser.add_argument('--regularized_fit', type=boolean_string, default=False)
    parser.add_argument('--temp_dir')
    parser.add_argument('--out_file', type=str)
    args = parser.parse_args()
    scmpra_counts = DataSet(args.scmpra_counts_file)
    formula = args.formula
    print('formula: %s' % formula)
    maxiter = args.maxiter
    print('maxiter: %s' % maxiter)
    reg_fit = args.regularized_fit
    print('reg_fit: %s' % reg_fit)
    temp_dir = args.temp_dir
    print('temp_dir: %s' % temp_dir)
    model_choice = args.model_choice
    print('model_choice: %s' % model_choice)
    out_file = args.out_file
    print('out_file: %s' % out_file)

    model_dict = {'poisson': poisson_model,
                  'zi_poisson' : zi_poisson_model, 
                  'negative_binomial' : negative_binomial_model,
                  'zi_negative_binomial' : zi_negative_binomial_model}
    

    # try:
    scmpra_model = model_dict[model_choice](scmpra_counts, formula)
    # except:
    #     print('Failed to build %s model' % model_choice)
    #     return
    
    print(model_choice)

    if reg_fit:
        scmpra_model_fit = scmpra_model.fit_regularized(maxiter = maxiter)

    else:
        scmpra_model_fit = scmpra_model.fit(maxiter = maxiter)

    scmpra_model_fit.save("%s/%s_fit_model.pickle" % (temp_dir, out_file))

    model_info = [model_choice, formula, maxiter, reg_fit]
    model_stats = get_stats(scmpra_model_fit)
    print(model_stats)
    print(model_info)
    out_list = model_info + model_stats

    with open("%s/%s_stats.txt" % (temp_dir, out_file), "w") as o:
        o.write("\t".join(str(x) for x in out_list))

    return










if __name__ == "__main__":
    main()