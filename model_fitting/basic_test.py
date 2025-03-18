import pyarrow.parquet as pq
import argparse
import statsmodels.discrete.count_model as smdc
import pandas as pd
import numpy as np

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
    # counts_model_poisson = smdc.GeneralizedPoisson.from_formula(formula = patsy_formula, data = counts_parq)

    # return counts_model_poisson
    return

def zi_poisson_model(counts_parq, patsy_formula):
    # counts_model_zi_poisson = smdc.ZeroInflatedPoisson.from_formula(formula = patsy_formula, data = counts_parq)

    # return counts_model_zi_poisson
    return

def negative_binomial_model(counts_parq, patsy_formula):
    # counts_model_negative_binomial = smdc.NegativeBinomialP.from_formula(formula = patsy_formula, data = counts_parq)

    # return counts_model_negative_binomial
    return

def zi_negative_binomial_model(counts_parq, patsy_formula, zi_param):
    counts_model_zi_negative_binomial = smdc.ZeroInflatedNegativeBinomialP.from_formula(formula = patsy_formula, data = counts_parq, exog_infl = pd.get_dummies(pd.Categorical(counts_parq.__getitem__(zi_param))))

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

    return [converged, aic, bic, loglike, llr_chi2, llr_chip, pseudr, list(resid)]

def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scmpra_counts_file', type=str)
    parser.add_argument('--model_choice', type=str)
    parser.add_argument('--formula', type=str)
    parser.add_argument('--maxiter', type=int, default=200)
    parser.add_argument('--zi_param', type=str)
    parser.add_argument('--temp_dir')
    parser.add_argument('--out_file', type=str)
    args = parser.parse_args()
    scmpra_counts = DataSet(args.scmpra_counts_file)
    formula = args.formula
    print('formula: %s' % formula)
    maxiter = args.maxiter
    print('maxiter: %s' % maxiter)
    zi_param = args.zi_param
    print('zi_param: %s' % zi_param)
    temp_dir = args.temp_dir
    print('temp_dir: %s' % temp_dir)
    model_choice = args.model_choice
    print('model_choice: %s' % model_choice)
    out_file = args.out_file
    print('out_file: %s' % out_file)
    count = out_file.split('_')[0]

    model_dict = {'poisson': poisson_model,
                  'zi_poisson' : zi_poisson_model, 
                  'negative_binomial' : negative_binomial_model,
                  'zi_negative_binomial' : zi_negative_binomial_model}
    

    # try:
    scmpra_model = model_dict[model_choice](scmpra_counts, formula, zi_param)
    # except:
    #     print('Failed to build %s model' % model_choice)
    #     return
    
    print(model_choice)
    
    n_count_params = scmpra_model.exog.shape[1]      # Count model parameters
    n_infl_params = scmpra_model.exog_infl.shape[1]    # Inflation model parameters
    n_total = n_count_params + n_infl_params + 1 # adding 1 for alpha
    start_params = np.full(n_total, 0.1)
    
    scmpra_model_fit = scmpra_model.fit(start_params=start_params, method="bfgs",maxiter=maxiter)

    scmpra_model_fit.save("%s/%s_fit_model.pickle" % (temp_dir, out_file))

    model_info = [count,model_choice, formula, maxiter, zi_param]
    model_stats = get_stats(scmpra_model_fit)
    print(model_stats)
    print(model_info)
    out_list = model_info + model_stats[:-1]

    with open("%s/%s_stats.txt" % (temp_dir, out_file), "w") as o:
        o.write("\t".join(str(x) for x in out_list))
        o.write("\n")


    return










if __name__ == "__main__":
    main()