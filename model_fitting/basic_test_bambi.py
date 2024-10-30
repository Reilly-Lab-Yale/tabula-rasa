
import argparse
import bambi as bmb
import arviz as az
import pymc
import matplotlib.pyplot as plt
import pandas as pd



def poisson_model(counts, patsy_formula):
    model = bmb.Model(patsy_formula, counts, family = "poisson")

    return model

def zi_poisson_model(counts, patsy_formula):
    model = bmb.Model(patsy_formula, counts, family = "zero_inflated_poisson")

    return model

def negative_binomial_model(counts, patsy_formula):
    model = bmb.Model(patsy_formula, counts, family = "negativebinomial")

    return model

def zi_negative_binomial_model(counts, patsy_formula):
    model = bmb.Model(patsy_formula, counts, family = "zero_inflated_negativebinomial")

    return model

def get_stats(model, fit_model, outpath, save_name):
    print('getting stats')
    model_sum = az.summary(fit_model)
    model.predict(fit_model, kind="response")
    model_fitting = fit_model.to_dataframe()

    model_sum.to_csv("%s/%s_summary.csv" % (outpath, save_name))
    model_fitting.to_csv("%s/%s_fit_model_df.csv" % (outpath, save_name))

    az.plot_ppc(fit_model, show=True)
    plt.savefig("%s/%s_ppc_plot.png" % (outpath, save_name))
    return model_sum, model_fitting


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scmpra_counts_file', type=str)
    parser.add_argument('--model_choice', type=str)
    parser.add_argument('--formula', type=str)
    parser.add_argument('--temp_dir')
    parser.add_argument('--out_file', type=str)
    args = parser.parse_args()
    scmpra_counts = pd.read_table(args.scmpra_counts_file)
    formula = args.formula
    print('formula: %s' % formula)
    maxiter = args.maxiter
    print('maxiter: %s' % maxiter)
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


    scmpra_model_fit = scmpra_model.fit(draws = maxiter)

    get_stats(scmpra_model, scmpra_model_fit, temp_dir, out_file)



    return










if __name__ == "__main__":
    main()