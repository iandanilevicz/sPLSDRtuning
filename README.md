# Challenging estimation of tuning parameters in sPLSDR model to predict the incidence of rare diseases
**Ian Meneghel Danilevicz**, **Sam Vidil**, **Séverine Sabia**  
*2025-10-06*

[![Creative Commons License](https://i.creativecommons.org/l/by/4.0/80x15.png)](https://creativecommons.org/licenses/by/4.0/)  
This work is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

---

### Authors

- **Ian Meneghel Danilevicz**  
  Université Paris Cité, Inserm, CRESS, France  

- **Sam Vidil**  
  Université Paris Cité, Inserm, CRESS, France  

- **Séverine Sabia**  
  Université Paris Cité, Inserm, CRESS, France  
  University College London, Division of Psychiatry, UCL Brain Sciences, London, United Kingdom  

---

### Abstract

We investigate some practical points of the sparse partial least squares deviance residuals (sPLSDR) model, 
  a machine learning (ML) method with important application potential and a powerful alternative for Cox regression 
  analysis when dimension reduction and variable selection are required. It can be used in a broad range of 
  contexts, among which the health domain. This method requires to tune two hyper-parameters: eta (sparsity level) 
  and D (number of components). Usually, the literature recommends cross-validation (CV) with a number of folds 
  that depends on the sample size to estimate the hyperparameters. However, in case of dichotomous outcome for 
  which the ratio of cases and non-cases is small this approach appears to be unstable. This is the case in rare 
  disease such as incidence of major neurocognitive disorder (NCD). In that case, even large samples may perform 
  as poorly as small samples. We specifically test this scenario of large samples with rare events using Monte 
  Carlo simulations with different performance criteria (C-index and AIC). We observed that five folds are 
  more conservative (biased) than one fold, and it presents more variability (less precise) in case of rare 
  event outcomes. Previous authors suggested that Bootstrap or leave-one-out (LOO) could be a wise alternative. 
  We checked Bootstrap with one fold in our real analysis application using the Whitehall II accelerometer 
  sub-study (N total = 3,969,  cases = 291 (7.3%)) as a training set and the UK Biobank as a validation set 
  (N total = 54,833,  cases = 750 (1.4\%)), both considerably large sample size studies, with incidence of 
  major neurocognitive disorders as an example of rare outcomes. This application shows that CV is sensitive 
  to the number of folds and the initial seed. The counterpart Bootstrap was more robust and displayed two 
  components with eta equal to 0.75 as a consistent optimum choice for the tuning parameters. However, this 
  approach requires external validation, which is not always feasible. Further studies can better map the 
  ideal number of folds for CV depending on the sample size and the ratio of cases/non-cases. As dichotomous 
  outcomes and survival time are classical statistical subjects, the ML literature provide clear and precise 
  recommendations for several available dawn methods.
