# Calibration of Breast Cancer Natural History Models Using Approximate Bayesian Computation
This is the repository for my master thesis "Calibration of Breast Cancer Natural History Models Using Approximate Bayesian Computation" in mathematical statistics at KTH Royal Institute of Technology. The master thesis considers lieklihood-free estimation of the breast cancer natural history models developed by Keith Humphreys, Linda Abrahamsson, Gabriel Isheden and Rickard Strandberg at Karolinska Institutet, MEB.  

## Abstract 
Natural history models for breast cancer describe the unobservable disease progression. These models can either be fitted using likelihood-based estimation to data on individual tumour characteristics, or calibrated to fit statistics at a population level. Likelihood-based inference using individual level data has the advantage of ensuring model parameter identifiability. However, the likelihood function can be computationally heavy to evaluate or even intractable.

In this thesis likelihood-free estimation using Approximate Bayesian Computation (ABC) will be explored. The main objective is to investigate whether ABC can be used to fit models to data collected in the presence of mammography screening. As a background, a literature review of ABC is provided.

As a first step an ABC-MCMC algorithm is constructed for two simple models both describing populations in absence of mammography screening, but assuming different functional forms of tumour growth. The algorithm is evaluated for these models in a simulation study using synthetic data, and compared with results obtained using likelihood-based inference.

Later, it is investigated whether ABC can be used for the models in presence of screening. The findings of this thesis indicate that ABC is not directly applicable to these models. However, by including a sub-model for tumour onset and assuming that all individuals in the population have the same screening attendance it was possible to develop an ABC-MCMC algorithm that carefully takes individual level data into consideration in the estimation procedure. Finally, the algorithm was tested in a simple simulation study using synthetic data.

Future research is still needed to evaluate the statistical properties of the algorithm (using extended simulation) and to test it on observational data where previous estimates are available for reference. 

## Contents of the repository
The report and the code is included in this repository. In the folder code/models-in-absence-of-screening the code for ABC for the models in absence of screening, described in Chapter 3 of the report, is included. The folder code/models-in-presence-of-screening contains the code for ABC for the models in presence of screening (Chapter 4).