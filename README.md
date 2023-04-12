# BIOS761-GROUP2

The following contains an R package titled "ZINB" which will be applied to a dataset on Kaggle containing information on a child's microbiome and whether or not a child has Autism Spectrum Disorder (ASD). The ultimate goal of the package is to build a model using microbiome differential abundance in species to predict whether or not someone has ASD.

This will be done through fitting a ZINB (Zero Inflated Negative Binomial Model) to find differential abundance in species using the EM algorithm, then using this differential abundance in logistic regression to predict ASD status.

Model performance will be compared with a Random Forest model which will also predict ASD status with the all microbiome species.