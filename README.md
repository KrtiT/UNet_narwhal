# Buzz_accelerometer_paper

Some codes in Python/Pytorch and R for the analysis in the paper " Detection of foraging behavior from accelerometer data using U-Net type convolutional networks" (https://www.sciencedirect.com/science/article/pii/S1574954121000662).
Free preprint availabe at https://arxiv.org/abs/2101.01992.

We use supervised machine learning models and logistic regression to predict the buzz,

File U_Time_7.ipynb is one of jupyter notebooks of U-Net deep learning models, using cross validation. 
There we trained on 3 whales, validate on one whale, test on the remaining whale to evaluate
the quality of the prediction vs ground truth buzzes.

File RF_all_OneHot.py consist the code of fitting Random Forest and Logisitic Regression model.

File jerk_buzz.R shows the code of jerk analysis, and plotting figures in the paper.
