# K-fold cross validation

from sklearn.linear_model import LogisticRegression
import numpy as np
import pandas as pd
import argparse
from scipy import stats

#############
#___SETUP___#
#############

parser = argparse.ArgumentParser(description=__doc__)
# add a positional writable file argument
parser.add_argument('input_csv_file', type=argparse.FileType('r'))
# parser.add_argument('output_csv_file', type=argparse.FileType('w')) 
args = parser.parse_args()

# Read in classification training data
csv_data = pd.read_csv(args.input_csv_file, parse_dates = ['datetime'], 
                        infer_datetime_format = True) #Read as DateTime obsject

# add "month of year" (MoY) to dataset 
csv_data['MoY'] = [int(x.day) for x in csv_data['datetime']]

# set randomness seed
np.random.seed(420)

# make training dataset model
# create data points for training algorithm
# 1 = classA (EAC), 0 = classB (TS)
var1 = list(csv_data['temp'])
var2 = list(csv_data['salt'])
MoY = list(csv_data['MoY'])
water_class = list(csv_data['class'])
# make data frame
train_data = {'var1': var1, 'var2': var2, 'MoY': MoY, 'class': water_class}
train_data = pd.DataFrame(data=train_data)
# replace current data strings with binary integers
train_data['class'] = train_data['class'].replace(to_replace='EAC', value=1)
train_data['class'] = train_data['class'].replace(to_replace='BS', value=0)
# shuffle dataset
train_data = train_data.iloc[np.random.permutation(np.arange(len(train_data)))]
train_data = train_data.reset_index(drop=True)

# split data set into 10
train_data = np.array_split(train_data, 10)

######################
#___Validate Model___#
######################

# result list
results = []

# run classification on each partition
for i in range(0,10):
    print('Model run ' + str(i+1) + ' of 10')
    model_dat = train_data[i]
    train_dat = train_data.remove(i)
    # fit logistic regression to the training data
    lr_model = LogisticRegression()
    lr_model = lr_model.fit(train_dat[['var1','var2','MoY']], np.ravel(train_dat[['class']]))
    # run predictive model
    prob = lr_model.predict_proba(model_dat[['var1','var2','MoY']])
    _, prob = zip(*probs)
    prob = list(probs)
    df = {'prob': prob, 'class': model_dat['class']}
    df = pd.DataFrame(data=df)
    df['result'] = [1 if x >= 0.5 else 0 for x in df['prob']]
    # test results
    valid_result = 0
    for idx, row in df.iterrows():
        if int(row['class']) == int(row['result']):
            valid_result = valid_result + 1
    valid_result = valid_result/len(df['class'])
    # add to resuklts list
    results.append(valid_result)

print('Here are the reuslts of the 10 runs...')
print('results')
print('Average accuary was ' + str(sum(results)/len(results)))























