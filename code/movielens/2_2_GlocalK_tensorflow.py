#!/usr/bin/env python
# coding: utf-8
# ============================================================================
# GLocal-K — Original Implementation
# Requires: Python 3.7, TensorFlow 1.15 
# Reference: https://github.com/usydnlp/Glocal_K/
# ============================================================================

import os
from time import time
from scipy.stats import spearmanr
import numpy as np
import pandas as pd
import tensorflow as tf

# Navigate to the project root directory
script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in dir() else os.getcwd()
project_root = os.path.abspath(os.path.join(script_dir, '..', '..', '..'))
os.chdir(project_root)
print(f"Working directory: {os.getcwd()}")
print(f"TensorFlow version: {tf.__version__}")

# Force TensorFlow 1 to use a single CPU core (equivalent to PyTorch benchmark setup)
config = tf.ConfigProto(
    intra_op_parallelism_threads=1,
    inter_op_parallelism_threads=1,
    device_count={'CPU': 1, 'GPU': 0}
)
# Common hyperparameter settings
n_hid = 500 # size of hidden layers
n_dim = 5 # inner AE embedding size
n_layers = 2 # number of hidden layers
gk_size = 3 # width=height of kernel for convolution

# Hyperparameters to tune for specific case
lambda_2 = 70.  # regularisation of number of parameters
lambda_s = 0.018 # regularisation of sparsity of the final matrix
iter_p = 50  # max iterations per L-BFGS-B call for pretraining
iter_f = 10  # max iterations per L-BFGS-B call for finetuning
epoch_p = 20 # number of pretraining epochs (outer loops)
epoch_f = 30 # number of finetuning epochs (outer loops)

try:
    import config
    if hasattr(config, 'epoch_p'): epoch_p = config.epoch_p
    if hasattr(config, 'epoch_f'): epoch_f = config.epoch_f
except ImportError:
    pass
dot_scale = 0.5 # dot product weight for global kernel


# ---------------------------------------------------------------------------
# Data Loader Functions
# ---------------------------------------------------------------------------
def load_data_1m(path='./', delimiter='::', frac=0.1, seed=1234):
    tic = time()
    print('reading data...')
    data = np.loadtxt(path + 'Movie_Y_Glocalk.dat', skiprows=0, delimiter=delimiter).astype('int32')
    print('taken', time() - tic, 'seconds')

    n_u = np.unique(data[:,0]).size
    n_m = np.unique(data[:,1]).size
    n_r = data.shape[0]

    udict = {}
    for i, u in enumerate(np.unique(data[:,0]).tolist()):
        udict[u] = i
    mdict = {}
    for i, m in enumerate(np.unique(data[:,1]).tolist()):
        mdict[m] = i

    np.random.seed(seed)
    idx = np.arange(n_r)
    np.random.shuffle(idx)

    train_r = np.zeros((n_m, n_u), dtype='float32')
    test_r = np.zeros((n_m, n_u), dtype='float32')

    for i in range(n_r):
        u_id = data[idx[i], 0]
        m_id = data[idx[i], 1]
        r = data[idx[i], 2]

        if i < int(frac * n_r):
            test_r[mdict[m_id], udict[u_id]] = r
        else:
            train_r[mdict[m_id], udict[u_id]] = r

    train_m = np.greater(train_r, 1e-12).astype('float32')
    test_m = np.greater(test_r, 1e-12).astype('float32')

    print('data matrix loaded')
    print('num of users: {}'.format(n_u))
    print('num of movies: {}'.format(n_m))
    print('num of training ratings: {}'.format(n_r - int(frac * n_r)))
    print('num of test ratings: {}'.format(int(frac * n_r)))

    return n_m, n_u, train_r, train_m, test_r, test_m, udict, mdict

def load_data_1m_test(path, udict, mdict, delimiter='::'):
    data = np.loadtxt(path + 'Movie_test_Glocalk.dat', delimiter=delimiter).astype('int32')
    n_m, n_u = len(mdict), len(udict)
    test_r = np.zeros((n_m, n_u), dtype='float32')
    test_m = np.zeros((n_m, n_u), dtype='float32')
    for u_id, m_id, r in data:
        if u_id in udict and m_id in mdict:
            i, j = mdict[m_id], udict[u_id]
            test_r[i, j] = r
            test_m[i, j] = 1
    return test_m, test_r

# ---------------------------------------------------------------------------
# Load Data
# ---------------------------------------------------------------------------
data_path = './data/movielens/raw/'

try:
    path = data_path + '/'
    n_m, n_u, train_r, train_m, test_r, test_m, udict, mdict = load_data_1m(path=path, delimiter='::', frac=0.1, seed=1234)
    final_test_m, final_test_r = load_data_1m_test(path, udict, mdict, delimiter='::')
except Exception:
    print('Error: Unable to load data')

# ---------------------------------------------------------------------------
# Network Functions
# ---------------------------------------------------------------------------
R = tf.placeholder("float", [n_m, n_u])

def local_kernel(u, v):
    dist = tf.norm(u - v, ord=2, axis=2)
    hat = tf.maximum(0., 1. - dist**2)
    return hat

def kernel_layer(x, n_hid=n_hid, n_dim=n_dim, activation=tf.nn.sigmoid, lambda_s=lambda_s, lambda_2=lambda_2, name=''):
    with tf.variable_scope(name, reuse=tf.AUTO_REUSE):
        W = tf.get_variable('W', [x.shape[1], n_hid])
        n_in = x.get_shape().as_list()[1]
        u = tf.get_variable('u', initializer=tf.random.truncated_normal([n_in, 1, n_dim], 0., 1e-3))
        v = tf.get_variable('v', initializer=tf.random.truncated_normal([1, n_hid, n_dim], 0., 1e-3))
        b = tf.get_variable('b', [n_hid])

    w_hat = local_kernel(u, v)

    sparse_reg = tf.contrib.layers.l2_regularizer(lambda_s)
    sparse_reg_term = tf.contrib.layers.apply_regularization(sparse_reg, [w_hat])

    l2_reg = tf.contrib.layers.l2_regularizer(lambda_2)
    l2_reg_term = tf.contrib.layers.apply_regularization(l2_reg, [W])

    W_eff = W * w_hat  # Local kernelised weight matrix
    y = tf.matmul(x, W_eff) + b
    y = activation(y)

    return y, sparse_reg_term + l2_reg_term

def global_kernel(input, gk_size, dot_scale):
    avg_pooling = tf.reduce_mean(input, axis=1)
    avg_pooling = tf.reshape(avg_pooling, [1, -1])
    n_kernel = avg_pooling.shape[1].value

    conv_kernel = tf.get_variable('conv_kernel', initializer=tf.random.truncated_normal([n_kernel, gk_size**2], stddev=0.1))
    gk = tf.matmul(avg_pooling, conv_kernel) * dot_scale
    gk = tf.reshape(gk, [gk_size, gk_size, 1, 1])

    return gk

def global_conv(input, W):
    input = tf.reshape(input, [1, input.shape[0], input.shape[1], 1])
    conv2d = tf.nn.relu(tf.nn.conv2d(input, W, strides=[1,1,1,1], padding='SAME'))
    return tf.reshape(conv2d, [conv2d.shape[1], conv2d.shape[2]])

# ---------------------------------------------------------------------------
# Network Instantiation — Pre-training
# ---------------------------------------------------------------------------
y = R
reg_losses = None

for i in range(n_layers):
    y, reg_loss = kernel_layer(y, name=str(i))
    reg_losses = reg_loss if reg_losses is None else reg_losses + reg_loss

pred_p, reg_loss = kernel_layer(y, n_u, activation=tf.identity, name='out')
reg_losses = reg_losses + reg_loss

# L2 loss
diff = train_m * (train_r - pred_p)
sqE = tf.nn.l2_loss(diff)
loss_p = sqE + reg_losses

optimizer_p = tf.contrib.opt.ScipyOptimizerInterface(loss_p, options={'disp': True, 'maxiter': iter_p, 'maxcor': 10}, method='L-BFGS-B')

# ---------------------------------------------------------------------------
# Network Instantiation — Fine-tuning
# ---------------------------------------------------------------------------
y = R
reg_losses = None

for i in range(n_layers):
    y, _ = kernel_layer(y, name=str(i))

y_dash, _ = kernel_layer(y, n_u, activation=tf.identity, name='out')

gk = global_kernel(y_dash, gk_size, dot_scale)
y_hat = global_conv(train_r, gk)

for i in range(n_layers):
    y_hat, reg_loss = kernel_layer(y_hat, name=str(i))
    reg_losses = reg_loss if reg_losses is None else reg_losses + reg_loss

pred_f, reg_loss = kernel_layer(y_hat, n_u, activation=tf.identity, name='out')
reg_losses = reg_losses + reg_loss

# L2 loss
diff = train_m * (train_r - pred_f)
sqE = tf.nn.l2_loss(diff)
loss_f = sqE + reg_losses

optimizer_f = tf.contrib.opt.ScipyOptimizerInterface(loss_f, options={'disp': True, 'maxiter': iter_f, 'maxcor': 10}, method='L-BFGS-B')

# ---------------------------------------------------------------------------
# Evaluation helpers
# ---------------------------------------------------------------------------
def dcg_k(score_label, k):
    dcg, i = 0., 0
    for s in score_label:
        if i < k:
            dcg += (2**s[1]-1) / np.log2(2+i)
            i += 1
    return dcg

def ndcg_k(y_hat, y, k):
    score_label = np.stack([y_hat, y], axis=1).tolist()
    score_label = sorted(score_label, key=lambda d:d[0], reverse=True)
    score_label_ = sorted(score_label, key=lambda d:d[1], reverse=True)
    norm, i = 0., 0
    for s in score_label_:
        if i < k:
            norm += (2**s[1]-1) / np.log2(2+i)
            i += 1
    dcg = dcg_k(score_label, k)
    return dcg / norm

def call_ndcg(y_hat, y):
    ndcg_sum, num = 0, 0
    y_hat, y = y_hat.T, y.T
    n_users = y.shape[0]

    for i in range(n_users):
        y_hat_i = y_hat[i][np.where(y[i])]
        y_i = y[i][np.where(y[i])]

        if y_i.shape[0] < 2:
            continue

        ndcg_sum += ndcg_k(y_hat_i, y_i, y_i.shape[0])
        num += 1

    return ndcg_sum / num

# ---------------------------------------------------------------------------
# Training and Test Loop
# ---------------------------------------------------------------------------
best_rmse_ep, best_mae_ep, best_ndcg_ep = 0, 0, 0
best_rmse, best_mae, best_ndcg = float("inf"), float("inf"), 0

time_cumulative = 0
init = tf.global_variables_initializer()

with tf.Session(config=config) as sess:
    sess.run(init)

    # Pre-Training
    for i in range(epoch_p):
        tic = time()
        optimizer_p.minimize(sess, feed_dict={R: train_r})
        pre = sess.run(pred_p, feed_dict={R: train_r})

        t = time() - tic
        time_cumulative += t

        error = (test_m * (np.clip(pre, 1., 5.) - test_r) ** 2).sum() / test_m.sum()
        test_rmse = np.sqrt(error)

        error_train = (train_m * (np.clip(pre, 1., 5.) - train_r) ** 2).sum() / train_m.sum()
        train_rmse = np.sqrt(error_train)

        print('.-^-._' * 12)
        print('PRE-TRAINING')
        print('Epoch:', i+1, 'test rmse:', test_rmse, 'train rmse:', train_rmse)
        print('Time:', t, 'seconds')
        print('Time cumulative:', time_cumulative, 'seconds')
        print('.-^-._' * 12)

    # Fine-Tuning
    for i in range(epoch_f):
        tic = time()
        optimizer_f.minimize(sess, feed_dict={R: train_r})
        pre = sess.run(pred_f, feed_dict={R: train_r})

        t = time() - tic
        time_cumulative += t

        error = (test_m * (np.clip(pre, 1., 5.) - test_r) ** 2).sum() / test_m.sum()
        test_rmse = np.sqrt(error)

        error_train = (train_m * (np.clip(pre, 1., 5.) - train_r) ** 2).sum() / train_m.sum()
        train_rmse = np.sqrt(error_train)

        test_mae = (test_m * np.abs(np.clip(pre, 1., 5.) - test_r)).sum() / test_m.sum()
        train_mae = (train_m * np.abs(np.clip(pre, 1., 5.) - train_r)).sum() / train_m.sum()

        test_ndcg = call_ndcg(np.clip(pre, 1., 5.), test_r)
        train_ndcg = call_ndcg(np.clip(pre, 1., 5.), train_r)

        if test_rmse < best_rmse:
            best_rmse = test_rmse
            best_rmse_ep = i+1

        if test_mae < best_mae:
            best_mae = test_mae
            best_mae_ep = i+1

        if best_ndcg < test_ndcg:
            best_ndcg = test_ndcg
            best_ndcg_ep = i+1

        print('.-^-._' * 12)
        print('FINE-TUNING')
        print('Epoch:', i+1, 'test rmse:', test_rmse, 'test mae:', test_mae, 'test ndcg:', test_ndcg)
        print('Epoch:', i+1, 'train rmse:', train_rmse, 'train mae:', train_mae, 'train ndcg:', train_ndcg)
        print('Time:', t, 'seconds')
        print('Time cumulative:', time_cumulative, 'seconds')
        print('.-^-._' * 12)

    # Final result
    print('Epoch:', best_rmse_ep, ' best rmse:', best_rmse)
    print('Epoch:', best_mae_ep, ' best mae:', best_mae)
    print('Epoch:', best_ndcg_ep, ' best ndcg:', best_ndcg)

    # -----------------------------------------------------------------------
    # Final Evaluation (inside the same session to keep trained weights)
    # -----------------------------------------------------------------------
    pre_final = sess.run(pred_f, feed_dict={R: train_r})

# Train Evaluation
error_new = (train_m * (np.clip(pre_final, 1., 5.) - train_r)**2).sum() / train_m.sum()
rmse_new = np.sqrt(error_new)
preds = pre_final[train_m == 1]
trues = train_r[train_m == 1]
train_spearman = spearmanr(preds, trues)[0]

print('Train rmse:', rmse_new, 'Train spearman:', train_spearman)

# -------------------------------------------------------------------

# Test Evaluation
error_new_test = (final_test_m * (np.clip(pre_final, 1., 5.) - final_test_r)**2).sum() / final_test_m.sum()
rmse_new_test = np.sqrt(error_new_test)
preds_test = pre_final[final_test_m == 1]
trues_test = final_test_r[final_test_m == 1]
test_spearman = spearmanr(preds_test, trues_test)[0]

print('New test rmse:', rmse_new_test, 'New test spearman:', test_spearman)

# Matrix Rank
rank_m = np.linalg.matrix_rank(pre_final)
print("Rank of the constructed matrix is:", rank_m)

# Compute relative RMSE (rrmse = rmse / std of true values)
train_rrmse = rmse_new / np.std(trues)
test_rrmse = rmse_new_test / np.std(trues_test)

print('Train rrmse:', train_rrmse, 'Test rrmse:', test_rrmse)

# Create DataFrame
results_df = pd.DataFrame({
    'model': ["Glocal-K"],
    'time': [time_cumulative/60],
    'test_rmse': [rmse_new_test],
    'test_spearman': [test_spearman],
    'train_rmse': [rmse_new],
    'rank_m': [rank_m],
    'rank_beta': [np.nan],
    'rank_gamma': [np.nan],
    'sparse_beta': [np.nan],
    'sparse_gamma': [np.nan],
    'test_rrmse': [test_rrmse],
    'train_rrmse': [train_rrmse],
    'train_spearman': [train_spearman]
})

# Save to disk
out_path = './output/movielens/model_fits/glocalk_results_tensorflow.csv'
results_df.to_csv(out_path, index=False)
print("Results saved to {}".format(out_path))
