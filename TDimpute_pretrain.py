import tensorflow as tf
import pandas as pd
import numpy as np
import time
import os

def train(drop_prob, dataset_test, normal_scale, sav=True, checkpoint_file='default.ckpt'):
    input_image = tf.placeholder(tf.float32, batch_shape_input, name='input_image')
    is_training = tf.placeholder(tf.bool)
    scale = 0.
    with tf.variable_scope('FCN') as scope:
        fc_1 = tf.layers.dense(inputs=input_image, units=4000,
                               kernel_regularizer=tf.contrib.layers.l2_regularizer(scale=scale))
        fc_1_out = tf.nn.sigmoid(fc_1)
        fc_1_dropout = tf.layers.dropout(inputs=fc_1_out, rate=drop_prob, training=is_training)
        fc_2_dropout = tf.layers.dense(inputs=fc_1_dropout, units=RNA_size)  # 46744
        fc_2_out = tf.nn.sigmoid(fc_2_dropout)  # fc_2_dropout #
        reconstructed_image = fc_2_out  # fc_2_dropout

    original = tf.placeholder(tf.float32, batch_shape_output, name='original')
    loss = tf.sqrt(tf.reduce_mean(tf.square(tf.subtract(reconstructed_image, original))))
    l2_loss = tf.losses.get_regularization_loss()
    optimizer = tf.train.AdamOptimizer(lr).minimize(loss + l2_loss)

    init = tf.global_variables_initializer()
    saver = tf.train.Saver()
    start = time.time()
    loss_val_list_train = 0
    loss_val_list_test = 0
    loss_test = 0

    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    with tf.Session(config=config) as session:
        session.run(init)
        print(("Loading variables from '%s'." % checkpoint_file))
        saver.restore(session, checkpoint_file)
        print('restored')

        ############### test the pretrain model for target dataset
        dataset_test = np.asarray(dataset_test).astype("float32")
        reconstruct = session.run(reconstructed_image,
                                  feed_dict={input_image: dataset_test[:, RNA_size:], is_training: False})

    end = time.time()
    el = end - start
    print(("Time elapsed %f" % el))
    return (loss_val_list_train, loss_val_list_test, loss_test, loss_test_pretrain, reconstruct)

#############################################################################################################
os.environ["CUDA_VISIBLE_DEVICES"] = '0'
original_dat_path_DNA = '/data0/zhoux/DNA_WT.csv'
imputed_dataset_path = '/data0/zhoux/imputed_RNA_WT.csv'

DNA_target = pd.read_csv(original_dat_path_DNA, delimiter=',',index_col=0, header=0)
DNA_target.index = [x[:19] for x in DNA_target.index.values]
DNA_size = DNA_target.shape[1]
save_ckpt = True
test_data = DNA_target.values

lr = 0.0001
drop_prob = 0.
batch_shape_input = (None, DNA_size)
batch_shape_output = (None, RNA_size)
tf.reset_default_graph()
loss_val_list_train, loss_val_list_test, loss_test, loss_test_pretrain, reconstruct = train(drop_prob,
                                                                                            test_data,
                                                                                            sav=save_ckpt,
                                                                                            checkpoint_file="./checkpoints/ref_general_model_quantiles.ckpt")


## only use the gene name list
original_dat_path_RNA = '/data0/zhoux/UCSC_RNA_WT.csv'
RNA_target = pd.read_csv(original_dat_path_RNA, delimiter=',', index_col=0, header=0)

RNA_txt = pd.DataFrame(reconstruct, index=DNA_target.index, columns=RNA_target.columns)
RNA_txt.to_csv(imputed_dataset_path)