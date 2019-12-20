import tensorflow as tf
import pandas as pd
import numpy as np
import time
import os


def get_next_batch(dataset1, batch_size_1, step, ind):
    start = step * batch_size_1
    end = ((step + 1) * batch_size_1)
    sel_ind = ind[start:end]
    newdataset1 = dataset1[sel_ind, :]
    return newdataset1


def train(drop_prob, dataset_train, dataset_test, normal_scale, sav=True, checkpoint_file='default.ckpt'):
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
        nz = dataset_test[:, :RNA_size].shape[0] * dataset_test[:, :RNA_size].shape[1]
        diff_mat = ((reconstruct - dataset_test[:, :RNA_size]) * normal_scale) ** 2
        loss_test_pretrain = np.sqrt(np.sum(diff_mat) / nz)
        print('RMSE loss at pretrain: ', loss_test_pretrain)

        ################## transfer learning -> in the target dataset
        dataset_train = target_data
        batch_size = 16
        num_epochs = 150

        dataset_size_train = dataset_train.shape[0]
        dataset_size_test = dataset_test.shape[0]
        print("Dataset size for training:", dataset_size_train)
        print("Dataset size for test:", dataset_size_test)
        num_iters = (num_epochs * dataset_size_train) // batch_size
        print("Num iters:", num_iters)
        ind_train = []
        for i in range(num_epochs):
            ind_train = np.append(ind_train, np.random.permutation(np.arange(dataset_size_train)))
        ind_train = np.asarray(ind_train).astype("int32")

        total_cost_train = 0.
        num_batchs = dataset_size_train // batch_size  # "//" => int()
        for step in range(num_iters):
            temp = get_next_batch(dataset_train, batch_size, step, ind_train)
            train_batch = np.asarray(temp).astype("float32")
            train_loss_val, _ = session.run([loss, optimizer],
                                            feed_dict={input_image: train_batch[:, RNA_size:],
                                                       original: train_batch[:, :RNA_size],
                                                       is_training: True})

            loss_val_list_train = np.append(loss_val_list_train, train_loss_val)
            total_cost_train += train_loss_val

            print_epochs = 10
            if step % (num_batchs * print_epochs) == 0:  # by epoch, num_batchs * batch_size = dataset_train_size
                ############### test section  ###########################
                dataset_test = np.asarray(dataset_test).astype("float32")
                test_loss_val = session.run(loss, feed_dict={input_image: dataset_test[:, RNA_size:],
                                                             original: dataset_test[:, :RNA_size],
                                                             is_training: False})
                loss_val_list_test = np.append(loss_val_list_test, test_loss_val)

                reconstruct = session.run(reconstructed_image,
                                          feed_dict={input_image: dataset_test[:, RNA_size:], is_training: False})
                nz = dataset_test[:, :RNA_size].shape[0] * dataset_test[:, :RNA_size].shape[1]
                diff_mat = ((reconstruct - dataset_test[:, :RNA_size]) * normal_scale) ** 2
                loss_test = np.sqrt(np.sum(diff_mat) / nz)

                print('RMSE loss by train_data_size: ', step, "/", num_iters,
                      total_cost_train / (num_batchs * print_epochs),
                      test_loss_val, loss_test)
                # print('RMSE loss by train_data_size: ', step, "/", num_iters, total_cost_train / num_batchs,
                #       total_cost_validation / num_batchs)
                #                 print('new loss: ', step, "/", num_iters, train_loss_val, valid_loss_val)
                total_cost_train = 0.
                total_cost_validation = 0.

        ####final test
        ############### test section  ###########################
        dataset_test = np.asarray(dataset_test).astype("float32")

        test_loss_val = session.run(loss, feed_dict={input_image: dataset_test[:, RNA_size:],
                                                     original: dataset_test[:, :RNA_size],
                                                     is_training: False})
        loss_val_list_test = np.append(loss_val_list_test, test_loss_val)

        reconstruct = session.run(reconstructed_image,
                                  feed_dict={input_image: dataset_test[:, RNA_size:], is_training: False})
        nz = dataset_test[:, :RNA_size].shape[0] * dataset_test[:, :RNA_size].shape[1]
        diff_mat = ((reconstruct - dataset_test[:, :RNA_size]) * normal_scale) ** 2
        loss_test = np.sqrt(np.sum(diff_mat) / nz)

        print('RMSE loss by train_data_size: ', step, "/", num_iters, total_cost_train / (num_batchs * print_epochs),
              test_loss_val, loss_test)

    end = time.time()
    el = end - start
    print(("Time elapsed %f" % el))
    return (loss_val_list_train, loss_val_list_test, loss_test, loss_test_pretrain, reconstruct)

#############################################################################################################
import getopt
import sys
print(sys.argv)

os.environ["CUDA_VISIBLE_DEVICES"] = '0' #sys.argv[1]
original_dat_path_DNA = '/data0/zhoux/DNA_WT.csv'
original_dat_path_RNA = '/data0/zhoux/UCSC_RNA_WT.csv'
imputed_dataset_path = '/data0/zhoux/imputed_RNA_WT.csv'
DNA_target = pd.read_csv(original_dat_path_DNA, delimiter=',',index_col=0, header=0)
RNA_target = pd.read_csv(original_dat_path_RNA, delimiter=',', index_col=0, header=0)
DNA_target.index = [x[:19] for x in DNA_target.index.values]
RNA_target.index = [x[:19] for x in RNA_target.index.values]
cancer_target = pd.merge(RNA_target, DNA_target, left_index=True, right_index=True, how = 'inner')

RNA_size = RNA_target.shape[1]
DNA_size = DNA_target.shape[1]
normal_scale = np.max(np.max(cancer_target.iloc[:, :RNA_size]))+0.001

aa = np.concatenate((cancer_target.values[:, :RNA_size] / normal_scale, cancer_target.values[:, RNA_size:]), axis=1)
cancer_target = pd.DataFrame(aa, index=cancer_target.index, columns=cancer_target.columns)


## train/test data split
missing_perc = 0.5
train_data = cancer_target.sample(frac = (1-missing_perc), random_state=sample_count, axis=0, replace=False)
test_data = cancer_target[~cancer_target.index.isin(train_data.index)]
new_dataset = pd.concat([test_data, train_data], axis=0)
train_data = train_data.values
test_data = test_data.values
print('train datasize:', train_data.shape[0], ' test datasize: ', test_data.shape[0])

lr = 0.0001
drop_prob = 0.
batch_shape_input = (None, DNA_size)
batch_shape_output = (None, RNA_size)
tf.reset_default_graph()
loss_val_list_train, loss_val_list_test, loss_test, loss_test_pretrain, reconstruct = train(drop_prob,
                                                                                            train_data,
                                                                                            test_data,
                                                                                            normal_scale,
                                                                                            sav=save_ckpt,
                                                                                            checkpoint_file="./checkpoints/ref_general_model_quantiles.ckpt")

imputed_data = np.concatenate([reconstruct * normal_scale, train_data[:, :RNA_size] * normal_scale], axis=0)
RNA_txt = pd.DataFrame(imputed_data[:, :RNA_size], index=new_dataset.index, columns=new_dataset.columns[:RNA_size])
RNA_txt.to_csv(imputed_dataset_path+'/TDimpute_WT_'+str(missing_perc*100)+'.csv')