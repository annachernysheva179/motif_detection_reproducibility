import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    warnings.warn("deprecated", DeprecationWarning)
    print("stop warning")
import os
import sys
sys.path.insert(0,"vConv/vConvbaseddiscovery/code/")
from VConvMDcore import runvConvC
import tensorflow as tf
import pdb
from keras.backend.tensorflow_backend import set_session
import subprocess, re
import numpy as np
from build_models import *
from seq_to_matrix import *
import glob
from keras import backend as K
from datetime import datetime
def tictoc():
    return datetime.now().minute * 60 + datetime.now().second + datetime.now().microsecond*(10**-6)
import gc
from sklearn import mixture
import scipy.stats
from VConvMDcore import *

# added sequence padding option
class GeneRateOneHotMatrix():

    def __init__(self, padding=0):

        self.OriSeq = None
        self.Seqlabel = None
        self.SeqMatrix = None
        self.TrainX = None
        self.TestX = None
        self.TrainY = None
        self.TestY = None
        self.TrainID = None
        self.TestID = None
        self.Trainlist = []
        self.Testlist = []
        self.seq_pos_matrix_out = None
        self.padding = padding
        
    def mkdir(self, path):
        """
        :param path:
        :return:
        """
        isExists = os.path.exists(path)
        if not isExists:
            os.makedirs(path)
            return (True)
        else:
            return False

    def seq_to_matrix(self,seq, seq_matrix, seq_order):
        '''
        change target 3D tensor according to sequence and order
        :param seq: Input single root sequence
        :param seq_matrix: Input initialized matrix
        :param seq_order: This is the first sequence
        :return:
        '''
        for i in range(len(seq)):
            if ((seq[i] == 'A') | (seq[i] == 'a')):
                seq_matrix[seq_order, i, 0] = 1
            if ((seq[i] == 'C') | (seq[i] == 'c')):
                seq_matrix[seq_order, i, 1] = 1
            if ((seq[i] == 'G') | (seq[i] == 'g')):
                seq_matrix[seq_order, i, 2] = 1
            if ((seq[i] == 'T') | (seq[i] == 't')):
                seq_matrix[seq_order, i, 3] = 1
        for i in range(self.padding):
            seq_matrix[seq_order, len(seq)+i, 0] = 0.25
            seq_matrix[seq_order, len(seq)+i, 1] = 0.25
            seq_matrix[seq_order, len(seq)+i, 2] = 0.25
            seq_matrix[seq_order, len(seq)+i, 3] = 0.25
        return seq_matrix

    def genarate_matrix_for_train(self, seq_shape, seq_series):
        """
        generate a large matrix composed of sequences.
        :param shape: (seq number, sequence_length, 4)
        :param seq_series: A file in dataframe format composed of seq
        :return:seq
        """
        seq_matrix = np.zeros(seq_shape)
        for i in range(seq_series.shape[0]):
            seq_tem = seq_series[i]
            seq_matrix = self.seq_to_matrix(seq_tem, seq_matrix, i)
        return seq_matrix

    def cross_validation(self, number_of_folds, total_number, random_seeds=233):
        """
        This function is used to generate the index of n-fold cross validation
        :param number_of_folds:
        :param total_number:
        :param random_seeds:
        :return:
        """
        x = np.zeros((total_number,), dtype=np.int)
        split_iterator = StratifiedKFold(n_splits=number_of_folds, random_state=random_seeds, shuffle=True)
        split_train_index_and_test_index_list = [
            (train_index, test_index)
            for train_index, test_index in split_iterator.split(x, x)
            ]
        return (split_train_index_and_test_index_list)

    def split_dataset(self, split_index_list, fold, data_x, data_y, data_id=None):
        """
        Return training set and test set as needed
        :param split_index_list:Index of each fold returned by cross_validation
        :param fold:
        :param data_id:
        :param data_x:X
        :param data_y:Y
        :return:
        """
        id_train = data_id[split_index_list[fold][0].tolist()]
        x_train = data_x[split_index_list[fold][0].tolist()]
        y_train = data_y[split_index_list[fold][0].tolist()]
        id_test = data_id[split_index_list[fold][1].tolist()]
        x_test = data_x[split_index_list[fold][1].tolist()]
        y_test = data_y[split_index_list[fold][1].tolist()]
        return [x_train, y_train, id_train, x_test, y_test, id_test]

    def StoreTrainSet(self, rootPath):
        """
        Store multi-fold  data
        :param rootPath:
        :param ValNum: Generate valNum different validation data sets
        :param RandomSeeds:
        :param allData:
        :return:
        """
        self.mkdir(rootPath)
        training_path = rootPath + "/training_set.hdf5"
        test_path = rootPath + "/test_set.hdf5"

        f_train = h5py.File(training_path)
        f_test = h5py.File(test_path)

        f_train.create_dataset("sequences", data=self.TrainX)
        f_train.create_dataset("labs", data=self.TrainY)
        f_train.create_dataset("seq_idx", data=self.TrainID)
        f_train.close()
        f_test.create_dataset("sequences", data=self.TestX)
        f_test.create_dataset("labs", data=self.TestY)
        f_test.create_dataset("seq_idx", data=self.TestID)
        f_test.close()

    def k_mer_shuffle(self,seq_shape, seq_series, k=2):
        """
        Return a shuffled negative matrix
        :param seq_series:sequence list
        :param seq_shape:
        :param k:kshuffle
        :return:
        """
        seq_shuffle_matrix = np.zeros(seq_shape)

        for i in range(seq_shape[0]):
            seq = seq_series[i]
            shuffler = Shuffler(seq, k)
            seqres = shuffler.shuffle()
            seq_shuffle_matrix = self.seq_to_matrix(seqres, seq_shuffle_matrix, i)

        return seq_shuffle_matrix

    def GeneRateTrain(self, allData, ValNum=10, RandomSeeds=233):
        """

        """
        dataNum = allData[1].shape[0]
        split_train_index_and_test_index_list = self.cross_validation(number_of_folds=ValNum, total_number=dataNum,
                                                                 random_seeds=RandomSeeds)
        i = 0
        outDataTem = self.split_dataset(split_train_index_and_test_index_list, fold=i, data_x=allData[0], data_y=allData[1],
                                   data_id=allData[2])

        self.TrainX = outDataTem[0]
        self.TestX = outDataTem[3]
        self.TrainY = outDataTem[1]
        self.TestY = outDataTem[4]
        self.TrainID = outDataTem[2]
        self.TestID = outDataTem[5]

    def runSimple(self, SeqPath, OutputDir, SaveData=False):
        """
        Call this function to directly generate data, here is a separate data set, including Training and test
        seq format requirements, odd-numbered lines are names, even-numbered lines are corresponding sequences
        :return:
        """
        SeqLen = 0
        ChipSeqlFileFa = pd.read_csv(SeqPath, sep='\t', header=None, index_col=None)
        self.OriSeq = np.asarray(ChipSeqlFileFa[1::2]).reshape(np.asarray(ChipSeqlFileFa[1::2]).shape[0], )
        seq_positive_shape = self.OriSeq.shape[0]
        for i in range(seq_positive_shape):
            SeqLen = max(SeqLen, len(self.OriSeq[i]))
        seq_positive_name = np.asarray(ChipSeqlFileFa[::2]).reshape(seq_positive_shape, )
        self.seq_pos_matrix_out = self.genarate_matrix_for_train((seq_positive_shape, SeqLen+self.padding, 4),
                                                                 self.OriSeq)
        seq_pos_label_out = np.ones(seq_positive_shape, )

        seq_negative_name = seq_positive_name
        seq_neg_matrix_out = self.k_mer_shuffle((seq_positive_shape, SeqLen+self.padding, 4), self.OriSeq)
        seq_neg_label_out = np.zeros(seq_positive_shape, )

        seq = np.concatenate((self.seq_pos_matrix_out, seq_neg_matrix_out), axis=0)
        label = np.concatenate((seq_pos_label_out, seq_neg_label_out), axis=0)
        id_tem = np.concatenate((seq_positive_name, seq_negative_name), axis=0)
        index_shuffle = range(seq_positive_shape + seq_positive_shape)
        random.shuffle(index_shuffle)
        self.SeqMatrix = seq[index_shuffle, :, :]
        self.Seqlabel = label[index_shuffle]
        id_out = id_tem[index_shuffle].astype("string_")
        outData = [self.SeqMatrix, self.Seqlabel, id_out]
        self.GeneRateTrain(allData=outData, ValNum=10, RandomSeeds=233)
        if SaveData:
            self.StoreTrainSet(rootPath=OutputDir)


# added initial_lr option
def build_vCNN(model_template, number_of_kernel, max_kernel_length, lr=1,rho=0.99, epsilon=1.0e-8, input_shape=(1000,4)):
    """

    Building a vCNN model
    :param model: Input model
    :param number_of_kernel:number of kernel
    :param kernel_size: kernel size
    :param k_pool: Former K  maxpooling
    :param input_shape: Sequence shape
    :return:
    """
    model_template.add(VConv1D(
        input_shape=input_shape,
        kernel_size=(max_kernel_length),
        filters=number_of_kernel,
        padding='same',
        strides=1))
    model_template.add(keras.layers.pooling.MaxPooling1D(pool_length=10, stride=None, border_mode='valid'))
    model_template.add(keras.layers.GlobalMaxPooling1D())
    model_template.add(keras.layers.core.Dropout(0.1))
    model_template.add(keras.layers.core.Dense(output_dim=1))
    model_template.add(keras.layers.Activation("sigmoid"))
    sgd = keras.optimizers.Adadelta(lr=lr, rho=rho, epsilon=epsilon, decay=0.0)

    return model_template, sgd


def train_vCNN_hyperparameter_tuning(input_shape, modelsave_output_prefix, data_set, random_seed, batch_size, epoch_scheme, hyperparameters):
    def SelectBestModel(models):
        val = [float(name.split("_")[-1].split(".c")[0]) for name in models]
        index = np.argmin(np.asarray(val))
        return models[index]

    mkdir(modelsave_output_prefix)
    training_set, test_set = data_set
    X_train, Y_train = training_set
    X_test, Y_test = test_set
    tf.set_random_seed(random_seed)
    random.seed(random_seed)

    best_model, best_test_auc, best_modelsave_output_filename = None, -1, None
    
    for number_of_kernel in hyperparameters['number_of_kernel']:
        for max_ker_len in hyperparameters['max_ker_len']:
            for kernel_init_size in hyperparameters['kernel_init_size']:
                for lr in hyperparameters['lr']:
                    for rho in hyperparameters['rho']:
                        for epsilon in hyperparameters['epsilon']:
                            for mu in hyperparameters['mu']:

                                init_ker_len_dict = {str(kernel_init_size): number_of_kernel}    
                                
                                model = keras.models.Sequential()
                                model, sgd = build_vCNN(model, number_of_kernel, max_ker_len, input_shape=input_shape, lr=lr,rho=rho,epsilon=epsilon)
                                model = init_mask_final(model, init_ker_len_dict, max_ker_len)

                                modelsave_output_filename = modelsave_output_prefix + "/model_KernelNum-" + str(number_of_kernel) + "_initKernelLen-" + \
                                                            init_ker_len_dict.keys()[0]+ "_maxKernelLen-" + str(max_ker_len) + "_seed-" + str(random_seed) \
                                                            +"_lr-" + str(lr).replace(".","")\
                                                            +"_rho-" + str(rho).replace(".","") +"_epsilon-" + str(epsilon).replace("-","").replace(".","") \
                                                            +"_mu-" + str(mu).replace(".","")+ ".hdf5"

                                tmp_path = modelsave_output_filename.replace("hdf5", "pkl")
                                test_prediction_output = tmp_path.replace("/model_KernelNum-", "/Report_KernelNum-")

                                auc_records = []
                                loss_records = []
                                earlystopper = keras.callbacks.EarlyStopping(monitor='val_loss', patience=50, verbose=1)
                                tmp_hist = Histories(data = [X_test,Y_test])
                                reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=np.sqrt(0.1),
                                                            patience=20, min_lr=0.0001)
                                CrossTrain = TrainMethod()
                                KernelWeights, MaskWeight = model.layers[0].LossKernel, model.layers[0].MaskFinal
                                mu = K.cast(mu, dtype='float32')
                                lossFunction = ShanoyLoss(KernelWeights, MaskWeight, mu=mu)
                                model.compile(loss=lossFunction, optimizer=sgd, metrics=['accuracy'])
                                checkpointer = keras.callbacks.ModelCheckpoint(
                                    filepath=modelsave_output_filename.replace(".hdf5", ".checkpointer.hdf5"), verbose=1, save_best_only=True)

                                model.fit(X_train, Y_train, batch_size=batch_size, nb_epoch=int(epoch_scheme), shuffle=True,
                                validation_split=0.1,
                                verbose=2, callbacks=[checkpointer, reduce_lr, earlystopper, tmp_hist, CrossTrain])

                                auc_records.append(tmp_hist.aucs)
                                loss_records.append(tmp_hist.losses)
                                model.load_weights(modelsave_output_filename.replace(".hdf5", ".checkpointer.hdf5"))
                                y_pred = model.predict(X_test)
                                test_auc = roc_auc_score(Y_test, y_pred)
                                print("test_auc = {0}".format(test_auc))
                                best_auc = np.array([y for x in auc_records for y in x]).max()
                                best_loss = np.array([y for x in loss_records for y in x]).min()
                                print("best_auc = {0}".format(best_auc))
                                print("best_loss = {0}".format(best_loss))
                                report_dic = {}
                                report_dic["auc"] = auc_records
                                report_dic["loss"] = loss_records
                                report_dic["test_auc"] = test_auc

                                tmp_f = open(test_prediction_output, "wb")
                                pickle.dump(np.array(report_dic), tmp_f)
                                tmp_f.close()

                                if best_test_auc < test_auc:
                                    best_model, best_test_auc, best_modelsave_output_filename = model, test_auc, modelsave_output_filename

    return best_model, best_test_auc, best_modelsave_output_filename

# added top subsequences option
def KernelSeqDive_top(tmp_ker, seqs, Pos=True, top=500):
    """
    The kernel extracts the fragments on each sequence and the corresponding volume points.
     At the same time retain the position information on the sequence fragments mined by the kernel [sequence number, sequence start position, end position]
    :param tmp_ker:
    :param seqs:
    :return:
    """
    ker_len = tmp_ker.shape[0]
    inputs = K.placeholder(seqs.shape)
    ker = K.variable(tmp_ker.reshape(ker_len, 4, 1))
    conv_result = K.conv1d(inputs, ker, padding="valid", strides=1, data_format="channels_last")
    max_idxs = K.argmax(conv_result, axis=1)
    max_Value = K.max(conv_result, axis=1)
    sort_idxs = tf.nn.top_k(tf.transpose(max_Value,[1,0]), top, sorted=True).indices
    f = K.function(inputs=[inputs], outputs=[max_idxs, max_Value, sort_idxs])
    ret_idxs, ret, seq_top_idxs = f([seqs])
    if Pos:
        seqlist = []
        SeqInfo = []
        for seq_idx in seq_top_idxs[0]:
            start_idx = ret_idxs[seq_idx]
            seqlist.append(seqs[seq_idx, start_idx[0]:start_idx[0] + ker_len, :])
            SeqInfo.append([seq_idx, start_idx[0], start_idx[0] + ker_len])
        del f
        return seqlist, ret[seq_top_idxs], np.asarray(SeqInfo)
    else:
        return ret[seq_top_idxs]