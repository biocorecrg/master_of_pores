#!/usr/bin/env python3
# coding: utf-8

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
import argparse
import sys
import os
from copy import deepcopy
import re
import csv
import time
import configparser
import h5py
import traceback
import math
import numpy as np
# from PIL import Image
import pyts
from pyts.image import MarkovTransitionField, GramianAngularField, RecurrencePlot
import tensorflow as tf
import keras
from keras.layers import Dense, Conv2D, BatchNormalization, Activation
from keras.layers import AveragePooling2D, Input, Flatten
from keras.optimizers import Adam
from keras.callbacks import ModelCheckpoint, LearningRateScheduler
from keras.callbacks import ReduceLROnPlateau
from keras.preprocessing.image import ImageDataGenerator
from keras.utils import multi_gpu_model
from keras.regularizers import l2
from keras import backend as K
from keras.models import Model
import pandas as pd
from sklearn import datasets, linear_model
from sklearn.model_selection import train_test_split
from tensorflow.python.client import device_lib
from keras.models import load_model



'''

    James M. Ferguson (j.ferguson@garvan.org.au)
    Genomic Technologies
    Garvan Institute
    Copyright 2019

    Tansel Ersevas (t.ersevas@garvan.org.au)

    script description



    ----------------------------------------------------------------------------
    version 0.0.0 - initial
    version 0.8.0 - CPU version Done
    version 0.9.0 - Fixed segment offset
    version 0.9.1 - added segment and squiggle output
    version 0.9.2 - separate segment output and code clean up
    version 1.0.0 - initial release

    So a cutoff of: 0.4958776 for high accuracy
    and another of 0.2943664 for high recovery

    TODO:
        - Remove leftover libraries
        - remove debug plots
        - Remove redundant code
        - create log files with information
        - take in fastq for dmux splitting
        - take in paf or bam for training splitting


    ----------------------------------------------------------------------------
    MIT License

    Copyright (c) 2019 James M. Ferguson

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
'''




class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def print_verbose(message):
    '''verbose printing'''
    sys.stderr.write('info: %s\n' % message)

def print_err(message):
    '''error printing'''
    sys.stderr.write('error: %s\n' % message)

def read_config(filename):
    config = configparser.ConfigParser()
    config.read(filename)
    return(config)

def _get_available_devices():
    local_device_protos = device_lib.list_local_devices()
    return [x.name for x in local_device_protos]

def _check_available_devices():
    available_devices = _get_available_devices()
    print_verbose(available_devices)
    # Make sure requested GPUs are available or at least warn if they aren't
    return(TRUE)

def read_model(model_name):
    # model = load_model('saved_models/' + model_name)
    model = load_model(model_name) # as a path
    model.compile(loss='categorical_crossentropy',
                     optimizer=Adam(),
                     metrics=['accuracy'])
    return(model)

squiggle_max = 1199
squiggle_min = 1
input_cut = 72000 #potenitall need to be configurable
image_size = 224
num_classes = 4
window = 2000

def main():
    '''
    Main function
    '''
    VERSION = "1.0.0"

    parser = MyParser(
        description="DeePlexiCon - Demultiplex direct RNA reads")
    #group = parser.add_mutually_exclusive_group()
    parser.add_argument("-p", "--path",
                        help="Top path of fast5 files to dmux")
    parser.add_argument("-f", "--form", default="multi", choices=["multi", "single"],
                        help="Multi or single fast5s")
    parser.add_argument("-c", "--config",
                        help="config file")
    # parser.add_argument("-g", "--gpu_list", default=1
    #                     help="list of gpus, 1, or [1,3,5], etc. of PCI_BUS_ID order")
    # parser.add_argument("-o", "--output",
    #                     help="Output directory")
    parser.add_argument("-s", "--threshold", type=float, default=0.50,
                        help="probability threshold - 0.5 hi accuracy / 0.3 hi recovery")
                        # populate choices with models found in saved_models/
    # parser.add_argument("-m", "--model", default="4_bc_normal.h5", choices=["4_bc_normal.h5", "model2"],
    parser.add_argument("-m", "--model",
                        help="Trained model name to use")
    parser.add_argument("--squiggle",
                        help="dump squiggle data into this .tsv file")
    parser.add_argument("--segment",
                        help="dump segment data into this .tsv file")
    parser.add_argument("-b", "--batch_size", type=int, default=4000,
                        help="batch size - for single fast5s")
    parser.add_argument("-V", "--version",
                        help="Prints version")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()
    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.verbose:
        print_verbose("Verbose mode active - dumping info to stderr")
        print_verbose("DeePlexiCon: {}".format(VERSION))
        print_verbose("arg list: {}".format(args))
        if tf.test.gpu_device_name():
            print_verbose('Default GPU Device: {}'.format(tf.test.gpu_device_name()))
        else:
            print_verbose("Please install GPU version of TF")

    if args.squiggle:
        squig_file = args.squiggle
        with open(squig_file, 'a') as f:
            f.write("{}\t{}\n".format("ReadID", "signal_pA"))
    else:
        squig_file = ''

    if args.segment:
        seg_file = args.segment
        with open(seg_file, 'a') as f:
            f.write("{}\t{}\t{}\n".format("ReadID", "start", "stop"))
    else:
        seg_file = ''

    # Globals
    if args.config:
        config = read_config(args.config) #TODO check config read error

    # gpu settings
    # Devices
    # os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"

    # os.environ["CUDA_VISIBLE_DEVICES"] = config[deeplexicon][gpu_list] if args.config else args.gpu_list

    # do check devices are available, else throw and error


    # main logic

    # read model

    model = read_model(config[deeplexicon][trained_model]) if args.config else read_model(args.model)
    barcode_out = {0: "bc_1",
                   1: "bc_2",
                   2: "bc_3",
                   3: "bc_4",
                   None: "unknown"
                   }
    labels = []
    images = []
    fast5s = {}
    stats = ""
    seg_dic = {}
    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format("fast5", "ReadID", "Barcode", "Confidence Interval", "P_bc_1", "P_bc_2", "P_bc_3", "P_bc_4"))
    # for file in input...
    for dirpath, dirnames, files in os.walk(args.path):
        for fast5 in files:
            if fast5.endswith('.fast5'):
                fast5_file = os.path.join(dirpath, fast5)
                if args.form == "single":
                    #everthing below this, send off in batches of N=args.batch_size
                    # The signal extraction and segmentation can happen in the first step
                    # read fast5 files
                    readID, seg_signal = get_single_fast5_signal(fast5_file, window, squig_file, seg_file)
                    if not seg_signal:
                        print_err("Segment not found for:\t{}\t{}".format(fast5_file, readID))
                        continue
                    # convert
                    sig = np.array(seg_signal, dtype=float)
                    img = convert_to_image(sig)
                    labels.append(readID)
                    fast5s[readID] = fast5
                    images.append(img)
                    # classify
                    if len(labels) >= args.batch_size:
                        C = classify(model, labels, np.array(images), False, args.threshold)
                        # save to output
                        for readID, out, c, P in C:
                            prob = [round(float(i), 6) for i in P]
                            cm = round(float(c), 4)
                            if args.verbose:
                                print_verbose("cm is: {}".format(cm))
                            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(fast5s[readID], readID, barcode_out[out], cm, prob[0], prob[1], prob[2], prob[3]))
                        labels = []
                        images = []
                        fast5s = {}


                elif args.form == "multi":
                    #everthing below this, send off in batches of N=args.batch_size
                    # The signal extraction and segmentation can happen in the first step
                    # read fast5 files
                    seg_signal = get_multi_fast5_signal(fast5_file, window, squig_file, seg_file)
                    sig_count = 0
                    for readID in seg_signal:
                        # convert
                        img = convert_to_image(np.array(seg_signal[readID], dtype=float))
                        labels.append(readID)
                        images.append(img)
                        fast5s[readID] = fast5
                        sig_count += 1
                        if len(labels) >= args.batch_size:
                            C = classify(model, labels, np.array(images), False, args.threshold)
                            # save to output
                            for readID, out, c, P in C:
                                prob = [round(float(i), 6) for i in P]
                                cm = round(float(c), 4)
                                if args.verbose:
                                    print_verbose("cm is: {}".format(cm))
                                print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(fast5s[readID], readID, barcde_out[out], cm, prob[0], prob[1], prob[2], prob[3]))
                            labels = []
                            images = []
                            fast5s = {}
                        elif args.verbose:
                            print_verbose("analysing sig_count: {}/{}".format(sig_count, len(seg_signal)))
                        else:
                            blah = 0 # clean
    #finish up
    C = classify(model, labels, np.array(images), False, args.threshold)
    # save to output
    for readID, out, c, P in C:
        prob = [round(float(i), 6) for i in P]
        cm = round(float(c), 4)
        if args.verbose:
            print_verbose("cm is: {}".format(cm))
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(fast5s[readID], readID, barcode_out[out], cm, prob[0], prob[1], prob[2], prob[3]))
    labels = []
    images = []
    fast5s = {}


    # final report/stats
    # print stats

# file handling and segmentation

def get_single_fast5_signal(read_filename, w, squig_file, seg_file):
    '''
    open sigle fast5 file and extract information
    '''

    # get readID and signal
    f5_dic = read_single_fast5(read_filename)
    if not f5_dic:
        print_err("Signal not extracted from: {}".format(read_filename))
        return 0, 0
    # segment on raw
    readID = f5_dic['readID']
    signal = f5_dic['signal']
    seg = dRNA_segmenter(readID, signal, w)
    if not seg:
        print_verbose("No segment found - skipping: {}".format(readID))
        return 0, 0
    # convert to pA
    pA_signal = convert_to_pA(f5_dic)
    if squig_file:
        with open(squig_file, 'a') as f:
            f.write("{}\t{}\n".format(readID, "\t".join(pA_signal)))
    if seg_file:
        with open(seg_file, 'a') as f:
            f.write("{}\t{}\t{}\n".format(readID, seg[0], seg[1]))
    # return signal/signals
    return readID, pA_signal[seg[0]:seg[1]]


def get_multi_fast5_signal(read_filename, w, squig_file, seg_file):
    '''
    open multi fast5 files and extract information
    '''
    pA_signals = {}
    seg_dic = {}
    f5_dic = read_multi_fast5(read_filename)
    seg = 0
    sig_count = 0
    for read in f5_dic:
        sig_count += 1
        print_verbose("reading sig_count: {}/{}".format(sig_count, len(f5_dic)))
        # get readID and signal
        readID = f5_dic[read]['readID']
        signal = f5_dic[read]['signal']

        # segment on raw
        seg = dRNA_segmenter(readID, signal, w)
        if not seg:
            seg = 0
            continue
        # convert to pA
        pA_signal = convert_to_pA(f5_dic[read])
        if squig_file:
            with open(squig_file, 'a') as f:
                f.write("{}\t{}\n".format(readID, "\t".join(pA_signal)))
        if seg_file:
            with open(seg_file, 'a') as f:
                f.write("{}\t{}\t{}\n".format(readID, seg[0], seg[1]))
        pA_signals[readID] = pA_signal[seg[0]:seg[1]]
        seg_dic[readID] = seg
    # return signal/signals
    return pA_signals


def read_single_fast5(filename):
    '''
    read single fast5 file and return data
    '''
    f5_dic = {'signal': [], 'readID': '', 'digitisation': 0.0,
              'offset': 0.0, 'range': 0.0, 'sampling_rate': 0.0}

    # open fast5 file
    try:
        hdf = h5py.File(filename, 'r')
    except:
        traceback.print_exc()
        print_err("extract_fast5():fast5 file failed to open: {}".format(filename))
        f5_dic = {}
        return f5_dic
    try:
        c = list(hdf['Raw/Reads'].keys())
        for col in hdf['Raw/Reads/'][c[0]]['Signal'][()]:
            f5_dic['signal'].append(int(col))

        f5_dic['readID'] = hdf['Raw/Reads/'][c[0]].attrs['read_id'].decode()
        f5_dic['digitisation'] = hdf['UniqueGlobalKey/channel_id'].attrs['digitisation']
        f5_dic['offset'] = hdf['UniqueGlobalKey/channel_id'].attrs['offset']
        f5_dic['range'] = float("{0:.2f}".format(hdf['UniqueGlobalKey/channel_id'].attrs['range']))
        f5_dic['sampling_rate'] = hdf['UniqueGlobalKey/channel_id'].attrs['sampling_rate']

    except:
        traceback.print_exc()
        print_err("extract_fast5():failed to extract events or fastq from: {}".format(filename))
        f5_dic = {}

    return f5_dic


def read_multi_fast5(filename):
    '''
    read multifast5 file and return data
    '''
    f5_dic = {}
    with h5py.File(filename, 'r') as hdf:
        for read in list(hdf.keys()):
            f5_dic[read] = {'signal': [], 'readID': '', 'digitisation': 0.0,
                            'offset': 0.0, 'range': 0.0, 'sampling_rate': 0.0}
            try:
                for col in hdf[read]['Raw/Signal'][()]:
                    f5_dic[read]['signal'].append(int(col))

                f5_dic[read]['readID'] = hdf[read]['Raw'].attrs['read_id'].decode()
                f5_dic[read]['digitisation'] = hdf[read]['channel_id'].attrs['digitisation']
                f5_dic[read]['offset'] = hdf[read]['channel_id'].attrs['offset']
                f5_dic[read]['range'] = float("{0:.2f}".format(hdf[read]['channel_id'].attrs['range']))
                f5_dic[read]['sampling_rate'] = hdf[read]['channel_id'].attrs['sampling_rate']
            except:
                traceback.print_exc()
                print_err("extract_fast5():failed to read readID: {}".format(read))
    return f5_dic


def dRNA_segmenter(readID, signal, w):
    '''
    segment signal/s and return coords of cuts
    '''
    def _scale_outliers(squig):
        ''' Scale outliers to within m stdevs of median '''
        k = (squig > 0) & (squig < 1200)
        return squig[k]


    sig = _scale_outliers(np.array(signal, dtype=int))

    s = pd.Series(sig)
    t = s.rolling(window=w).mean()
    # This should be done better, or changed to median and benchmarked
    # Currently trained on mean segmented data
    mn = t.mean()
    std = t.std()
    # Trained on 0.5
    bot = mn - (std*0.5)

    # main algo

    begin = False
    # max distance for merging 2 segs
    seg_dist = 1500
    # max length of a seg
    hi_thresh = 200000
    # min length of a seg
    lo_thresh = 2000

    start = 0
    end = 0
    segs = []
    count = -1
    for i in t:
        count += 1
        if i < bot and not begin:
            start = count
            begin = True
        elif i < bot:
            end = count
        elif i > bot and begin:
            if segs and start - segs[-1][1] < seg_dist:
                segs[-1][1] = end
            else:
                segs.append([start, end])
            start = 0
            end = 0
            begin = False
        else:
            continue

    # offset = -1050
    # buff = 150
    offset = -1000
    buff = 0

    x, y = 0, 0

    for a, b in segs:
        if b - a > hi_thresh:
            continue
        if b - a < lo_thresh:
            continue
        x, y = a, b

        # to be modified in next major re-training
        return [x+offset-buff, y+offset+buff]
        break
    print_verbose("dRNA_segmenter: no seg found: {}".format(readID))
    return 0


def convert_to_pA(d):
    '''
    convert raw signal data to pA using digitisation, offset, and range
    float raw_unit = range / digitisation;
    for (int32_t j = 0; j < nsample; j++) {
        rawptr[j] = (rawptr[j] + offset) * raw_unit;
    }
    '''
    digitisation = d['digitisation']
    range = d['range']
    offset = d['offset']
    raw_unit = range / digitisation
    new_raw = []
    for i in d['signal']:
        j = (i + offset) * raw_unit
        new_raw.append("{0:.2f}".format(round(j,2)))
    return new_raw


def pyts_transform(transform, data, image_size, show=False, cmap='rainbow', img_index=0):
    try:
        t_start=time.time()
        X_transform = transform.fit_transform(data)
        if (show):
            plt.figure(figsize=(4, 4))
            plt.grid(b=None)
            plt.imshow(X_transform[0], cmap=cmap, origin='lmtfower')
            plt.savefig(transform.__class__.__name__ + "_image_" + str(img_index) + ".svg", format="svg")
            plt.show()
        return(X_transform)
    except Exception as e:
        print_err(str(e))
        return([])


def mtf_transform(data, image_size=500, show=False, img_index=0):
    transform = MarkovTransitionField(image_size)
    return(pyts_transform(transform, data, image_size=image_size, show=show, cmap='rainbow', img_index=img_index))

def rp_transform(data, image_size=500 ,show=False ,img_index=0):
    # RP transformationmtf
    transform = RecurrencePlot(dimension=1,
                    threshold='percentage_points',
                    percentage=30)
    return(pyts_transform(transform, data, image_size=image_size, show=show, cmap='binary', img_index=img_index))

def gasf_transform(data, image_size=500, show=False, img_index=0):
    # GAF transformation
    transform = GramianAngularField(image_size, method='summation')
    return(pyts_transform(transform, data, image_size=image_size, show=show, cmap='rainbow', img_index=img_index))

def gadf_transform(data, image_size=500, show=False ,img_index=0):
    # GAF transformation
    transform = GramianAngularField(image_size, method='difference')
    return(pyts_transform(transform, data, image_size=image_size, show=show, cmap='rainbow', img_index=img_index))


def labels_for(a_file_name):
    segments=re.split(r'[_\-\.]+', a_file_name)
    return(segments)

def max_in_sequence(sequence):
    return(max(np.amax([list(d.values()) for d in sequence]), 0.01))

def compress_squiggle(squiggle, compress_factor):
    squiggle_len = len(squiggle)
    rem = squiggle_len % compress_factor
    if rem > 0:
        return(np.mean(squiggle[0:squiggle_len - rem].reshape(-1,compress_factor), axis=1))
    return(squiggle)

def convert_to_image(signal):
   transformed_squiggle = gasf_transform(signal.reshape(1,-1), image_size=image_size, show=False)
   return(transformed_squiggle)


def confidence_margin(npa):
    sorted = np.sort(npa)[::-1]    #return sort in reverse, i.e. descending
    # sorted = np.sort(npa)   #return sort in reverse, i.e. descending
    d = sorted[0] - sorted[1]
    return(d)

def classify(model, labels, image, subtract_pixel_mean, threshold):
    input_shape = image.shape[1:]
    # x = image.astype('float32') / 255
    x = image.astype('float32') + 1
    x = x / 2

    # If subtract pixel mean is enabled
    if subtract_pixel_mean:
        x_mean = np.mean(x, axis=0)
        x -= x_mean
    x=[x]
    y = model.predict(x, verbose=0)
    res = []
    for i in range(len(y)):
        cm = confidence_margin(y[i])
        if y[i][np.argmax(y[i])] >= threshold:
            res.append([labels[i], np.argmax(y[i]), cm, y[i]])
        else:
            res.append([labels[i], None, cm, y[i]])
    return res

if __name__ == '__main__':
    main()
