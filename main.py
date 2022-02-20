#from predict.ensemble import ensemble
from numpy import extract
from utils.argparser import argparser
import os


from data_processing.generate_dataset import extract_dataset
from predict.predict import predict
from predict.visgrid_pred import visgrid_pred
from data_processing.generate_dataset import extract_dataset

from predict.ensemble import ensemble

from predict.rank_pockets import ranking


if __name__ == "__main__":
    params = argparser()
    os.environ['CUDA_VISIBLE_DEVICES'] = params['gpu']

    if params['mode']==1: #
        input_dir = params['test_dir']
        output_dir = params['test_odir']
        #output_dir = "../dataset/test_processed/"
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
            extract_dataset(input_dir,output_dir)
        predict(params)
        ensemble(params)
        bad_case = ranking()
        if len(bad_case) > 0:
            visgrid_pred(bad_case,params)