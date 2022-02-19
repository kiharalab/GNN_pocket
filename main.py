#from predict.ensemble import ensemble
from utils.argparser import argparser
import os
# from data_processing import extract_edges, extract_face_feature,extract_label,extract_pocket_feature,extract_vertex_feature
# from data_processing import prepare_input
# from train import train

# from data_processing.extract_edges import extract_edges
# from data_processing.extract_face_feature import extract_face_feature
# from data_processing.extract_label import extract_label
# from data_processing.extract_pocket_feature import extract_pocket_feature
# from data_processing.extract_vertex_feature import extract_vertex_feature
# from data_processing.extract_atom_type import extract_atoms
# from data_processing.prepare_input import prepare_input

from data_processing.generate_dataset import extract_dataset
from train.train import train
from predict.predict import predict
from predict.visgrid_pred import visgrid_pred

from predict.ensemble import ensemble

from predict.rank_pockets import ranking


if __name__ == "__main__":
    params = argparser()
    os.environ['CUDA_VISIBLE_DEVICES'] = params['gpu']
    if params['mode']==0: #train
        input_dir = params['train_dir']
        output_dir = params['train_odir']
        # input_dir = "../dataset/train/"
        # output_dir = "../dataset/train_processed/"
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
            extract_dataset(input_dir,output_dir)
        
        train(params)
    elif params['mode']==1: #
        input_dir = params['test_dir']
        output_dir = params['test_odir']
        #output_dir = "../dataset/test_processed/"
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
            extract_dataset(params,input_dir,output_dir)
        #predict(params)
        ensemble(params)
        bad_case = ranking()
        print(len(bad_case))
        if len(bad_case) > 0:
            visgrid_pred(bad_case,params)