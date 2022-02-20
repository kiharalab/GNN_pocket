import argparse



def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode',type=int,default=1,help='1: predicting')
    parser.add_argument('--test_dir',type=str,default="./dataset/test")
    parser.add_argument('--test_odir',type=str,default="./dataset/test_processed")
    parser.add_argument('--gpu',type=str,default='0',help='Choose gpu id, example: \'1,2\'(specify use gpu 1 and 2)')
    args = parser.parse_args()
    params = vars(args)
    return params
