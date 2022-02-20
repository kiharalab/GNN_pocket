import os
import numpy as np

import numpy as np
from sklearn.cluster import AgglomerativeClustering



def ranking():
    cluster_model = AgglomerativeClustering(n_clusters=None,distance_threshold=12,linkage='single')
    label_dir = '../dataset/test_processed/label'
    input_dir = "../dataset/test/"
    pred_dir = '../test_prediction/avg_pdb_pred'
    dirs = os.listdir(label_dir)
    #dirs = ['949.npy','973.npy','994.npy']
    bad_case = []
    for id in dirs:
        output_dir = os.path.join("../final_pred",id[:-4])
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        pred_path = os.path.join(pred_dir,id[:-4]+'.pdb')
        with open(pred_path,'r') as file:
            lines = file.readlines()
        if len(lines)==0:
            bad_case.append(id[:-4])
            print(id,0)
            input_path = os.path.join(input_dir,id[:-4]+'/structure.pqr')
            output_path = os.path.join(output_dir,'structure.pqr')
            with open(input_path,'r') as file:
                rows = file.read()
            with open(output_path,'w') as f:
                f.write(rows)
            continue
        elif len(lines)==1:
            bad_case.append(id[:-4])
            print(id,1)
            input_path = os.path.join(input_dir,id[:-4]+'/structure.pqr')
            output_path = os.path.join(output_dir,'structure.pqr')
            atom_id = lines[0][:11]
            new_lines = ''
            with open(input_path,'r') as file:
                rows = file.readlines()
            for row in rows:
                if atom_id in row:
                    new_row = row[:59]+'1'+row[60:]
                else:
                    new_row =row
                new_lines += new_row
            with open(output_path,'w') as f:
                f.write(new_lines)
            continue
        label_path = os.path.join(label_dir,id)
        label = np.load(label_path)
        coords = np.zeros((len(lines),3))
        prob = []
        atoms = []
        for i,line in enumerate(lines):
            line = line.strip('\n').split(' ')
            line  = [i for i in line if len(i)!=0]
            atoms.append(int(line[1]))
            coords[i][0] = float(line[5])
            coords[i][1] = float(line[6])
            coords[i][2] = float(line[7])
            prob.append(float(line[10]))
        clustering = cluster_model.fit_predict(coords)
        cluster_label = cluster_model.labels_
        #print(cluster_model.labels_)
        cluster_label_set = set(cluster_label)

        cluster_score = {}
        for i in cluster_label_set:
            score = 0
            index_i = [j for j,x in enumerate(list(cluster_label)) if x == i ]
            for idx in index_i:
                score += prob[idx]
            cluster_score[i] = score
        cluster_score = sorted(cluster_score.items(), key=lambda x: x[1], reverse=True)

        if len(cluster_label)>10:
            cluster_score = cluster_score[:10]
        ranked_cluster_score = [i[0] for i in cluster_score]
        TP_list = np.zeros(len(label))
        for idx,rank in enumerate(ranked_cluster_score):
            idx_lst = [(j,idx+1) for j,x in enumerate(list(cluster_label)) if x == rank ]
            for item in idx_lst:
                TP_list[atoms[item[0]]-1] = item[1]
        num = len(np.argwhere(TP_list>0))
        if num < 10:
            bad_case.append(id[:-4])
        output_lines = ''
        input_path = os.path.join(input_dir,id[:-4]+'/structure.pqr')
        output_path = os.path.join(output_dir,'structure.pqr')

        with open(input_path,'r') as file:
            rows = file.readlines()
        for idx,row in enumerate(rows):
            new_row = row
            if TP_list[idx]>0 and TP_list[idx]<10:
                new_row = new_row[:59]+str(TP_list[idx])[0]+new_row[60:]
            elif TP_list[idx]==10:
                new_row = new_row[:58]+str(TP_list[idx])[:2]+new_row[60:]
            output_lines += new_row
        with open(output_path,'w') as file:
            file.write(output_lines)
    return bad_case

        
