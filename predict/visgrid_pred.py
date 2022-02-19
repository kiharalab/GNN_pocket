import os


def pqr2pdb_vis(ifile,ofile):
    with open(ifile,'r') as file:
        lines = file.readlines()
    new_lines = ''
    spaces = '     '
    for idx,line in enumerate(lines):
        id = spaces+str(idx+1)
        new_lines += line[:22]+id[-4:]+line[26:]
    with open(ofile,'w') as file:
        file.write(new_lines)

def visgrid_pred(bad_case,params):
    dir = params['test_odir']
    vis_pdb_dir = "../dataset/test_processed/vis_pdb"
    output_vis_dir = "../dataset/test_processed/output_vis"
    if not os.path.exists(vis_pdb_dir):
        os.mkdir(vis_pdb_dir)
    if not os.path.exists(output_vis_dir):
        os.mkdir(output_vis_dir)
    tools_dir = './tools'
    vis_pred_exe = os.path.join(tools_dir,'VisGrid/VisGrid')
    os.system("chmod 777 "+ vis_pred_exe)
    for id in bad_case:
        pdb_path = os.path.join(dir,"pdb/"+str(id)+".pdb")
        vis_pdb_path = os.path.join(vis_pdb_dir,str(id)+'.pdb')
        pqr2pdb_vis(pdb_path, vis_pdb_path) # adding residual number to pdb
        output_vis_path = os.path.join(output_vis_dir,str(id)+'.pdb')
        command = (vis_pred_exe+" "+ vis_pdb_path + " >"+ output_vis_path)
        os.system(command) # running visgrid
        with open(output_vis_path,'r') as file:
            lines = file.readlines()
        pocket_1 = lines[2].strip().split(',')
        pocket_2 = lines[4].strip().split(',')
        pocket_3 = lines[6].strip().split(',')
        input_dir = "../dataset/test/"
        output_dir = os.path.join("../final_pred",str(id))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_lines = ''
        pqr_path = os.path.join(input_dir, id+'/structure.pqr')
        output_path = os.path.join(output_dir,'structure.pqr')

        with open(pqr_path,'r') as file:
            rows = file.readlines()
        output_lines = ''
        for idx,row in enumerate(rows):
            if str(idx+1) in pocket_1:
                output_lines += row[:59]+'1'+row[60:]
            elif str(idx+1) in pocket_2:
                output_lines += row[:59]+'2'+row[60:]
            elif str(idx+1) in pocket_3:
                output_lines += row[:59]+'3'+row[60:]
            else:
                output_lines +=row
        with open(output_path,'w') as file:
            file.write(output_lines)




