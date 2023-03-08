import pandas as pd
import os
import time
import concurrent.futures # for IO intensive work
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
import json
import urllib.request
import numpy as np
from multiprocessing import Lock, Process, Manager #for computation intensive work
import subprocess
import shutil

def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record

def get_assemblies(gca_id, download=False, path = "Data/", email = "a@email.address"):
    """Download genbank assembly metadata and ftp address according to a GCA number.
    Args:
        gca_id: a GCA numbers
        download: whether to download the assembly. Don't use this function to download a large amount of genome. 
        path: folder to save data
    """
    try:
        Entrez.email = email

        handle = Entrez.esearch(db="assembly", term= gca_id, retmax=1)
        record = Entrez.read(handle)
        ids = record['IdList']

        #get summary
        summary = get_assembly_summary(ids)
        
        #three path to save output
        path_jsons = path + "jsons/"
        path_ftp = path + "ftp/"
        path_genome = path + "download_genome/"

        for p in [path_jsons, path_ftp, path_genome]:
            if os.path.exists(p):
                pass
            else:
                os.mkdir(p)

        # save raw summary to json
        with open(path_jsons + gca_id + ".json", "w") as outfile:
            json.dump(summary, outfile)
        
        #get ftp link
        url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
        label = os.path.basename(url)
        #get the fasta link - change this to get other formats
        link = os.path.join(url,label+'_genomic.fna.gz')

        #save ftp link
        op = open(path_ftp + gca_id + ".txt", "w")
        op.write(link)
        op.close()

        if download == True:
            #download link
            urllib.request.urlretrieve(link, path_genome + f'{label}_genomic.fna.gz')
    
    #error handling
    except Exception as ex:
        error_report = gca_id + str((type(ex))) + " " + str(ex.args) + "\n"
        op = open("link_download_error.txt", "a")
        op.write(error_report)
        op.close()

def ftp_download(gca_list, worker = 2, download_path = "Data/", email = "a@email.address"):
    """Download genbank assembly metadata and ftp address according to a list of GCA number.
    Args: 
        worker: threads number
    """
    if os.path.exists(download_path):
        pass
    else:
        os.mkdir(download_path)

    op = open("link_download_error.txt", "w")
    op.write("")
    op.close()

    start_time_external = time.time()
    print(time.asctime())

    if __name__ == 'blast_at_local_tools':
        
        exe=concurrent.futures.ThreadPoolExecutor(max_workers=worker) 


        for gca in gca_list:
            exe.submit(get_assemblies, gca, False, download_path, email)
            #exe.shutdown()

        exe.shutdown()
    
    end_time_external = time.time()
    print("final time usage " + str(end_time_external - start_time_external))
    print(time.asctime())

def ftp_re_download(error_file = "link_download_error.txt", worker = 2, download_path = "Data/", email = "a@email.address"):
    """Download again genbank assembly metadata and ftp address according to a list of GCA number.
    Args: 
        worker: threads number
    """
    
    #check these assembly IDs of which have not been successfully downloaded
    error_1 = pd.read_csv(error_file, sep="<",names = ["gca", "rest"])
    print(str(len(error_1)) + " genome metadata still need to be downloaded.")

    op = open("link_download_error.txt", "w")
    op.write("")
    op.close()

    start_time_external = time.time()
    print(time.asctime())

    if __name__ == 'blast_at_local_tools':
        exe=concurrent.futures.ThreadPoolExecutor(max_workers=worker) 


        for gca in error_1.gca:
            exe.submit(get_assemblies, gca, False, download_path, email)
            #exe.shutdown()

        exe.shutdown()

    end_time_external = time.time()
    print("final time usage " + str(end_time_external - start_time_external))
    print(time.asctime())

def ftp_modify(readins, batch, ftp_path = "Data/ftp/", rsync_path = "Data/rsync/"):
    """function for changing ftp address to rsync address
    """
    
    if os.path.exists(rsync_path):
        pass
    else:
        os.mkdir(rsync_path)
    
    
    for x in readins:
        with open (ftp_path + x, "r") as f1:
            address = f1.read().replace('ftp', 'rsync', 1).replace('\\', '/').replace("_genomic.fna.gz","_genomic.fna.gz\n")
            f1.close()    
        with open (rsync_path + "address" + str(batch) +".txt", "a") as f2:# write according to batch No, batch No. should larger than threads number
            f2.write(address)
            f2.close()

def ftp_to_rsync(ftp_path = "Data/ftp/", rsync_path = "Data/rsync/", process_num = 1):
    """change ftp address to rsync address in muilti threads"""
    
    #take file names
    file_names = os.listdir(ftp_path)
    #splite file names based on process number
    mission_lst = list(np.array_split(file_names, process_num))

    start_time_external = time.time()
    print(time.asctime())

    if __name__ == "blast_at_local_tools":

        jobs = []

        for x in range(0, process_num):

            p = Process(target = ftp_modify,
                        args = (list(mission_lst[x]), 
                                x,
                                ftp_path,
                                rsync_path
                               )
                       )
            p.start()
            jobs.append(p)


        for z in jobs:
            z.join()

    end_time_external = time.time()
    print("final time usage " + str(end_time_external - start_time_external))
    print(time.asctime())

def rsync_one(address, save_dir):
    """Download one genome data from address and save to save_dir
    Args:
        address: the rsync address for one genome
        save_dir: the local dir to save the genome
    """
    
    try:
        subprocess.run(["rsync", "-avzq", address, save_dir]) #use python to run the code in the list in terminal
    except Exception as ex:
        error_report = address + " <> " + str((type(ex))) + " " + str(ex.args) + "\n"
        op = open("rsync_download_error.txt", "a")
        op.write(error_report)
        op.close()

def genome_download(genome_path = "Data/download_genome/", rsync_path = "Data/rsync/", worker = 45):
    """use multithreads to download genome data"""

    if os.path.exists(genome_path):
        pass
    else:
        os.mkdir(genome_path)

    #merge all the addresses into one mission list
    ls_mission = [] # the list of all the address 
    for x in os.listdir(rsync_path):
        f = open(rsync_path + x)
        addresses = f.readlines()
        f.close()
        addresses = [address.rstrip('\n') for address in addresses]

        ls_mission.extend(addresses)    
    
    start_time_external = time.time()
    print(time.asctime())
    
    if __name__ == 'blast_at_local_tools':
        exe=concurrent.futures.ThreadPoolExecutor(max_workers=worker)         #45 threads


        for address in ls_mission:
            exe.submit(rsync_one, address, genome_path)

        exe.shutdown(wait=True)
    
    end_time_external = time.time()
    print("final time usage " + str(end_time_external - start_time_external))
    print(time.asctime())

def genome_re_download(genome_path = "Data/download_genome/", rsync_path = "Data/rsync/", worker = 45):
    """Download again genbank assembly metadata and ftp address according to a list of GCA number.
    Args: 
        worker: threads number
    """
    #merge all the addresses into one mission list
    ls_mission = [] # the list of all the address 
    for x in os.listdir(rsync_path):
        f = open(rsync_path + x)
        addresses = f.readlines()
        f.close()
        addresses = [address.rstrip('\n') for address in addresses]

        ls_mission.extend(addresses)
    
    
    #check the genome that haven't been successfully downloaded
    ls_download = os.listdir(genome_path)#take the name that already been downloaded

    dict_name_address = {} #make a dict to match genome file name and genome rsync address
    for address in ls_mission:
        name = address.split('/')[-1]
        dict_name_address[name] = address

    ls_not_download = set(dict_name_address.keys()).difference(set(ls_download)) #genome have not been download
    print(str(len(ls_not_download)) + " genome data still need to be downloaded.")

    op = open("rsync_download_error.txt", "w")
    op.write("")
    op.close()

    start_time_external = time.time()
    print(time.asctime())

    #use multithreads to download genome data that are not successful previously
    if __name__ == 'blast_at_local_tools':
        exe=concurrent.futures.ThreadPoolExecutor(max_workers=worker)         #45 threads


        for name in ls_not_download:
            address = dict_name_address[name]
            exe.submit(rsync_one, address, genome_path)

        exe.shutdown(wait=True)
    
    end_time_external = time.time()
    print("final time usage " + str(end_time_external - start_time_external))
    print(time.asctime())

def md5_address(readins, batch, ftp_path = "Data/ftp/", md5_address_path = "Data/md5_address/"):
    """function to change ftp address to md5 address"""
    if os.path.exists(md5_address_path):
        pass
    else:
        os.mkdir(md5_address_path)

    for x in readins:
        with open (ftp_path + x, "r") as f1:
            address = f1.read().replace('ftp', 'rsync', 1).replace('\\', '/').split("/")
            address = address[0:-1]
            address.append('md5checksums.txt\n')
            address = "/".join(address)
            f1.close()    
        with open (md5_address_path + str(batch) +".txt", "a") as f2:# write according to batch No, batch No. should larger than threads number
            f2.write(address)
            f2.close()

def ftp_to_md5(ftp_path = "Data/ftp/", md5_address_path = "Data/md5_address/", process_num = 1):
    """change ftp address to md5 address in muilti threads"""


        
    #take file names
    file_names = os.listdir(ftp_path)
    #splite file names based on process number
    mission_lst = list(np.array_split(file_names, process_num))
    
    
    #change ftp address to md5 address in muilti threads
    start_time_external = time.time()
    print(time.asctime())

    if __name__ == "blast_at_local_tools":

        jobs = []

        for x in range(0, process_num):

            p = Process(target = md5_address,
                        args = (list(mission_lst[x]), 
                                x,
                                ftp_path,
                                md5_address_path
                               )
                       )
            p.start()
            jobs.append(p)


        for z in jobs:
            z.join()

    end_time_external = time.time()
    print("final time usage " + str(end_time_external - start_time_external))
    print(time.asctime())

def rsync_one_md5(address, save_dir):
    """Download one MD5 file from address and save to save_dir
    Args:
        address: the rsync address for one MD5 file
        save_dir: the local dir to save the MD5 file
    """
    
    try:
        name = address.split("/")[-2] + "_md5.txt"
        subprocess.run(["rsync", "-avzq", address, save_dir + name]) #use python to run the code in the list in terminal
    except Exception as ex:
        error_report = address + str((type(ex))) + " " + str(ex.args) + "\n"
        op = open("md5_download_error.txt", "a")
        op.write(error_report)
        op.close()

def md5_download(md5_download_path = "Data/download_md5/", md5_address_path = "Data/md5_address/", worker = 45):
    """use multithreads to download md5 file"""

    if os.path.exists(md5_download_path):
        pass
    else:
        os.mkdir(md5_download_path)

    op = open("md5_download_error.txt", "w")
    op.write("")
    op.close()
    
    #merge all the md5 addresses into one mission list
    ls_mission_md5 = [] # the list of all the address 
    for x in os.listdir(md5_address_path):
        f = open(md5_address_path + x)
        addresses = f.readlines()
        f.close()
        addresses = [address.rstrip('\n') for address in addresses]

        ls_mission_md5.extend(addresses)

    start_time_external = time.time()
    print(time.asctime())

    #use multithreads to download md5 address
    if __name__ == 'blast_at_local_tools':
        exe=concurrent.futures.ThreadPoolExecutor(max_workers=worker)         #45 threads


        for address in ls_mission_md5:
            exe.submit(rsync_one_md5, address, md5_download_path)

        exe.shutdown(wait=True)

    end_time_external = time.time()
    print("final time usage " + str(end_time_external - start_time_external))
    print(time.asctime())

def md5_re_download(md5_download_path = "Data/download_md5/", md5_address_path = "Data/md5_address/", worker = 45):
    """Download again md5 files 
    Args: 
        worker: threads number
    """
    
    #check the md5 that haven't been successfully downloaded
    ls_download_md5 = os.listdir(md5_download_path)#take the name that already been downloaded
    ls_download_md5 = [x.replace("_md5.txt", "") for x in ls_download_md5 ]

    #merge all the md5 addresses into one mission list
    ls_mission_md5 = [] # the list of all the address 
    for x in os.listdir(md5_address_path):
        f = open(md5_address_path + x)
        addresses = f.readlines()
        f.close()
        addresses = [address.rstrip('\n') for address in addresses]

        ls_mission_md5.extend(addresses)

    dict_name_address = {} #make a dict to match md5 address
    for address in ls_mission_md5:
        name = address.split('/')[-2]
        dict_name_address[name] = address

    ls_not_download_md5 = set(dict_name_address.keys()).difference(set(ls_download_md5)) #md5 have not been download
    print(str(len(ls_not_download_md5)) + " md5 files still need to be downloaded.")
    
    start_time_external = time.time()
    print(time.asctime())

    #use multithreads to download md5 files that are not successful previously
    if __name__ == 'blast_at_local_tools':
        exe=concurrent.futures.ThreadPoolExecutor(max_workers=worker)         #45 threads


        for name in ls_not_download_md5:
            address = dict_name_address[name]
            exe.submit(rsync_one_md5, address, md5_download_path)

        exe.shutdown(wait=True)

    end_time_external = time.time()
    print("final time usage " + str(end_time_external - start_time_external))
    print(time.asctime())

def md5_sum(addresses, batch, genome_path = "Data/download_genome/", md5_generate_path = "Data/generated_md5/"):
    """function to generate md5 file"""
    for address in addresses:
        try:
            address = genome_path + address
            md5_value = subprocess.run(["md5sum", address], capture_output=True, text=True)
            #save to file
            with open(md5_generate_path + "md5_" + str(batch) +  ".txt", "a") as f:
                        f.write(str(md5_value.stdout))
        except Exception as ex:
            error_report = address + " <> " + str((type(ex))) + " " + str(ex.args) + "\n"
            op = open("md5_generated_error.txt", "a")
            op.write(error_report)
            op.close()

def md5_generate(md5_generate_path = "Data/generate_md5/", genome_path = "Data/download_genome/", process_num = 1):
    """use multithreads to generate md5 file"""
    
    if os.path.exists(md5_generate_path):
        pass
    else:
        os.mkdir(md5_generate_path)
       
    #take file names
    genome_names = os.listdir(genome_path)
    #splite genome files based on process number
    mission_lst = list(np.array_split(genome_names, process_num))

    #generate md5 files in muilti threads
    start_time_external = time.time()
    print(time.asctime())

    if __name__ == "blast_at_local_tools":

        jobs = []

        for x in range(0, process_num):

            p = Process(target = md5_sum,
                        args = (list(mission_lst[x]), 
                                x,
                                genome_path,
                                md5_generate_path
                               )
                       )
            p.start()
            jobs.append(p)


        for z in jobs:
            z.join()

    end_time_external = time.time()
    print("final time usage " + str(end_time_external - start_time_external))
    print(time.asctime())
        
def md5sum_check(the_tup_lst, shared_lst, md5_download_path = "Data/download_md5/"):
    """a function to check the downloaded md5 is equal to locally generated md5
    Arg:
        the_tup_lst: a list of tuples
            tuple[0] is the path of download genome
            tuple[1] is the locally generated md5
        shared_lst: a list shared by multiprocess
    """
    for the_tup in the_tup_lst:
        genome_ID = the_tup[0].split('/')[-1].replace('_genomic.fna.gz', '')

        download_md5 = pd.read_csv(md5_download_path + genome_ID + '_md5.txt',
                              sep = " ",
                              header = None).set_index([2])

        if download_md5.loc['./' + genome_ID + '_genomic.fna.gz', 0] == the_tup[1]:
            pass
        else:
            shared_lst.append(genome_ID)

def md5_check(md5_generate_path = "Data/generate_md5/", 
              md5_download_path = "Data/download_md5/",
              process_num = 1, show_not_match=False):
    
    #merge all the addresses into one mission list
    generated_md5 = pd.DataFrame()  #all the md5 value
    for x in os.listdir(md5_generate_path):
        f = pd.read_csv(md5_generate_path + x,
                        sep = " ",
                        header = None).set_index([2])



        generated_md5 = pd.concat([generated_md5, f])
        
    #splite generated md5 tuple based on process number
    mission_lst = list(generated_md5.itertuples(index=True, name=None))    
    mission_lst = list(np.array_split(mission_lst, process_num))
    
    
    #check md5
    if __name__ == "blast_at_local_tools":
        with Manager() as manager:
            lst = manager.list() #the shared list

            jobs = []

            for x in range(0, process_num):

                p = Process(target = md5sum_check,
                            args = (mission_lst[x], lst, md5_download_path)
                           )
                p.start()
                jobs.append(p)


            for z in jobs:
                z.join()

            lst_ = [] # the list to show

            for x in lst: #get item from the shared list and put to the normal list
                lst_.append(x)

        # check how many md5 cannot match 
        print(str(len(lst_)) + " MD5 values cannot match.")
        if show_not_match:
            print(lst_)

def g_unzip(genome_path = "Data/download_genome/"):
    """unzip the downloaded genome data"""
    gunzip_command = "gunzip " + genome_path + "*.gz"
    subprocess.run(gunzip_command, shell = True)

def make_a_db(in_abs_pth, 
              db_dir_abs_pth, 
              dbtype,
              makeblastdb_bin):
    """make a BLAST database based on one genome data
    Args:
        in_abs_pth: path of input genome data
        db_dir_abs_pth: path to save database
        dbtype: Molecule type of target db
        makeblastdb_bin: location of makeblastdb
    """
    strain_naam = os.path.basename(in_abs_pth)
    strain_naam = os.path.splitext(strain_naam)[0]
    
    save_dest = os.path.join(db_dir_abs_pth, strain_naam)
    
    if os.path.exists(save_dest):
        pass
    else:
        os.mkdir(save_dest)
        
    save_naam = os.path.join(save_dest, strain_naam)
        
    subprocess.run([makeblastdb_bin, 
                    "-dbtype", dbtype, 
                    "-in", in_abs_pth,
                    "-out", save_naam,
                    "-parse_seqids",
                    "-logfile",
                    "/dev/null"
                   ])
    
def make_db_by_ls(in_abs_pth_lst,
                 db_dir_abs_pth_,
                 dbtype_,
                 makeblastdb_bin_):
    """make a list of BLAST database based on a list of genome data
    Args:
        in_abs_pth_ls: path of input genome data list
        db_dir_abs_pth_: path to save each database
        dbtype_: Molecule type of target db
        makeblastdb_bin_: location of makeblastdb
    """
    for s in in_abs_pth_lst:
        make_a_db(s, 
                  db_dir_abs_pth_, 
                  dbtype_,
                  makeblastdb_bin_)
        
def make_database(genome_path = "Data/download_genome/",
                  blastdb_path = "Data/blast_db/",
                  dbtype = "nucl",
                  makeblastdb_bin = "makeblastdb", 
                  process_num = 1):
    """use multithreads to make BLAST database"""

    start_time_external = time.time()
    print(time.asctime())

    if __name__ == "blast_at_local_tools":

        # mkdir or not
        if os.path.exists(blastdb_path):
            pass
        else:
            os.mkdir(blastdb_path)

        mission_lst = os.listdir(genome_path)

        mission_lst = [os.path.join(genome_path, x) for x in mission_lst]

        mission_lst = list(np.array_split(mission_lst, 
                                          process_num
                                         )
                          )

        # mp
        jobs = []

        for sub_mission_lst in mission_lst:

            p = Process(target = make_db_by_ls,
                        args = (sub_mission_lst, 
                                blastdb_path, 
                                dbtype,
                                makeblastdb_bin
                               )
                       )
            p.start()
            jobs.append(p)


        for z in jobs:
            z.join()

    end_time_external = time.time()
    print("final time usage " + str(end_time_external - start_time_external))
    print(time.asctime())

def sequential_blast_high(catted_fasta_pth,
                          db_pth, # alread abs
                          blast_res_dir, # already abs
                          blast_bin,
                          e_value,
                          t_num):
    """Run blast for query against one db
        catted_fasta_pth: path to query file 
        db_pth: path to db
        blast_res_dir: path to blast result
        blast_bin: path to blastn
        e_value: Expectation value (E) threshold for saving hits
        t_num: Number of threads (CPUs) to use in the BLAST search
    """
#   get db name to blast
    strain_naam = db_pth.split('/')[-1]
    
    db_name = os.path.join(db_pth, strain_naam)
    save_dest = os.path.join(blast_res_dir, strain_naam + '.tab')
        
    subprocess.run([blast_bin, 
                    "-query", catted_fasta_pth,
                    "-db", db_name,
                    "-out", save_dest,
                    "-evalue", str(e_value),
                    "-max_hsps", "1",
                    "-dust", "no",
                    "-outfmt", str("6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore slen"),
                    "-mt_mode", "1",
                    "-num_threads", str(t_num)
                   ])
    
def sequential_blast_high_s(catted_fasta_pth_,
                            db_pth_lst, # alread abs
                            blast_res_dir_, # already abs
                            blast_bin_,
                            e_value_,
                            pro_id_, 
                            s_t_d,
                            st_lock_,
                            available_threads):
    """Run blast for query against a list of db"""
    for db in db_pth_lst:
        st_lock_.acquire()
        if sum(s_t_d.values()) >= available_threads:
            useable_th = s_t_d[pro_id_]
        else:
            s_t_d[pro_id_] = s_t_d[pro_id_] + (available_threads - sum(s_t_d.values()
                                                                      )
                                              )
            useable_th = s_t_d[pro_id_]

        st_lock_.release()
                                               
        sequential_blast_high(catted_fasta_pth_,
                              db,
                              blast_res_dir_, 
                              blast_bin_, 
                              e_value_,
                              useable_th)
        
    st_lock_.acquire()
    s_t_d[pro_id_] = 0
    
    st_lock_.release()

def blast(blastdb_path = "Data/blast_db/",
          query_path = 'Data/example_query.fas', 
          blast_output_path = 'Data/blast_output/',
          blast_bin = 'blastn', 
          E_value = 1e-5, 
          process_num = 1):
    """"use multithreads to run BLAST"""


    start_time_external = time.time()
    print(time.asctime())

    # mkdir or not
    if os.path.exists(blast_output_path):
        pass
    else:
        os.mkdir(blast_output_path)
    
    
    mission_lst = os.listdir(blastdb_path)
    mission_lst = [os.path.join(blastdb_path, x) for x in mission_lst]

    mission_lst = list(np.array_split(mission_lst, 
                                      process_num
                                     )
                      )

    if __name__ == "blast_at_local_tools":
        with Manager() as manager:
            smart_threads_d = manager.dict()

            st_lock = Lock()

            for x in range(process_num):
                smart_threads_d[x] = 1

            # mp
            jobs = []

            pro_id = 0
            for db_lst in mission_lst:

                p = Process(target = sequential_blast_high_s,
                            args = (query_path, 
                                    db_lst, 
                                    blast_output_path,
                                    blast_bin,
                                    E_value,
                                    pro_id,
                                    smart_threads_d, 
                                    st_lock,
                                    process_num
                                   )
                           )
                p.start()
                jobs.append(p)

                pro_id += 1


            for z in jobs:
                z.join()

    end_time_external = time.time()
    print("final time usage " + str(end_time_external - start_time_external))
    print(time.asctime())

def database_remove_old(GCA_list_remove,
                       ftp_path = "Data/ftp/",
                       blastdb_path = "Data/blast_db/",
                       achive_removed_blastdb_path =  "Data/blast_db_removed/"):
    """Remove a list of database based on a list of gca id"""

    if os.path.exists(achive_removed_blastdb_path):
        pass
    else:
        os.mkdir(achive_removed_blastdb_path)

    for i in GCA_list_remove:
        with open (ftp_path + i + ".txt", "r") as f:
            db_name = f.read().replace(".fna.gz", "").split("/")[-1]
            db_from_path = blastdb_path + db_name
            db_to_path = achive_removed_blastdb_path + db_name

            shutil.move(db_from_path, db_to_path)

def extract_seq(df, process_id, 
                blast_op_path, 
                blast_db_path, 
                extract_seq_path, 
                dtypes,
                blastdbcmd_bin ):
    """function to extract sequence based on one blast output df
    Note: the hit sequences that reach the start or end of the subject sequences are filtered out
    """ 

    df_path = blast_op_path + "/" + df
    db_path = blast_db_path + "/" + df[0:-4] + "/" + df[0:-4]
    output_path = extract_seq_path
    
    try: 
        blast_df = pd.read_csv(df_path, sep="\t", header=None, dtype=dtypes)

        #for loop of each query
        for i in blast_df.index:
            s_start = blast_df.loc[i, 8]
            s_end = blast_df.loc[i, 9]
            s_tot_len = blast_df.loc[i, 12]
            s_ls = [s_start, s_end]
            if (1 not in s_ls) & (s_tot_len not in s_ls): #check the s_start or s_end equal to 1 or s_total_length 
                #take the sequence from db, use blastcmd

                #check forward or reverse
                strand = "plus"
                out_range = str(s_start) + "-" + str(s_end)
                if s_start > s_end: #if reverse
                    strand = "minus"
                    out_range = str(s_end) + "-" + str(s_start) 

                out_seq_record = subprocess.run([
                                blastdbcmd_bin,
                                "-db", db_path,
                                "-entry", blast_df.loc[i, 1],
                                "-line_length", "60",
                                "-range", out_range,
                                "-strand", strand
                            ], capture_output=True, text=True)

                #save to file
                with open(output_path + str(process_id) + "_" + blast_df.loc[i, 0] + ".fas", "a") as f:
                    f.write(str(out_seq_record.stdout))
    except Exception as ex:
        error_report = df_path + " <> " + str((type(ex))) + " " + str(ex.args) + "\n"
        op = open("extract_seq_error.txt", "a")
        op.write(error_report)
        op.close()
                
def extract_seq_list(df_list, process_id, blast_op_path, blast_db_path, extract_seq_path, dtypes, blastdbcmd_bin):
    """function to extract sequence based on a list of blast output df""" 
    for x in df_list:
        extract_seq(x, process_id, blast_op_path, blast_db_path, extract_seq_path, dtypes, blastdbcmd_bin)

def merge_seq(query, extract_seq_path):
    """function to merge extract sequence output from different threads"""
    
    seq_file_ls = extract_seq_path  + "*_" + query + ".fas"
    total_extract_seq = extract_seq_path  + query + "_total.fas"
    cat_command = "cat " + seq_file_ls + " > " + total_extract_seq 
    rm_command = "rm " + seq_file_ls
    
    subprocess.run(cat_command, shell = True)
    subprocess.run(rm_command, shell = True)
        
def blast_result_seq(blastdb_path = "Data/blast_db/",
                     blast_output_path = 'Data/blast_output/',
                     query_path = 'Data/example_query.fas',
                     extract_seq_path = "Data/extract_seq/",
                     dtypes = {0: str, 1: str, 2: float, 3: int, 4: int, 5: int, 6: int, 7: int, 8: int, 9: int, 10: float, 11: float, 12:int},
                     blastdbcmd_bin = "blastdbcmd",
                     process_num = 1):
    """Extract the sequences from the BLAST results in several process
    Note: the hit sequences that reach the start or end of the subject sequences are filtered out
    Args: 
        blastdb_path:
        blast_output_path: path to blast output tab file
        query_path: path to query file for blast
        extract_seq_path: path to save the extract sequences
        dtypes: data type for each column in blast output tab file
        blastdbcmd_bin: blastdbcmd location
        process_num: threads number
    """
    
    if os.path.exists(extract_seq_path):
        pass
    else:
        os.mkdir(extract_seq_path)
    
    op = open("extract_seq_error.txt", "w")
    op.write("")
    op.close()
    
    #take output file names
    file_names = os.listdir(blast_output_path)
    
    #splite blast result based on process number
    mission_lst = list(np.array_split(file_names, process_num))

    #get query list
    dict_1 = SeqIO.to_dict(SeqIO.parse(query_path,
                               "fasta"))
    query_list = dict_1.keys()
    
    #extract sequence in muilti threads
    start_time_external = time.time()
    print(time.asctime())

    if __name__ == "blast_at_local_tools":

        jobs = []

        for x in range(0, process_num):

            p = Process(target = extract_seq_list,
                        args = (list(mission_lst[x]), 
                                x, 
                                blast_output_path, 
                                blastdb_path, 
                                extract_seq_path,
                                dtypes,
                                blastdbcmd_bin
                               )
                       )
            p.start()
            jobs.append(p)


        for z in jobs:
            z.join()

    for q in list(query_list):
        merge_seq(q, extract_seq_path)

    end_time_external = time.time()
    print("final time usage " + str(end_time_external - start_time_external))
    print(time.asctime())

def extract_tab(df,
                process_id,
                blast_op_path,
                extract_tab_path,
                dtypes): 
    """extract BLAST result from one output file"""
    
    df_path = blast_op_path + "/" + df
    output_path = extract_tab_path

    try:
        blast_df = pd.read_csv(df_path, sep="\t", header=None, dtype=dtypes)
        
        for i in set(blast_df.index):

            q = blast_df.loc[i, 0]

            save_df_dir = output_path + str(process_id) + "_" + q + ".tab"
            blast_df.loc[[i]].to_csv(save_df_dir, header=False, sep="\t", mode = "a")
    except Exception as ex:
        error_report = df_path + " <> " + str((type(ex))) + " " + str(ex.args) + "\n"
        op = open("extract_tab_error.txt", "a")
        op.write(error_report)
        op.close()

def extract_tab_list(df_list, process_id, blast_op_path, extract_tab_path, dtypes):
    """function to extract result tab file based on a list of blast output df""" 
    for x in df_list:
        extract_tab(x, process_id, blast_op_path, extract_tab_path, dtypes)

def merge_tab(query, extract_tab_path, columns):
    """function to merge extract tab output from different threads"""
    
    tab_file_ls = extract_tab_path  + "*_" + query + ".tab"
    total_extract_tab = extract_tab_path  + query + "_total.tab"
    cat_command = "cat " + tab_file_ls + " > " + total_extract_tab
    rm_command = "rm " + tab_file_ls
    
    head_line = "\t".join(columns)
    sed_command = "sed -i '1i\\" + head_line + "' " + total_extract_tab
    
    subprocess.run(cat_command, shell = True)
    subprocess.run(rm_command, shell = True)
    subprocess.run(sed_command, shell = True) # add head to each tab file


def blast_result_df(blast_output_path = 'Data/blast_output/',
                     query_path = 'Data/example_query.fas',
                     extract_tab_path = "Data/extract_tab/",
                     dtypes = {0: str, 1: str, 2: float, 3: int, 4: int, 5: int, 6: int, 7: int, 8: int, 9: int, 10: float, 11: float, 12:int},
                     columns = ["query_id", "subject_id", "pct_identity", "alignment_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score", "s_total_length"],
                     process_num = 1):
    """Extract the BLAST result tab file from the BLAST results in several process
    Args: 
        blast_output_path: path to blast output tab file
        query_path: path to query file for blast
        extract_tab_path: path to save the extract tab files
        dtypes: data type for each column in blast output tab file
        columns: the name of each column in the tab files
        process_num: threads number
    """

    if os.path.exists(extract_tab_path):
        pass
    else:
        os.mkdir(extract_tab_path)
        
    op = open("extract_tab_error.txt", "w")
    op.write("")
    op.close()

    #take output file names
    file_names = os.listdir(blast_output_path)

    #splite blast result based on process number
    mission_lst = list(np.array_split(file_names, process_num))

    #get query list
    dict_1 = SeqIO.to_dict(SeqIO.parse(query_path,
                            "fasta"))
    query_list = dict_1.keys()

    #extract sequence in muilti threads
    start_time_external = time.time()
    print(time.asctime())

    if __name__ == "blast_at_local_tools":

        jobs = []

        for x in range(0, process_num):
            p = Process(target = extract_tab_list,
                        args = (list(mission_lst[x]), 
                                x, 
                                blast_output_path, 
                                extract_tab_path,
                                dtypes
                            )
                    )
            p.start()
            jobs.append(p)


        for z in jobs:
            z.join()

    for q in list(query_list):
        merge_tab(q, extract_tab_path, columns)

    end_time_external = time.time()
    print("final time usage " + str(end_time_external - start_time_external))
    print(time.asctime())
