from Bio import SeqIO
import os
from datetime import datetime
import pandas as pd

def running_message(function):
    def wrapper(*args, **kwargs):
        def format_argument(arg):
            if isinstance(arg, pd.DataFrame):
                return f"DataFrame({len(arg)} rows x {len(arg.columns)} columns)"
            elif isinstance(arg, (list, dict)) and len(arg) > 10:
                return f"{type(arg).__name__}({len(arg)} items)"
            return repr(arg)

        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        arg_names = function.__code__.co_varnames[:function.__code__.co_argcount]
        args_repr = [f"{arg}={format_argument(a)}" for arg, a in zip(arg_names, args)]
        kwargs_repr = [f"{k}={format_argument(v)}" for k, v in kwargs.items()]
        signature = ", ".join(args_repr + kwargs_repr)
        print("----------------------------------------")
        print(f"Time: {current_time} - Running {function.__name__}\n\nUsing inputs:\n{function.__name__}({signature})\n")
        result = function(*args, **kwargs)

        now2 = datetime.now()
        current_time2 = now2.strftime("%H:%M:%S")
        print(f"\nTime: {current_time2} - {function.__name__} Completed\n\n")
        return result
    return wrapper

def read_fasta(fastafile: str)->list:
    '''
    Reads fastafile
    :param fastafile: file path
    :return: list of fasta record
    '''
    recordlist = []
    # open file
    with open(fastafile) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            recordlist.append(record)
    return recordlist


def write_fasta(outpath: str, recordlist: list)->None:
    '''
    Writes a fasta file to a given location
    :param path: location to write fasta
    :param recordlist: list of fasta records
    :return: None
    '''
    with open(outpath, "w") as output_handle:
        SeqIO.write(recordlist, output_handle, "fasta")


def makedir(path: str, force_make: bool=False)->str:
    '''
    Makes a directory if the given direcotry doesn't exist. If force_make true, makes a directory with a number
    :param path: location of new directory
    :param force_make: To make a new directory with a number if a directory already exists at given path
    :return: path to new dir
    '''
    try:
        os.mkdir(path)
        return path
    except:
        if not force_make:
            return path
    i = 1
    if force_make:
        while True:
            new_name = f"{path}_{i}"
            try:
                os.mkdir(new_name)
                break
            except:
                i += 1
        return new_name


