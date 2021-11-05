# coding=utf-8
import os
import sys
import fileinput
import argparse

# parser = argparse.ArgumentParser(
#     description='''This script will append a string to all headers in a fasta file
#     taking a file of strings(desired suffix) to a directory of fasta files ''')
# parser.add_argument('--inputFasta', type=str, help='directory path to fasta files')
# parser.add_argument('--headerFile', type=str, help='File listing the strings to append - on to each fasta file')
# args=parser.parse_args()


def append_header_text(dir_path_fasta, path_headers):
    f = open(path_headers, "r")
    headers = f.read()
    headers = headers.split("\n")

    if (len(headers[-1]) < 2):
        headers = headers[:-1]

    files = []
    for filename in os.listdir(dir_path_fasta):
        if filename.endswith(".fna") or filename.endswith(".fasta"):
            files.append(filename)

    combinedList = [files, headers]
    check_lists_are_same_length(combinedList)

    files.sort(reverse=True)

    files_headers = dict(zip(files, headers))

    for keys, values in files_headers.items():
        print(keys)
        input_file = fileinput.input(dir_path_fasta + keys, inplace=1)
        for line in input_file:
            if line.startswith(">"):
                print(line.strip() + str(' ') +str(values))
            else:
                print (line),
        input_file.close()


def check_lists_are_same_length(lists):
    it = iter(lists)
    the_len = len(next(it))
    if not all(len(l) == the_len for l in it):
        raise ValueError('not all lists have same length!')

if __name__ == '__main__':

    # Path_to_fasta_files = sys.argv[1]
    # Path_to_header_index = sys.argv[2]

    # append_header_text(dir_path_fasta=str(Path_to_fasta_files), path_headers=str(Path_to_header_index))
    append_header_text(dir_path_fasta = "/Users/Maja/Documents/DTU/thesis/Andet/sample_test_set/", path_headers ="/Users/Maja/Documents/DTU/thesis/Andet/sample_test_set/ref_headers.txt")