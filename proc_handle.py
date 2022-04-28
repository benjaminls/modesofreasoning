from enum import unique
import os
import sys

if __name__ == "__main__":
    path = os.getcwd()
    print(f'Using {path} as current working directory.')

    # check if "archive" dir doesn't exists, if so, make it
    if not os.path.isdir(os.path.join(path, 'archive')):
        os.mkdir('archive')
        print(
            'Created new directory "archive" to store old data files.\n',
            'If you accidentally overwrite something, try checking there!'
        )
    archive = os.path.join(path, 'archive')

    # check if "proc" dir doesn't exist, if so, make it
    if not os.path.isdir(os.path.join(path, 'proc')):
        os.mkdir('proc')
        print(
            'Created "proc" directory.\n'
            'This will be used to store temporary data files.'
        )
    
    # check if proc is not empty, if so, then move to archive
    proc = os.path.join(path, 'proc')
    if len(os.listdir(proc)) != 0:
        # find untaken name in archive
        basename = 'old_'
        unique_num = 0
        try_name = basename + str(unique_num)
        if try_name in os.listdir(archive):
            unique_num += 1
            try_name = basename + str(unique_num)
        
        os.rename(proc, os.path.join(archive, try_name))
        os.mkdir(proc)
        
