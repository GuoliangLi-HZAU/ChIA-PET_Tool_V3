#!/bin/bash

MAX_LINES_PER_CHUNK=10000000
ORIGINAL_FILE=$1
SORTED_FILE=$2
CHUNK_FILE_PREFIX=$ORIGINAL_FILE.split.
SORTED_CHUNK_FILES=$CHUNK_FILE_PREFIX*.sorted

        usage ()
        {
                echo Parallel sort
                echo usage: psort file1 file2
                echo Sorts text file file1 and stores the output in file2
                echo Note: file1 will be split in chunks up to $MAX_LINES_PER_CHUNK lines
                echo  and each chunk will be sorted in parallel
                exit;
        }

        # test if we have two arguments on the command line
        if [ $# != 2 ]
        then
                usage
        fi

        #Cleanup any lefover files
        rm -f $SORTED_CHUNK_FILES > /dev/null
        rm -f $CHUNK_FILE_PREFIX* > /dev/null
        rm -f $SORTED_FILE

        #Splitting $ORIGINAL_FILE into chunks ...
        split -l $MAX_LINES_PER_CHUNK $ORIGINAL_FILE $CHUNK_FILE_PREFIX

        for file in $CHUNK_FILE_PREFIX*
        do
                LC_ALL=C sort -k1,1V -k4,4V -k2,2n $file > $file.sorted &
        done
        wait
        
        #Merging chunks to $SORTED_FILE ...
        LC_ALL=C sort -m -k1,1V -k2,2n $SORTED_CHUNK_FILES > $SORTED_FILE

        rm -f $SORTED_CHUNK_FILES > /dev/null
        rm -f $CHUNK_FILE_PREFIX* > /dev/null