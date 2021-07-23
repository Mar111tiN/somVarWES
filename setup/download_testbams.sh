#!/bin/sh


DATATXT="https://www.dropbox.com/s/36ye4pdhva7hrrd/bamlinks.txt"


# store current folder
CURRENT=$(pwd);
if [ -z "$1" ]; then
    echo "Please provide a folder as argument!"
else
    FOLDER="$1";
    echo "Downloading data to $FOLDER";
    cd $FOLDER;
    mkdir -p bam && cd bam;
    # download the links txt
    wget $DATATXT;
    # download the files from the link list
    cat bamlinks.txt | xargs wget;
    # check the files and remove md5 files for OK files
    for MD in *.md5; do md5sum -c $MD && rm $MD; done && echo "Testdata downloaded to $FOLDER" && rm bamlinks.txt;
    cd $CURRENT;
fi