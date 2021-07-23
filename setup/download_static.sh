#!/bin/sh

STATICTXT="https://www.dropbox.com/s/urmw86jw3xv7mtb/static_links.txt"

# store current folder
CURRENT=$(pwd);
if [ -z "$1" ]; then
    echo "Please provide a static folder as argument!"
else
    FOLDER="$1";
    echo "Downloading static data to $FOLDER";
    cd $FOLDER;
    # download the links txt
    wget $STATICTXT;
    # download the files from the link list
    cat static_links.txt | xargs wget;

    # check and expand static folder
    STATIC="hg38_static.tar.gz"
    echo "Expanding static tar.gz files and checking integrity";
    md5sum -c ${STATIC}.md5 && tar -xzvf $STATIC && rm ${STATIC}* && echo "Static data downloaded to $FOLDER";
    # make a PON folder and move the PON.tars into that folder
    mkdir -p hg38_static/PON && mv PON.tar* hg38_static/PON && cd hg38_static/PON;

    # check and expand the TESTPON into tatic folder
    echo "Expanding PON.tar files and checking integrity";
    md5sum -c PON.tar.md5 && tar -xvf PON.tar && rm PON.tar* && echo "PON data downloaded to ${FOLDER}/hg38_static/PON/HAEv7_hg38_NovaSeq";
    cd $FOLDER && rm static_links.txt;
    cd $CURRENT;
fi