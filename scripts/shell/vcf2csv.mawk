#!/bin/sh

# only use the header
mawk '
BEGIN {
  ################ VCF2CSV TRANSLATOR ###################
  ########## SPECS ###############
  # here go the standard 8 tab-separated fields as described in the VCF-specs
  vcfspex="Chr,Pos,ID,Ref,Alt,Qual,Filter,Info";
  split(vcfspex,SPEX,",");
  # get the col number for the Spex
  for (col in SPEX) {
    SPEXCOL[SPEX[col]] = col;
  }

  # here come the actual output fields. omitting some specs cols
  # can be flexibly changed
  spexOut="Chr,Pos,Ref,Alt";
  # split the outputSpex into the FIELDS array
  spexCount = split(spexOut,FIELDS,",");

  ######### TAGS #####################
  # get extra tags from command argument $1 (with default for dbSNP)
  # provide in format Tag:Header,Tag[:Header],..)
  # optional Header is used in HeaderColumns
  # different defaults are: 
  common_all="G5A:MinorAlleleAllPop>5,G5:MinorAllele+1";
  dbSNP="FREQ:AlleleFreq,VC:VariantClass";
  varscan="SS:somaticStatus,GPV:variantP,SPV:somaticP,DP:Depth,RD:R1,AD:R2,DP4:readCounts";
  tags="'${1:-SS:somaticStatus,GPV:variantP,SPV:somaticP,DP:Depth,RD:R1,AD:R2,DP4:readCounts}'";
  tagsCount = split(tags, TAGS, ",");

  # split the tags into the TAG array and translations into TAGHEADER
  for (i=0; i++<tagsCount;) {
    # split TagName and Header translation
    split(TAGS[i], TAG, ":");
    TAGS[i]=TAG[1];
    if (TAG[2]) {
        TAGNAME[TAG[1]] = TAG[2];
    } else {
        TAGNAME[TAG[1]] = TAG[1];
    }
    # TAGS [1-> "TAG"]
    # TAGHEADER ["TAG" -> "TAGNAME"]
    tagsFound = 0;
  }
}

readData {
  # SPEX
  # first line extra
  printf("%s", $SPEXCOL[FIELDS[1]]);
  for (i=1; i++<spexCount;) {
    printf("\t%s",$SPEXCOL[FIELDS[i]]);
  }
  ######## INFO FIELDS ##########
  for (i=spexCount; i++< spexCount+IL;) {
    tag=FIELDS[i];
    tagLen=length(tag)+1;
    pattern=tag "=[^\\t;$]+";
    if (match($8, pattern)) {
      printf("\t%s", substr($8,RSTART+tagLen,RLENGTH-tagLen));
    } else {
      printf(".");
    }
  ####### FORMAT FIELDS #########
  }
  if (FL) { # FORMAT IS USED
    if (firstLine) { # find the order of the format tags
      fCount=split($9,FORMAT,":");
      for (i=0;i++<fCount;) {
        FCOL[FORMAT[i]]=i;
      }
      firstLine=0;
    }
    for (s=0; s++<sampleCount;) {
      for (i=spexCount+IL; i++<spexCount+IL+FL;) {
        split($(9+s),FDATA,":");
        printf("\t%s",FDATA[FCOL[FIELDS[i]]]);
      }
    }
  }
  printf("\n");
  next;
}
/^##/ { 
  # get INFO on FIELDS and allocate to info and format arrays

  if (length(tags) == 0) next; # skip if all tags are found
  for (i=0; i++<tagsCount;) {
    tag=TAGS[i];
    tagLen=length(tag);
    # make regex pattern
    pattern="##[INFORMAT]+=<ID=" tag ",";
    if (match($0, pattern)) {
      # INFO Field
      informat=substr($0,RSTART+2,4);
      if (informat == "INFO") {
        INFO[++infoCount]=tag;
      } 
      if (informat == "FORM") {
        FORMAT[++formatCount]=tag;
      }
    } 
  }
  next;
}
/^#/ {
  # get the number and names of the samples defined in format
  for (col=9; col++<NF;) {
    sampleCount++;
    SAMPLES[++s]=$col;
  }

  IL=length(INFO);
  FL=length(FORMAT); # FL is switch for processing of FORMAT data in the read part
  # allocate the found tags to the FIELD array in order INFO --> FORMAT
  
  for (t=0; t++<length(TAGS);) {
    tag = TAGS[t];
    # look in the INFO Tags
    for (i=0;i++<IL;) {
      if (tag == INFO[i]) {
        # append Info tag to FIELD array
        FIELDS[spexCount+t] = tag;
        break;
      }
    }
    # look in the FORMAT Tags
    for (f=0;f++<FL;) {
      if (tag == FORMAT[f]) {
        for (s=0;s++<sampleCount;){
          fieldPos=spexCount+IL+f+(s-1)*FL;
          FIELDS[fieldPos]=tag;
        }
      }
    }
  }
  ############ PRINT HEADER ##########################
  printf("%s",FIELDS[1]);
  for (i=1; i++<spexCount;) {
    printf("\t%s",FIELDS[i]);
  }
  for (i=spexCount; i++< spexCount+IL;) {
    printf("\t%s",TAGNAME[FIELDS[i]]);
  }
  if (FL) {
    for (s=0; s++<sampleCount;) {
      for (i=spexCount+IL; i++<spexCount+IL+FL;) {
        line=SAMPLES[s] "-" TAGNAME[FIELDS[i]];
        printf("\t%s",line);
      }
    }
  }
  printf("\n");
  readData=1;
  firstLine=1;
}'