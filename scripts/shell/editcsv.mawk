#!/bin/sh

mawk '
#################################################
################## DATALINES ####################
NR > 1 { # data line
  ##### INCOLUMNS #######
  for (col=0;col++<NF;) {
    DATA[col]=$col;
  }
  ##### EXPAND DATA #####
  for (xp in XP) {
    col=XP[xp];
    split($col,XPDATA,",");
    for (d in XPDATA) {
      DATA[col "-" d] = XPDATA[d];
    }
  }
  
  ###### TRANSLATE ######
  for (tr in TR) {
    if (TRANS[tr ":" DATA[tr]]) {
      DATA[tr] = TRANS[tr ":" DATA[tr]];
    }
  }

  ####### OUTPUT #############
  # first line
  printf("%s",DATA[COLPOS[1]]);
  for (col=1;col++<pointer;) {
    printf("\t%s",DATA[COLPOS[col]]);
  }
  printf("\n");
  next;
}
NR == 1 {
  ######## EDITS #############
  substitute="TUMOR-:T,NORMAL-:N,Status:_status";
  expand="TreadCounts>>TR1+,-,TR2+,-;NreadCounts>>NR1+,-,NR2+,-";
  outCols="=>,somatic_status,TR1,TR1+,TR2,TR2+,NR1,NR1+,NR2,NR2+,somaticP,variantP";
  translate="somatic_status>>1:Germline,2:Somatic,3:LOH";

  ######## SUBSTITUTION in the header
  header=$0;
  split(substitute,SUB,",");
  for (s in SUB) {
    split(SUB[s],SUBB,":");
    gsub(SUBB[1],SUBB[2]);
  }

  # get the original col order into COLPOS
  for (col=0;col++<NF;) {
    INCOLS[$col]=col;
    # INCOLS["Chr"] = 1;
    # ...
    HEADER[col]=$col;
  }

  ######### EXPAND ######################
  split(expand,XP1,";");
  for (xpnd in XP1) { # split each entry
    split(XP1[xpnd],XP2,">>");
    # get the col position of the expand
    col = INCOLS[XP2[1]];
    # make XP array containing the cols that need expanding
    XP[xpnd]=col;
    # XP[1] = 4
    # XP[2] = 8
    # the 4th and 8th will be expanded

    xpColCount=split(XP2[2],XPCOLS,",");
    # make the association array for the expand fields
    for (c=0;c++<xpColCount;) {
      EXP[XPCOLS[c]] = col "-" c;
      # EXP[TR2+] = 4-3;

      # also populate the HEADER array for the header output
      HEADER[col "-" c] = XPCOLS[c];
      # HEADER[4-3] = "TR2+"
    }
  }
  ########## DBG ###########
  # for (xp in EXP) {
  #   print(xp, EXP[xp]);
  # }
  ##########################

  ####### OUTCOLS ##########################
  # outCols --> FIELDS array
  fieldCount = split(outCols,FIELDS,",");

  # cycle through FIELDS and assign the original input cols from ORGCOLS
  pointer=0;
  for (col=0;col++<fieldCount;) {
    pointer++;
    # for =>, fill up the COLPOS array from the next FIELD element
    if (FIELDS[col] == "=>") {
      fillUp=1;
      col++; # skip to next outCol
    } # normal field
    if (FIELDS[col] in INCOLS) {
      # assign to OUTCOLS the proper col position
      COLPOS[pointer] = INCOLS[FIELDS[col]];
      # COLPOS[INCOLS["somatic_status"]] = COLPOS[6] = 6
    } else { # is an expanded field
      COLPOS[pointer] = EXP[FIELDS[col]];
    }
    # for "=>" case 
    if (fillUp) { # fill up the previous cols
      # get the position also from expand position
      pos = COLPOS[pointer];
      if (match(pos, "-")) {
        pointer=substr(pos,1,RSTART-1);
      } else {
        pointer=pos
      }
      for (i=0; i++<pointer;) {
        COLPOS[i] = i;
        # COLPOS[1]=1
        # COLPOS[6]=6
      }
      fillUp=0;
    }
  }
  ########## DBG ###########
  # for (c=0; c++<length(COLPOS);) {
  #   print("COLPOS", c, COLPOS[c]);
  # }
  # printf("\n");
  ##########################

  ################ TRANSLATE #####################
  split(translate,TRANS1,";");
  # TRANS1[1] = somatic_status>>1:Germline,2:Somatic,3:LOH
  for (trans in TRANS1) {
    split(TRANS1[trans],TRANS2,">>");
    # TRANS2[1] = somatic_status
    # TRANS2[2] = 1:Germline,2:Somatic,3:LOH
    # lookup the right colpos of the datafield (incl. expanded)
    for (col=0;col++<length(HEADER);) {
      if (HEADER[COLPOS[col]] == TRANS2[1]) {
        col = COLPOS[col];
        # TR array for identifying transable columns
        TR[col] = col;
        # col = 1-4
        break;
      } 
    }
    split(TRANS2[2],TRANS3,",");
    # TRANS3[1] = "1:Germline"
    for (trans in TRANS3){
      split(TRANS3[trans],TRANS4,":");
      TRANS[col ":" TRANS4[1]] = TRANS4[2];
    } 
  }

  ######## HEADER OUTPUT
  # first line
  printf("%s",HEADER[COLPOS[1]]);
  for (col=1;col++<pointer;) {
    printf("\t%s",HEADER[COLPOS[col]]);
  }
  printf("\n");
}'