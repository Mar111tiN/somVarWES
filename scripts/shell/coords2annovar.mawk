#!/bin/sh

mawk '
################## DATALINES ####################
NR > 1 { # data line
  lines++;
  start=$2;
  R=$3
  A=$4
  RL=length(R);
  AL=length(A);
  ######### CHANGE INDEL COORDS ############
  # variant is an insertion/complex substitution 
  if (AL > 1) {
    if (RL == 1) { # an insertion
      # A AC
      end=start;
      # if first base is the same 
      # remove the redundant base in Ref and Alt
      if (R == substr(A,1,1)) {               
          R="-";
          A=substr(A,2);
      }
    } else { # is a complex substitution AC -> TG  
      end=start+RL-1;
    }
  # variant is a deletion
  } else if (RL > 1) {
    # is a deletion
    # make the pos of deleted bases to start and end
    start=start+1;
    end=start + RL -2;
    # remove the redundant base in Ref and Alt
    R=substr(R,2);
    A="-";
  } else { # it is a simple SNP
    end = start;
  }

  ######## OUTPUT #################
  printf("%s\t%s\t%s\t%s\t%s",$1,start,end,R,A);
  for (i=4;i++<NF;) {
    printf("\t%s",$i);
  }
  printf("\n");
  next;
}
############# HEADER ###################
NR == 1 {
  printf("%s\t%s\t%s",$1,"Start","End");
  for (i=2;i++<NF;) {
    printf("\t%s",$i);
  }
  printf("\n");
}'