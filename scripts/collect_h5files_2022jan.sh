#!/bin/bash

export COPY_PATH=/net/rijn2/data2/Haoyang/ALICE/selfcal_1b1/
export PASTE_PATH=/net/rijn2/data2/Haoyang/ALICE/selected_h5/h5_collection_2022jan/

for i in {8,12,16,19,22,24,26,30,32,36,50,54,61,70,73,78,85,87,88}; do
  echo "Copying h5 file from folder "$i
  cp $COPY_PATH/$i/merged_selfcalcyle008*.h5 $PASTE_PATH
done

for i in {7,25,31,45,46,48,57,79,81}; do
  echo "Copying h5 file from folder "$i
  cp $COPY_PATH/$i/merged_selfcalcyle003*.h5 $PASTE_PATH
done

