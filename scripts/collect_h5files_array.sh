#!/bin/bash

export COPY_PATH=/net/rijn2/data2/Haoyang/ALICE/selfcal_array/
export PASTE_PATH=/net/rijn2/data2/Haoyang/ALICE/selected_h5/h5_collection_array/

echo "Coping h5 files from selfcal_array and save them to folder h5_collection_array"

for i in {8,12,16,19,22,24,26,30,32,36,43,47,50,54,61,70,73,78,85,87,88}; do
  echo "Copying h5 file from folder "$i
  cp $COPY_PATH/$i/merged_selfcalcyle008*.h5 $PASTE_PATH
done

for i in {7,14,25,31,45,46,48,57,79,81}; do
  echo "Copying h5 file from folder "$i
  cp $COPY_PATH/$i/merged_selfcalcyle003*.h5 $PASTE_PATH
done

#Not so good, but still better than no self-calibration, mostly visually
for i in {4,9,23,27,29,34,35,37,49,53,56,59,66,71,76,80,82,84,86,90}; do
  echo "Copying h5 file from folder "$i
  cp $COPY_PATH/$i/merged_selfcalcyle008*.h5 $PASTE_PATH
done

for i in {51,83}; do
  echo "Copying h5 file from folder "$i
  cp $COPY_PATH/$i/merged_selfcalcyle003*.h5 $PASTE_PATH
done
