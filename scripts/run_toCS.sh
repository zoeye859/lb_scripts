#!/bin/bash

filename="merged_selfcalcyle003_ILTJ160335.24+540514.4_concat.ms.avg.h5
merged_selfcalcyle003_ILTJ160511.12+553857.0_concat.ms.avg.h5
merged_selfcalcyle003_ILTJ160622.16+541336.7_concat.ms.avg.h5
merged_selfcalcyle003_ILTJ160744.52+545022.0_concat.ms.avg.h5
merged_selfcalcyle003_ILTJ161002.77+555243.0_concat.ms.avg.h5
merged_selfcalcyle003_ILTJ161057.72+553527.9_concat.ms.avg.h5
merged_selfcalcyle003_ILTJ161116.64+534419.1_concat.ms.avg.h5
merged_selfcalcyle003_ILTJ161129.74+532709.2_concat.ms.avg.h5
merged_selfcalcyle003_ILTJ161228.96+550637.3_concat.ms.avg.h5
merged_selfcalcyle003_ILTJ161716.32+551546.9_concat.ms.avg.h5
merged_selfcalcyle003_ILTJ161815.18+555350.9_concat.ms.avg.h5
merged_selfcalcyle003_ILTJ161834.06+551751.2_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ160232.60+545406.6_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ160421.92+550545.3_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ160435.15+535938.4_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ160454.73+555949.7_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ160538.33+543922.6_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ160559.99+545405.5_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ160606.04+555400.0_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ160607.60+552135.5_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ160609.66+540325.5_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ160623.53+540555.9_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ160633.74+535616.5_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ160702.74+534128.6_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ160721.60+534639.2_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ160744.79+554841.5_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ160820.72+561355.7_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ160828.33+541031.0_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ160834.19+553241.4_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ160835.65+533822.0_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ160958.90+550026.6_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161100.30+544204.1_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161119.68+552846.4_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161129.52+544606.6_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161142.05+561952.0_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161212.30+552303.6_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161224.83+555437.1_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161240.15+533558.3_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161327.82+555158.0_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161406.23+560926.4_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161507.56+554540.7_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161537.86+534646.4_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161623.37+545738.5_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161640.37+535813.3_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161710.15+551746.6_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161724.13+555610.7_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161832.97+543146.3_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161859.41+545246.3_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161900.61+542936.8_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161919.70+553556.7_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161920.16+544817.4_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ161929.70+551315.3_concat.ms.avg.h5
merged_selfcalcyle008_ILTJ162055.04+544854.6_concat.ms.avg.h5"


for i in $filename; do
  echo "Processing file "$i
  python gains_toCS_h5parm.py ../$i '../../h5_collection/sub6asec_L686962_SB001_uv_12CFFDDA8t_121MHz.ms.sub.shift.avg.apply_infield.avg4ch.ms'

done
