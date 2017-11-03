#!/bin/bash

for (( separation = 94; separation <= 95; separation++ ))
do
  ./HF-gaussian $separation 5; result=$?
  if [ $result -eq 9 ];
  then
    echo "Too many loops, program skipped for separation=".$separation
    echo $separation >> data/skipped_points.txt
  elif [ $result -eq 0 ];
   then
    echo "program completed successfully for separation=".$separation
  fi

done

newfilename="data/energy_correction_data "$(date)".txt"
mv data/correction.txt "$newfilename"

$newfilename="data/skipped_points "$(date)".txt"
mv data/skipped_points.txt "$newfilename"
