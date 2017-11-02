#!/bin/bash

for (( separation = 35; separation <= 36; separation++ ))
do
   ./HF $separation
   echo "program completed for separation=".$separation
done

newfilename="data/energy_correction_data "$(date)".txt"

mv data/correction.txt "$newfilename"
