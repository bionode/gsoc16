#!/bin/bash

# Download
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
# Extract
unzip Trimmomatic-0.36.zip
# Move files
cp Trimmomatic-0.36/trimmomatic-0.36.jar ./bin
cp -r Trimmomatic-0.36/adapters ./adapters
# Clean
rm Trimmomatic-0.36.zip
rm -rf Trimmomatic-0.36
