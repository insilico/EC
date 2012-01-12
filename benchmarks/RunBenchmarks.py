#! /usr/bin/python
#
# runBenchmarks.py
# Bill White - 1/10/11

import os
import sys
import random
import subprocess
import re

def makeDataset(numAttributes, numCases, numControls):
  newFileHandle = open("testData.txt", "w")
  for attributeIndex in range(1, numAttributes+1):
    newFileHandle.write("X" + str(attributeIndex) + "\t")
  newFileHandle.write("Class\n")
  for lineIndex in range(1, numCases+1):
    for attributeIndex in range(1, numAttributes+1):
      newFileHandle.write(str(random.randint(0, 2)) + "\t")
    newFileHandle.write("1\n")
  for lineIndex in range(1, numControls+1):
    for attributeIndex in range(1, numAttributes+1):
      newFileHandle.write(str(random.randint(0, 2)) + "\t")
    newFileHandle.write("0\n")
  newFileHandle.close()

def main(numIterations):
  numAttributes = 100
  numCases = 500
  numControls = 500
  print "Cases: %d, Controls: %d" % (numCases, numControls)
  print "Attributes\tTarget\tMax RAM\tElapsed Time"
  print "----------\t------------\t------------"
  for iter in range(1, int(numIterations)+1):
    # make a new data set
    makeDataset(numAttributes, numCases, numControls)
    # run the generated data set, gathering runtime stats
    numTargetAttributes = numAttributes * 0.01
    runCommands = ["ec", "--snp-data",  "testData.txt", 
                   "--ec-iter-remove-percent", "50", 
                   "--ec-num-target",  str(int(numTargetAttributes))]
    # print runCommands
    process = subprocess.Popen(runCommands, shell=False, stdout=subprocess.PIPE)
    processOutput = process.communicate()
    #print processOutput[0]
    matchGroup = re.search("EC RAM used:(.*)\\n", processOutput[0])
    if matchGroup:
        matchString = matchGroup.group().strip()
        ramMsg, ramUsed = matchString.split(':')
    matchGroup = re.search("EC elapsed time (.*)\\n", processOutput[0])
    if matchGroup:
        ecSeconds = matchGroup.group().strip()        
        ecTimeMsg, ecTimeUsed = ecSeconds.split('time')
        ecTimeUsed.strip()
    print "%10d\t%s\t%s" % (numAttributes, numTargetAttributes, ramUsed, ecTimeUsed)
    numAttributes *= 2
  # remove data set
  os.unlink("testData.txt")
  
if __name__ == "__main__":
  if len(sys.argv) != 2:
    print "Usage: " + sys.argv[0] + " [number of iterations]"
    sys.exit(1)
  else:
    main(sys.argv[1]);
