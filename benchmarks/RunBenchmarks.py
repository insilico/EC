#! /usr/bin/python
#
# runBenchmarks.py
# Bill White - 1/10/11

import os
import sys
import random
import subprocess
import re
import time

def makeDataset(numAttributes, numCases, numControls):
  newFileName = "testData_" + str(numCases) + "_" + str(numControls) + "_" + str(numAttributes) + ".txt"
  newFileHandle = open(newFileName, "w")
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

def main(numIterations, numCases, numControls, numAttributes, ompChunkSize):
  print "Cases: %d, Controls: %d, Chunk size: %s" % (numCases, numControls, ompChunkSize)
  print "Attributes Target Max RAM     Elapsed Time Wall Time"
  print "---------- ------ ----------- ------------ ---------"
  for iter in range(1, int(numIterations)+1):
    # make a new data set
    makeDataset(numAttributes, numCases, numControls)
    # run the generated data set, gathering runtime stats
    numTargetAttributes = numAttributes * 0.1
    newFileName = "testData_" + str(numCases) + "_" + str(numControls) + "_" + str(numAttributes) + ".txt"
    runCommands = ["ec", "--snp-data",  newFileName,
                   "--ec-iter-remove-percent", "50",
                   "--ec-num-target",  str(int(numTargetAttributes))]
    os.putenv("OMP_SCHEDULE", "dynamic, " + ompChunkSize)
    # print runCommands
    t0 = time.time()
    process = subprocess.Popen(runCommands, shell=False, stdout=subprocess.PIPE)
    processOutput = process.communicate()
    tdone = time.time()
    #print processOutput[0]
    matchGroup = re.search("EC RAM used:(.*)\\n", processOutput[0])
    if matchGroup:
        matchString = matchGroup.group().strip()
        ramMsg, ramUsed = matchString.split(':')
    else:
      print "Error matching EC RAM used"
      sys.exit(1)
    matchGroup = re.search("EC elapsed time (.*)\\n", processOutput[0])
    if matchGroup:
        ecSeconds = matchGroup.group().strip()
        ecTimeMsg, ecTimeUsed = ecSeconds.split('time')
        ecTimeUsed.strip()
    else:
      print "Error matching EC elapsed time"
      sys.exit(1)
    print "%9d %6d %11s %12s %9d" % (numAttributes, numTargetAttributes, ramUsed, ecTimeUsed, tdone - t0)
    numAttributes *= 2
  # remove data set
  # os.unlink("testData.txt")

if __name__ == "__main__":
  if len(sys.argv) != 6:
    print "Usage: " + sys.argv[0] + " [number of iterations] [number of cases] [number of controls] [number of starting attributes] [OpenMP chunk size]"
    sys.exit(1)
  else:
    main(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]);
