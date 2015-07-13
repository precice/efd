#!/usr/bin/env python
# encoding: utf-8
globalSourceData = 0
globalTargetData = 0
globalTime = 0.0

# This function is called first at configured timing. It can be omitted, if not
# needed. Its parameters are time, the source data, followed by the target data.
# Source and target data can be omitted (selectively or both) by not mentioning
# them in the preCICE XML configuration (see the configuration reference).
def performAction(time, sourceData, targetData):
  global globalSourceData
  global globalTargetData
  global globalTime
  globalSourceData = sourceData
  globalTargetData = targetData
  timeDiff = time - globalTime
  for i in range(globalSourceData.size):
    if timeDiff <= 0.0:
      globalTargetData[i] = 0.0
    else:
      globalTargetData[i] = (globalTargetData[i] - globalSourceData[i]) / timeDiff
      # globalTargetData[i] = (0.0 - globalSourceData[i]) / timeDiff
  globalTime = time

# This function is called for every vertex in the configured mesh. It is called
# after performAction, and can also be omitted.
def vertexCallback(id, coords, normal):
  global globalSourceData
  global globalTargetData

# This function is called at last, if not omitted.
def postAction():
  global globalSourceData
  global globalTargetData
