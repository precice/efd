#!/usr/bin/env python
#!/us# encoding: utf-8
mySourceData = 0
myTargetData = 0

#
# This function is called first at configured timing. It can be omitted, if not 
# needed. Its parameters are time, the source data, followed by the target data. 
# Source and target data can be omitted (selectively or both) by not mentioning 
# them in the preCICE XML configuration (see the configuration reference).
#
def performAction(time, sourceData, targetData):
  # Usage example 1:
  global mySourceData
  global myTargetData
  mySourceData = sourceData # store (reference to) sourceData for later use
  myTargetData = targetData # store (reference to) targetData for later use
  for i in range(mySourceData.size):
    myTargetData[i] = myTargetData[i] - mySourceData[i]
    
#
# This function is called for every vertex in the configured mesh. It is called
# after performAction, and can also be omitted.
#
def vertexCallback(id, coords, normal):
  # Usage example:
  global mySourceData # Make global data set in performAction visible
  global myTargetData
  # myTargetData[id] += coords[0] + mySourceData[id] # Add data to vertex coords
    
#
# This function is called at last, if not omitted.
#
def postAction():
  global mySourceData # Make global data set in performAction visible
  global myTargetData
  # Do something ..
