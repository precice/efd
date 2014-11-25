# To Do
**Priority:**
Low  
**Definition:**
PETSc produce warnings to unused arguments.  
**Date:** 24.11.14  
**Notes:**

    WARNING! There are options you set that were not used!
    WARNING! could be spelling mistake, etc!


**Priority:**
High  
**Definition:**
No velocities in volume visualization mode in Paraview.  
**Date:** 18.11.14

# In Progress

# To Verify

# Done
**Priority:**
Major  
**Definition:**
PreCICE initialization fails.  
**Date:** 24.11.14 / 25.11.14  / 25.11.14  / 25.11.14   
**Notes:**

    1416813203   08:13:23     [wscrichy]   error
    precice::utils::XMLTag::getDynVectorAttributeValue()    (0)  [PRECICE] ERROR:
    Vector attribute "value" of tag <length> has less dimensions than required (1
    instead of 2)!

Solved by using ';' instead of ',' in vector definitions.
