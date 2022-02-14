#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
To load paralell folders it is required to append the folder utils to
the path. This is done by this function.

Normally the path is removed after reloading the kernel, so it should
not be necessary to manually remove it. In any case, I have put the 
if statement to ensure there is no additional instances of the path.

To manually remove the path, check below.
"""

#%%
# Adding functions folder to the path, so you can access it.
import os.path
import sys
curdir = os.path.dirname(__file__)
cdir = os.path.abspath(os.path.join(curdir,'..'))

# Check if current folder in path to avoid adding many instances
if cdir not in sys.path:
    sys.path.append(cdir)
    print('Path has been appended. Package utils can be loaded')

# To remove from path
# sys.path.remove(cdir)

# To check path
# for p in sys.path:
    # print(p)
