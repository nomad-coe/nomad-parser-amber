import sys
import os
import os.path
basedir = os.path.dirname(os.path.abspath(__file__))
commondir = os.path.normpath(os.path.join(basedir,"../../../../python-common/common/python"))

if not commondir in sys.path:
    sys.path.insert(0, commondir)
