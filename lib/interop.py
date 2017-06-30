import os

def is_unix():
    return( os.name == "posix" )

def is_windows():
    return( os.name == "nt" )
    