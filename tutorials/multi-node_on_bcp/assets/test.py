import platform
try:
    import pytorch_lightning as plt
    VER=plt.__version__
except:
    VER="FAIL"
print(platform.node(), VER)