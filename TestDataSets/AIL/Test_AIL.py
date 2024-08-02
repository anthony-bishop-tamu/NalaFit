from Masfit import runMasFit
import shutil
shutil.rmtree("output")
runMasFit("kr042.txt","output", "False")