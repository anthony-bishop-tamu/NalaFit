from Masfit2.Masfit2 import runMasFit
import shutil

shutil.rmtree("../output",ignore_errors=True)
runMasFit("../gtm0375t2.txt","../output", False)