from Masfit2.Masfit2 import runMasFit
import shutil

shutil.rmtree("../output",ignore_errors=True)
runMasFit("../IL1Ra_gtm0375t2.txt","../output", False)