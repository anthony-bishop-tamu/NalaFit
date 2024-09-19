from Masfit import runMasFit
import shutil

shutil.rmtree("output",ignore_errors=True)
runMasFit("gtm0375t2.txt","output", "False")