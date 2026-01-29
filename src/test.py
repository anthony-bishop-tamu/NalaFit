from NalaFit.NalaFit import runNalaFit
import shutil

shutil.rmtree("../output",ignore_errors=True)
runNalaFit("../IL1Ra_gtm0375t2.txt", "../output", False)