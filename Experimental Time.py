'''
time = time of interest; can be a single time point, panda dataframe, or list
fmt = format of the list, refer to datetime string formatting to determine appropriate formatting
'''

from datetime import datetime
import pandas as pd

def Exp_Time(time, fmt='%m/%d/%Y %H:%M:%S'):
    j = []
    for i in time:
        x=datetime.strptime(i, fmt)
        x=x.timestamp()
        j.append(x)
    j[:] = [k - j[0] for k in j]
    y = {'Experimental Time': j}
    z = pd.DataFrame(y)
    return(z)