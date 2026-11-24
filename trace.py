import numpy as np, pandas as pd
SRC='../trace/wc98/invocation_count.csv'
def days(start='1998-06-01',end='1998-07-13',mean=1450.0):
    df=pd.read_csv(SRC,index_col=0,parse_dates=True)['count']
    df=df.groupby(level=0).sum().asfreq('min',fill_value=0)
    out=[]; names=[]
    for d in pd.date_range(start,end,freq='D'):
        x=df[d:d+pd.Timedelta('1D')-pd.Timedelta('1min')].values.astype(float)
        if len(x)!=1440: continue
        x=x/x.mean()*mean
        out.append(np.repeat(x,60)); names.append(str(d.date()))
    return out,names
