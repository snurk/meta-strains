

import sys
import pandas as pd

df = pd.read_csv(sys.argv[1], index_col=0)
df.index = df.index.astype(str).str[6:].astype(int)
df = df.sort_index()
df = df.transpose()
df = df.round(3)

df.to_csv(sys.argv[2], header=None)

