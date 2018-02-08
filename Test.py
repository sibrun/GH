import pandas as pd

v1 = {(1,2):(2,3,None)}

df = pd.DataFrame()
df.update(v1)
print(df)

