# %%
import pandas as pd
import seaborn as sns
import numpy as np

# %%
# Read and prepare data
df = pd.read_csv("log.csv")

df.loc[df['matrix'] == 'matrices/A_n1e6_d4.mtx', 'matrix'] = 'n=1e6,d=4'
df.loc[df['matrix'] == 'matrices/A_n2e6_d4.mtx', 'matrix'] = 'n=2e6,d=4'
df.loc[df['matrix'] == 'matrices/A_n3e6_d4.mtx', 'matrix'] = 'n=3e6,d=4'
df.loc[df['matrix'] == 'matrices/A_n5e6_d4.mtx', 'matrix'] = 'n=5e6,d=4'
# print(df)
# df.matrix[df.matrix == 'matrices/A_n1e6_d4.mtx'] = 'n=1e6,d=4'

# %%
type = 'hybrid'

filt = df[(df['type'] == type)]
time = filt[filt['matrix'] == 'n=1e6,d=4']
time = time.sort_values(by=['time'], ascending=True)
time = time.groupby('tasks')['time'].mean()
base = time.values[0]
print(time)
speedup = base / time.values
print(speedup)
filt = filt[['tasks', 'time', 'matrix']]
h = sns.lineplot(data=filt, x='tasks', y='time',
                 style='matrix', hue='matrix', markers=True, dashes=False)
h.set_xticks(np.unique(filt['tasks'].values))
h.set_xlabel("Threads")
h.set_ylabel("Time (s)")

if type == 'hybrid':
    h.set_xticklabels(['1x1','2x2','3x3','4x4','5x5','6x6'])

fig = h.get_figure()
fig.savefig("../report/assets/" + str(type) + ".png", dpi=300)
# %%
# SPEEDUP

matrices = ['n=1e6,d=4', 'n=2e6,d=4', 'n=3e6,d=4', 'n=5e6,d=4']
type = 'hybrid'

filt = df[df['type'] == type]
filt = filt[['tasks', 'time', 'matrix']]
for i in matrices:
    base = filt[(filt['matrix'] == i) & (filt['tasks'] == 1)]
    base = np.mean(base['time'])
    filt.loc[filt['matrix'] == i, 'time'] = base / \
        filt.loc[filt['matrix'] == i, 'time']

#time = filt[filt['matrix'] == 'n=1e6,d=4']
#time = time.sort_values(by=['time'], ascending=True)
#print(time)

h = sns.lineplot(data=filt, x='tasks', y='time', style='matrix', hue='matrix', markers=True, dashes=False)
h.set_xticks(np.unique(filt['tasks'].values))
h.set_xlabel("Processes x Threads")
h.set_ylabel("Speedup")

if type == 'hybrid':
    h.set_xticklabels(['1x1','2x2','3x3','4x4','5x5','6x6'])

fig = h.get_figure()
fig.savefig("../report/assets/" + str(type) + "_speedup" + ".png", dpi=300)

# %%
type = 'openmpi'

#filt = df[(df['type'] == type) & (df['tasks'] >= 4) | (df['tasks'] == 8) | (df['tasks'] == 12) | (df['tasks'] == 16) | df['tasks'] == 20]
filt = df[(df['type'] == type) & (df['tasks'] == 4) | ((df['tasks'] == 8) | (df['tasks'] == 12) | (df['tasks'] == 16) | (df['tasks'] == 20))]
filt = filt[['tasks', 'time', 'matrix']]
h = sns.lineplot(data=filt, x='tasks', y='time',
                 style='matrix', hue='matrix', markers=True, dashes=False)
h.set_xticks(np.unique(filt['tasks'].values))
h.set_xlabel("Threads")
h.set_ylabel("Time (s)")

fig = h.get_figure()
fig.savefig("../report/assets/" + str(type) + "_special.png", dpi=300)

# %%
# SPEEDUP

matrices = ['n=1e6,d=4', 'n=2e6,d=4', 'n=3e6,d=4', 'n=5e6,d=4']
type = 'openmpi'

filt = df[(df['type'] == type) & (df['tasks'] == 4) | ((df['tasks'] == 8) | (df['tasks'] == 12) | (df['tasks'] == 16) | (df['tasks'] == 20))]
filt = filt[['tasks', 'time', 'matrix']]
for i in matrices:
    base = filt[(filt['matrix'] == i) & (filt['tasks'] == 4)]
    base = np.mean(base['time'])
    filt.loc[filt['matrix'] == i, 'time'] = base / \
        filt.loc[filt['matrix'] == i, 'time']

h = sns.lineplot(data=filt, x='tasks', y='time', style='matrix', hue='matrix', markers=True, dashes=False)
h.set_xticks(np.unique(filt['tasks'].values))
h.set_xlabel("Processes x Threads")
h.set_ylabel("Speedup")

fig = h.get_figure()
fig.savefig("../report/assets/" + str(type) + "_speedup" + "_special.png", dpi=300)