# https://www.online-python.com/E8c6nwAQP2
# https://chat.openai.com/share/c17ff428-1155-4a09-9671-f368ef252609
# https://www.incidentalfindings.org/posts/2020-04-10_finding-minecraft-strongholds-with-math/#math

import numpy as np
import pandas as pd

# Create DataFrame
df = pd.DataFrame({
    'x': [96, -304, -281],
    'z': [298, 127, 1844],
    'heading': [35.8, 20.2, 109.4]
})

# Convert degrees to radians and calculate unit vectors
df['radians'] = np.radians(df['heading'])
df['unit_x'] = -np.sin(df['radians'])
df['unit_z'] = np.cos(df['radians'])

# Number of lines and dimension
k = df.shape[0]
dimension = 2

# Extract matrices a and n
a = df[['x', 'z']].values.T
n = df[['unit_x', 'unit_z']].values.T

# Initialize matrices R and q
R = np.zeros((dimension, dimension))
q = np.zeros(dimension)

# Generating a system of linear equations, with Rp = q
for i in range(k):
    R += np.eye(dimension) - np.outer(n[:, i], n[:, i].T)
    q += (np.eye(dimension) - np.outer(n[:, i], n[:, i].T)) @ a[:, i]

# Calculate the least squares fit best point
p_hat = np.linalg.pinv(R) @ q

# Convert the solution to a DataFrame for plotting
df_p = pd.DataFrame(p_hat.reshape(1, -1), columns=['x', 'z'])

print(df_p)
