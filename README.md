# `Minecraft stronghold`

![image](https://github.com/imvickykumar999/minecraft_stronghold/assets/50515418/bf7dc07b-e7c7-46e3-9952-cada71fe6449)
![image](https://github.com/imvickykumar999/minecraft_stronghold/assets/50515418/47b0a8f3-fe00-454b-92af-3990a5e9afed)
![image](https://github.com/imvickykumar999/minecraft_stronghold/assets/50515418/320527c8-1802-4115-950a-09c80da36f29)
![image](https://github.com/imvickykumar999/minecraft_stronghold/assets/50515418/861c3aa4-083b-469e-a777-587bf5931fe8)
![image](https://github.com/imvickykumar999/minecraft_stronghold/assets/50515418/9e48647a-1325-40f3-b8de-b5f2ea556b8c)

    Locate Minecraft stronghold from Eye of Ender throws by least-square intersection of lines

![image](https://github.com/imvickykumar999/minecraft_stronghold/assets/50515418/e8b9cbea-39d2-40a0-b85f-1de40c8c7c0e)

Repo includes:

* stronghold.R
  * R script where coordinates and headings can be manually edited in the code
  * run in RStudio console
  * df_p will contain the [x,z] coordinates for the least-squares solution for the stronhold
  
* app.R (in /stronhold_triangulator/)
  * R shiny app for interactive calculations
  * can open app.R in RStudio and then Run App for interactive use
  * enter rows of X coordinate, Z coordinate, heading and click Submit button
  * will plot the rays, least-sqaures solution, and text annotation of the solution coordinates

```python
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
```
