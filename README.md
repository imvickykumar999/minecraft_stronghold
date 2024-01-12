# [`Minecraft stronghold`](https://imvickykumar999.github.io/minecraft_stronghold/)

<tr>
<td>Facing <i>direction</i> (Towards <i>axis</i>) (<strong>yaw</strong>/<strong>pitch</strong>)
</td>
<td>The orientation of the player:
<br><br>
<ul><li><i>direction</i> is the cardinal direction the player is facing;</li>
<li><i>axis</i> is the horizontal axis and the direction along this axis that the player is tooking towards ("negative Z", "positive Z", "negative X", or "positive X");</li>
<br>
<li><i>yaw</i> is the player's <a href="https://en.wikipedia.org/wiki/Aircraft_principal_axes#Vertical_axis_(yaw)" class="extiw" title="wikipedia:Aircraft principal axes">yaw</a>, 0 meaning facing South (positive is to the West, negative is to the East);</li>
<li><i>pitch</i> is the player's <a href="https://en.wikipedia.org/wiki/Aircraft_principal_axes#Transverse_axis_(pitch)" class="extiw" title="wikipedia:Aircraft principal axes">pitch</a>, 0 meaning looking horizontally (positive is down, negative is up).</li></ul>
</td>
</tr>

<br>
<table>
<tr>
<td>Roll</td>
<td>Pitch</td>
<td>Yaw</td>
</tr>
<tr>
<td><img src="https://upload.wikimedia.org/wikipedia/commons/thumb/c/cc/Aileron_roll.gif/223px-Aileron_roll.gif" alt="3" height = '100%'></td>
<td><img src="https://upload.wikimedia.org/wikipedia/commons/thumb/e/ec/Aileron_pitch.gif/224px-Aileron_pitch.gif" alt="3" height = '100%'></td>
<td><img src="https://upload.wikimedia.org/wikipedia/commons/thumb/9/96/Aileron_yaw.gif/225px-Aileron_yaw.gif" alt="3" height = '100%'></td>
</tr>
</table>
<br>

    Locate Minecraft stronghold from Eye of Ender throws
    
![image](https://github.com/imvickykumar999/minecraft_stronghold/assets/50515418/bf7dc07b-e7c7-46e3-9952-cada71fe6449)
![image](https://github.com/imvickykumar999/minecraft_stronghold/assets/50515418/47b0a8f3-fe00-454b-92af-3990a5e9afed)
![image](https://github.com/imvickykumar999/minecraft_stronghold/assets/50515418/320527c8-1802-4115-950a-09c80da36f29)

    Eye of Ender Throws: (X, Z, Yaw)
    Yaw is Direction in XZ plane.

![image](https://github.com/imvickykumar999/minecraft_stronghold/assets/50515418/861c3aa4-083b-469e-a777-587bf5931fe8)
![image](https://github.com/imvickykumar999/minecraft_stronghold/assets/50515418/9e48647a-1325-40f3-b8de-b5f2ea556b8c)

    Locate Minecraft stronghold from Eye of Ender throws by least-square intersection of lines

![image](https://github.com/imvickykumar999/minecraft_stronghold/assets/50515418/e8b9cbea-39d2-40a0-b85f-1de40c8c7c0e)

<!--
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
-->

## Implement using 3D [`vixtor.py`](https://github.com/imvickykumar999/100th-Repository-First-Milestones/blob/main/PYPI%20python%20package/multivicks/vixtor.py#L108) python library.

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

## `Reading` `X, Y, Z`, `Roll`, `Yaw` ***and*** `Pitch`.
    
```python
import mcpi.minecraft as minecraft
mc = minecraft.Minecraft.create()

# Get the player's position
pos = mc.player.getTilePos()

# Print the player's x,y,z
print("X:", pos.x)
print("Y:", pos.y)
print("Z:", pos.z)

# Get the player's rotation
rot = mc.player.getRotation()

# Print the player's roll, yaw and pitch
print("Roll:", rot.roll)
print("Yaw:", rot.yaw)
print("Pitch:", rot.pitch)

# We need only X, Z, and Yaw
print(pos.x, pos.z, rot.yaw)
```

### `Minecraft` uses angle as `Degree`, while `Numpy` uses it as `Radian`.

```python
import numpy as np

def convert_minecraft_to_numpy(yaw, pitch):
  """Converts Minecraft yaw and pitch to a numpy array.

  Args:
    yaw: The Minecraft yaw angle in degrees.
    pitch: The Minecraft pitch angle in degrees.

  Returns:
    A numpy array containing the yaw and pitch angles in radians.
  """

  yaw = yaw * np.pi / 180.0
  pitch = pitch * np.pi / 180.0
  return np.array([yaw, pitch])

def convert_numpy_to_minecraft(yaw, pitch):
  """Converts numpy yaw and pitch to Minecraft angles.

  Args:
    yaw: The numpy yaw angle in radians.
    pitch: The numpy pitch angle in radians.

  Returns:
    A tuple containing the yaw and pitch angles in degrees.
  """

  yaw = yaw * 180.0 / np.pi
  pitch = pitch * 180.0 / np.pi
  return yaw, pitch
```
