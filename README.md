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

```R
library(tidyverse)
library(corpcor) # for pseudo-inverse

df <- tribble( # data frame of origin points and Minecraft 'heading'
    ~x, ~z, ~heading,
    96, 298, 35.8,
    -304, 127, 20.2,
    -281, 1844, 109.4
  ) %>%
  mutate(
    radians = pi * heading / 180.0,
    unit_x = -sin(radians),
    unit_z = cos(radians)
  )

# Determine the 'best' point for the intersection of the lines,
# by minimizing the perpendicular distances of the point to the lines
#
# From: "Least-Squares Intersection of Lines, by Johannes Traa - UIUC 2013"
# - http://cal.cs.illinois.edu/~johannes/research/LS_line_intersect.pdf (link broken as of 2020-04)

k <- nrow(df) # number of lines
dimension <- 2
a <- df[, 1:2] %>% as.matrix() %>% t() # *columns* of origin points
n <- df[, 5:6] %>% as.matrix() %>% t() # *columns* of the points' unit direction vectors
R = matrix(data = 0, nrow = dimension, ncol = dimension) # initialize an empty matrix
q = vector(mode = 'numeric', length = dimension) # initialize an empty vector

# Generating a system of linear equations, with Rp = q, where p will be the 'best' point
for (i in 1:k) {
  R <- R + (diag(dimension) - n[, i] %*% t(n[, i]))
  q <- q + (diag(dimension) - n[, i] %*% t(n[, i])) %*% a[, i]
}
# So p_hat = pseudoinverse(R) x q
p_hat <- pseudoinverse(R) %*% q # column vector of the least squares fit best point

# Turn solution into proper data frame to plot optimal point
df_p <- t(p_hat) %>% data.frame # need to convert it to a row vector for data frame
names(df_p) <- c('x', 'z') # and match the column names to plot

df_p

#           x        z
#   -864.3211 1637.285
```
