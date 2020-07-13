using Amaru

lc = 0.0025
r  = 0.025   # radius
b  = 0.002
h  = sqrt(r^2-b^2) # contact length

coords = [ 
          0.0 0.0 lc
          -r  0.0 lc
          +r  0.0 lc
           b  h   lc/4
          -b  h   lc/5
         ]

conns = [ 
         [2, 3],
         [3, 1, 4],
         [4, 5],
         [5, 1, 2],
        ]

mesh = UMesh(coords, conns)

