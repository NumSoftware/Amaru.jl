using Amaru

lc = 0.0025
r  = 0.025   # radius
b  = 0.002
h  = sqrt(r^2-b^2) # contact length


# 2D
coords = [ 
          0.0 0.0 lc
          -r  0.0 lc
          +r  0.0 lc
           b  h   lc/4
          -b  h   lc/5
         ]

lines = [ 
         [2, 3],
         [3, 1, 4],
         [4, 5],
         [5, 1, 2],
        ]

#mesh = UMesh(coords, lines)

# 3D
coords = [ -1.0 -1.0 -1.0 0.2
            1.0 -1.0 -1.0 0.2
            1.0  1.0 -1.0 0.2
           -1.0  1.0 -1.0 0.2
           -1.0 -1.0  1.0 0.2
            1.0 -1.0  1.0 0.2
            1.0  1.0  1.0 0.2
           -1.0  1.0  1.0 0.2 ]

lines = [ 
         [1, 2], #1
         [2, 3], #2
         [3, 4], #3
         [4, 1], #4
         [5, 6], #5
         [6, 7], #6
         [7, 8], #7
         [8, 5], #8
         [1, 5], #9
         [2, 6], #10
         [3, 7], #11
         [4, 8], #12
        ]

faces = [
            [ 1, 10, 5, 9 ],
            [ 3, 12, 7, 11 ],
            [ 4, 9, 8, 12 ], 
            [ 2, 11, 6, 10 ],
            [ 1, 4, 3, 2 ],
            [ 5, 6, 7, 8 ]
           ]

#mesh = UMesh(coords, lines, faces)
#save(mesh, "unstructured.vtu")

poly = Polygon([0 0; 1 0; 0 1])
polym = extrude(poly)
#coords = Amaru.get_coords(polym)
#lines = Amaru.line_indexes(polym)
#faces = Amaru.poly_indexes(polym)
#coords = [ coords 0.1*ones(6) ]


poly = Polygon([0 0; 1 0; 1 1; 0 1])
polym = extrude(poly)

tag!(polym.points[:(x==0)], "left") 
tag!(polym.points[:(x==1)], "right") 

sizes = [
         "left" => 0.1
         "right" => 0.05
        ]

embed = [
         :(z==0) => [0.5, 0.5, 0.0]
        ]

mesh = UMesh(polym.polygons[1], sizes, embed, size=0.1)
#mesh = UMesh(polym, sizes, embed)

#mesh = UMesh(coords, lines, faces)
# elemsize, sizes=[1=>0.1, 2=>1.1]
# tag!(polym.points[:(x==5)], "tag1")
# sizes = [
#          "tag1" => 0.1
#          "tag2" => 0.2
#          ]
# embed = [
#          "facetag" => [ [], [] ]
#         ]
save(mesh, "unstructured.vtu")
