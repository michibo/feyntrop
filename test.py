
import feyntrop

#edges = [((0,1),1), ((0,1),2)]
edges = [((0,1),1), ((1,2),1), ((2,3),1), ((3,0),1)]

g = feyntrop.graph(edges)
print(g)

D = 6

#pipj = [ [  5, -5 ],
#         [ -5,  5 ] ]

pipj = [ [  0,  1,  1, -2  ],
         [  1,  0, -2,  1  ], 
         [  1, -2,  0,  1  ], 
         [ -2,  1,  1,  0  ], 
]

#pipj = [ [  .5,  .25,  1, -2  ],
#         [  .25,  0, -2,  1  ], 
#         [  1, -2,  .5,  1  ], 
#         [ -2,  1,  1,  0  ], 
#]

masses_sqr = [ 0.0, 0.0, 0.0, 0.0 ]

lamb = 1

N = 1000000

res, Itr = feyntrop.integrate_graph( g, D, pipj, masses_sqr, 6, lamb, N )

def num_err_str( number_error ):
    number, error = number_error

    return "%f +/- %f" % (number, error)

for i,term in enumerate(res):
    real, imag = term

    print( "( (%s) + I * (%s) ) * eps^%d + " % ( num_err_str(real), num_err_str(imag), i ), )

print("...")
  
    

