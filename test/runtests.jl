using ChebyshevNodes
using Test

@test intoab(5,10,fromab(5,10,8))  == 8
@test intoab(5,10,fromab(5,10,5))  == 5
@test fromab(5,10,intoab(5,10,8))  == 8
@test fromab(5,10,intoab(5,10,10)) == 10
