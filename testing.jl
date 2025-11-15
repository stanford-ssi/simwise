using LinearAlgebra

# using Simwise.Math: hamilton_product, Quat

x0 = [0.1, 100.4, 5.0]
v0 = [1.1, 1.2, 1.3]
a = [1.0, 10.0, 100.0]/17.0
t = 1.0

println(x0 + v0*t + 0.5*a *t*t)
println(v0 + a * t)