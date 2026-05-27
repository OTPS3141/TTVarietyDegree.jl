using TTVarietyDegree

D = [1, 2, 2, 2, 1]
d = [3, 2, 2, 2]

P = compute_tail_polynomial(D, d)
fP = compute_fP(P, D[end - 1], d[end])
deg = degree_TT_variety(D, d)

println("D = ", D)
println("d = ", d)
println("P(CC*) = ", P)
println("f(P) = ", fP)
println("degree = ", deg)
