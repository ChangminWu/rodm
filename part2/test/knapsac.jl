using JuMP
using CPLEX

include("data.dat")

m = Model(solver = CplexSolver())
@variable(m, x[i in 1:n], Bin)
@constraint(m, sum(x[i] * w[i] for i = 1:n) <= K)
@objective(m, Max, sum(x[i] * p[i] for i in 1:n))
solve(m)
