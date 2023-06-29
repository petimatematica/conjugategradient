using LinearAlgebra
include("conjugadoDY.jl")
include("conjugado.jl")
include("testfunctionsII.jl")

x = rand(2,1)

iter1, iter2, iter3, iter4, iter5, iter6 = 0, 0, 0, 0, 0, 0

  (x1, fx1, iter1) = conjugado(x, rosenbrock, gradrosenbrock, 2, wolfe, 1);
#  (x2, fx2, iter2) = conjugado(x, powell, gradpowell, 2, armijo, 1);
#  (x3, fx3, iter3) = conjugado(x, powell, gradpowell, 2, goldstein, 0);
#  (x4, fx4, iter4) = conjugado(x, powell, gradpowell, 2, goldstein, 1);
#  (x5, fx5, iter5) = conjugado(x, powell, gradpowell, 2, wolfe, 0);
#  (x6, fx6, iter6) = ([0; 1], powell, gradpowell, 2, armijo, 1, maxiter = 1000);

# v = [iter1; iter2; iter3; iter4; iter5; iter6]

n = 100;
p = Float64[]
for j in 1:n
    push!(p, 1-(j/n))
end

x = rand(2)

# x1, fx1, k, counter = conjugadoPRP(x, rosenbrock, gradrosenbrock, 2);
 x, fx, k, counter = conjugadoPRP(p, vardim, gradvardim, n);

