using DelimitedFiles

include("solver.jl")

function hilbert(m,n)
    return [1/(i+j-1) for i in 1:m, j in 1:n]
end

table2 = [(1000,1000),(500,500),(200,200),(100,100),(1000,2000),(500,2000),
(200,2000),(100,2000),(50,2000),(20,2000),(10,2000),(200,10000),(200,20000),
(200,50000),(100,20000),(100,100000),(100,200000),(50,100000),(50,200000),
(50,400000),(20,200000),(20,400000),(20,1000000),(10,100000),(10,200000),
(10,400000),(10,1000000),(10,2000000)]

for (m, n) in table2
    A = rand(m, n) .- 0.5
    z = ones(n)
    x0 = zeros(n)
    b = A * z
    local t0 = time_ns()
    y = systemSolver(A,b,x0)
    elapsed = round((time_ns()-t0)/1e9, digits=2)
    open("Table2.csv", "a") do io
        writedlm(io, [m n norm(y-x0,Inf) norm(b-A*y,Inf)/norm(b,Inf) elapsed], ',')
    end
end

table3 = [(200,200),(100,100),(50,2000),(10,2000),(200,10000),(200,20000),
(100,20000),(100,100000),(50,100000),(50,200000),(20,200000),(20,400000),
(10,100000),(10,200000),(10,400000),(10,1000000),(10,2000000)]

for (m, n) in table3
    A = hilbert(m, n)
    z = ones(n)
    x0 = zeros(n)
    b = A * z
    local t0 = time_ns()
    y = systemSolver(A,b,x0)
    elapsed = round((time_ns()-t0)/1e9, digits=2)
    open("Table3.csv", "a") do io
        writedlm(io, [m n norm(y-x0,Inf) norm(b-A*y,Inf)/norm(b,Inf) elapsed], ',')
    end
end

table4 = [(100,100),(200,200),(1000,1000),(2000,2000),(100,100000),(10,1000000)]

for (m, n) in table4
    A = rand(m, n) .- 0.5
    z = ones(n)
    x0 = zeros(n)
    b = A * z
    local t0 = time_ns()
    y = systemSolver(A,b,x0)
    elapsed = round((time_ns()-t0)/1e9, digits=2)
    open("Table4.csv", "a") do io
        writedlm(io, [m n norm(y-x0,Inf) norm(b-A*y,Inf)/norm(b,Inf) elapsed], ',')
    end
end

table5 = [(10,10),(20,20),(50,50),(100,100),(200,200),(100,100000),(10,1000000)]

for (m, n) in table5
    A = hilbert(m, n)
    z = ones(n)
    x0 = zeros(n)
    b = A * z
    local t0 = time_ns()
    y = systemSolver(A,b,x0)
    elapsed = round((time_ns()-t0)/1e9, digits=2)
    open("Table5.csv", "a") do io
        writedlm(io, [m n norm(y-x0,Inf) norm(b-A*y,Inf)/norm(b,Inf) elapsed], ',')
    end
end
