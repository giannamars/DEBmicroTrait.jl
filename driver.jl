using Pkg
Pkg.activate(".")

include("test/test_output.jl")

for arg in ARGS
  println(arg)
end
