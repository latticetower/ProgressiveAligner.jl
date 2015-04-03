using ProgressiveAligner
using Base.Test

# write your own tests here
@test 1 == 1
tests = [
    "DataReaderTests",
    #"DataWriterTests",
    "ClusteringTests",
    "ProfileAlignerTests"
    ]

println("Running tests:")
for t in tests
    println(" * $(t)")
    include("$(t).jl")
end
