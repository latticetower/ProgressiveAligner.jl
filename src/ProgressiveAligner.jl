module ProgressiveAligner
  push!(LOAD_PATH, dirname(@__FILE__()))

  #export DataReader,
  #      DataWriter,
  #      ProfileAligner,
  #      Clustering

  include("DataReader.jl")
  include("DataWriter.jl")
  include("ProfileAligner.jl")
  include("Clustering.jl")
end # module
