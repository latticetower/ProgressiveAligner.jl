module Clustering
  export upgma, wpgma
  import ProfileAligner.Profile

  typealias DistanceHash{T} Dict{Profile{T}, Dict{Profile{T}, T}}

  function getUPGMAWeight{T}(distancesInfo :: Dict{Profile{T}, Dict{Profile{T}, T}},
                          k :: Profile{T}, f :: Profile{T}, g :: Profile{T})
    (f.numberofstrings*get(distancesInfo, f, k) + g.numberofstrings*get(distancesInfo, g, k))/(g.numberofstrings + f.numberofstrings)
  end

  function getWPGMAWeight{T}(distancesInfo :: Dict{Profile{T}, Dict{Profile{T}, T}}, k :: Profile{T}, f :: Profile{T}, g :: Profile{T})
    (get(distancesInfo, f, k) + get(distancesInfo, g, k))/(one(T) + one(T))
  end
  
  upgma{T}(startvertices::Vector{Profile{T}}, scoreFunc, mergeFunc) = buildTree(startvertices, scoreFunc, mergeFunc, getUPGMAWeight)
  wpgma{T}(startvertices::Vector{Profile{T}}, scoreFunc, mergeFunc) = buildTree(startvertices, scoreFunc, mergeFunc, getWPGMAWeight)


  function buildTree{T}(startvertices :: Vector{Profile{T}}, scoreFunc, mergeFunc, weightFunc)
    distancesInfo = initDistanceMatrix(startvertices, scoreFunc)
    while(length(startvertices) > 1)
      simplifyTree!(distancesInfo, startvertices, mergeFunc, weightFunc)
    end
    #println(length(startvertices))
    startvertices[1]
  end

  function simplifyTree!{T}(distancesInfo :: Dict{Profile{T}, Dict{Profile{T}, T}},
                               vertices :: Vector{Profile{T}}, mergeFunc, weightFunc)
    n = length(vertices)
    Q = Array(T, n, n)
    K1 = zero(Array(T, n))
    K2 = zero(Array(T, n))
    #println("in simplify")
    #1. calculate Q

    (mini, minj) = (1, 2)
    for i = 1 : n - 1
      for j = i + 1 : n
        #println(vertices[i], vertices[j])
        Q[i, j] = get(distancesInfo, vertices[i], vertices[j])
        if Q[mini, minj] > Q[i, j]
          (mini, minj) = (i, j)
        end
      end
    end
    #
    newT = mergeFunc(vertices[mini], vertices[minj])
    computeDistances!(distancesInfo, vertices, mini, minj, newT, weightFunc)
    #
    deleteat!(vertices, [mini, minj])
    push!(vertices, newT)
  end

  function computeDistances!{T}(distancesInfo :: Dict{Profile{T}, Dict{Profile{T}, T}},
                               vertices :: Vector{Profile{T}}, mini :: Int64, minj :: Int64,
                               newT :: Profile{T}, weightFunc)
    #compute new distances from new to previous 2
    distancesInfo[newT] = Dict{Profile{T}, T}()
    f = vertices[mini]
    g = vertices[minj]
    t2 = one(T) + one(T)
    distancesInfo[newT][f] = get(distancesInfo, f, g) / t2
    distancesInfo[newT][g] = get(distancesInfo, f, g) - distancesInfo[newT][f]
    distancesInfo[newT][newT] = zero(T) #assume that profile aligns to profile with zero distance
    #compute to all others
    for i in 1 : length(vertices)
      i == mini && continue
      i == minj && continue
      k = vertices[i]
      distancesInfo[newT][k] = weightFunc(distancesInfo, k, f, g)
    end

  end

  function initDistanceMatrix{T}(startvertices :: Vector{Profile{T}}, ScoreFunc)
    distancesInfo = Dict{Profile{T}, Dict{Profile{T}, T}}()
    for vi in 1 : length(startvertices)
      for vj in vi : length(startvertices)
        if !haskey(distancesInfo, startvertices[vi])
          distancesInfo[startvertices[vi]] = Dict{Profile{T}, T}()
        end
        distancesInfo[startvertices[vi]][startvertices[vj]] =
            ScoreFunc(startvertices[vi], startvertices[vj])
      end
    end
    distancesInfo
  end

  function get{T}(distancesInfo :: Dict{Profile{T}, Dict{Profile{T}, T}},
         k1 :: Profile{T},
         k2 :: Profile{T})
    haskey(distancesInfo, k1) && haskey(distancesInfo[k1], k2) && return distancesInfo[k1][k2]
    haskey(distancesInfo, k2) && haskey(distancesInfo[k2], k1) && return distancesInfo[k2][k1]
    haskey(distancesInfo, k1) && error("invalid k2 key")
    haskey(distancesInfo, k2) && error("invalid k1 key")
    error("invalid keys in distance hash")
  end
end
