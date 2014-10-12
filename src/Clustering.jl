module Clustering
  export UPGMA, WPGMA, NeighbourJoining
  import ProfileAligner.Profile

  typealias DistanceHash{T} Dict{Profile{T}, Dict{Profile{T}, T}}

  function getUPGMAWeight{T}(distancesInfo :: DistanceHash{T},
                          k :: Profile{T}, f :: Profile{T}, g :: Profile{T})
    (f.numberofstrings*get(distancesInfo, f, k) + g.numberofstrings*get(distancesInfo, g, k))/(g.numberofstrings + f.numberofstrings)
  end

  getCurrentPGMADistance{T} (
          distancesInfo :: DistanceHash{T},
          vertices :: Vector{Profile{T}},
          K1 :: Vector{T},
          K2 :: Vector{T},
          i :: Int,
          j :: Int
          ) = get(distancesInfo, vertices[i], vertices[j])

  function getCurrentNJDistance{T}(
          distancesInfo :: DistanceHash{T},
          vertices :: Vector{Profile{T}},
          K1 :: Vector{T},
          K2 :: Vector{T},
          i :: Int,
          j :: Int
          )
    (length(vertices) - one(T) - one(T))*get(distancesInfo, vertices[i], vertices[j]) - K1[i] - K2[j]
  end

  getWPGMAWeight{T}(
          distancesInfo :: DistanceHash{T},
          k :: Profile{T},
          f :: Profile{T},
          g :: Profile{T}
          ) = (get(distancesInfo, f, k) + get(distancesInfo, g, k))/(one(T) + one(T))


  UPGMA{T}(startvertices :: Vector{Profile{T}}, scoreFunc, mergeFunc) =
          buildTree(startvertices,
              scoreFunc, mergeFunc, getUPGMAWeight)

  WPGMA{T}(startvertices :: Vector{Profile{T}}, scoreFunc, mergeFunc) = buildTree(
          startvertices, scoreFunc, mergeFunc, getWPGMAWeight)

  NeighbourJoining{T}(startvertices :: Vector{Profile{T}}, scoreFunc, mergeFunc) =
    buildTree(startvertices, scoreFunc, mergeFunc,
        getNJWeight, getNJParentDist,
        precomputeK, getCurrentNJDistance)


  function getNJWeight{T}(distancesInfo :: DistanceHash{T},
          k :: Profile{T},
          f :: Profile{T},
          g :: Profile{T})
    (
      get(distancesInfo, k, f) +
      get(distancesInfo, k, g) -
      get(distancesInfo, f, g)
    ) / (one(T) + one(T))
  end


  function precomputeK{T}(distancesInfo :: DistanceHash{T}, vertices :: Vector{Profile{T}})
    n = length(vertices)
    K1 = zero(Array(T, n))
    K2 = zero(Array(T, n))
    for i = 1 : n
      for j = 1 : n
        K1[i] += get(distancesInfo, vertices[i], vertices[j])
        K2[j] += get(distancesInfo, vertices[i], vertices[j])
      end
    end
    (K1, K2)
  end

  function precomputeEmptyK{T}(distancesInfo :: DistanceHash{T}, vertices ::Vector{Profile{T}})
    (Array(T, 0), Array(T, 0))
  end


  function getPGMAParentDist{T}(
          distancesInfo :: DistanceHash{T},
          vertices :: Vector{Profile{T}},
          i :: Int,
          j :: Int,
          K1 :: Vector{T})
    get(distancesInfo, vertices[i], vertices[j]) / (one(T) + one(T))

  end

  function getNJParentDist{T}(
          distancesInfo :: DistanceHash{T},
          vertices :: Vector{Profile{T}},
          i :: Int,
          j :: Int,
          K1 :: Vector{T})
    t2 = one(T) + one(T)
    get(distancesInfo, vertices[i], vertices[j]) / t2 + (K1[i] - K1[j])/(t2*(length(vertices) - t2))
  end

  function buildTree{T}(startvertices :: Vector{Profile{T}},
          scoreFunc, mergeFunc, weightFunc,
          parentWeightFunc = getPGMAParentDist,
          getK = precomputeEmptyK,
          currentDistFunc = getCurrentPGMADistance
          )
    println("in build")
    distancesInfo = initDistanceMatrix(startvertices, scoreFunc)
    while(length(startvertices) > 1)
      simplifyTree!(
          distancesInfo,
          startvertices,
          mergeFunc,
          parentWeightFunc,
          weightFunc,
          getK,
          currentDistFunc
          )
    end
    startvertices[1]
  end

  function simplifyTree!{T}(
          distancesInfo :: DistanceHash{T},
          vertices :: Vector{Profile{T}},
          mergeFunc,
          parentWeightFunc,
          weightFunc,
          getK,
          currentDistFunc
          )
    n = length(vertices)
    (K1, K2) = getK(distancesInfo, vertices)
    #println("in simplify")
    #1. calculate Q

    (mini, minj) = (1, 2)
    Qmin = currentDistFunc(distancesInfo, vertices, K1, K2, 1, 2)

    for i = 1 : n - 1
      for j = i + 1 : n
        #println(vertices[i], vertices[j])
        Qcurrent = currentDistFunc(distancesInfo, vertices, K1, K2, i, j)
        if Qmin > Qcurrent
          (mini, minj) = (i, j)
          Qmin = Qcurrent
        end
      end
    end
    #
    newT = mergeFunc(vertices[mini], vertices[minj])
    computeDistances!(distancesInfo, vertices, mini, minj, newT, K1, parentWeightFunc, weightFunc)
    #
    deleteat!(vertices, [mini, minj])
    push!(vertices, newT)
  end

  function computeDistances!{T}(
          distancesInfo :: DistanceHash{T},
          vertices :: Vector{Profile{T}},
          mini :: Int,
          minj :: Int,
          newT :: Profile{T},
          K1::Vector{T},
          parentWeightFunc,
          weightFunc
          )
    #compute new distances from new to previous 2
    distancesInfo[newT] = Dict{Profile{T}, T}()
    f = vertices[mini]
    g = vertices[minj]
    t2 = one(T) + one(T)
    distancesInfo[newT][f] = parentWeightFunc(distancesInfo, vertices, mini, minj, K1)
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

  function initDistanceMatrix{T}(startvertices :: Vector{Profile{T}}, scoreFunc)
    distancesInfo = DistanceHash{T}()
    for vi in 1 : length(startvertices)
      for vj in vi : length(startvertices)
        if !haskey(distancesInfo, startvertices[vi])
          distancesInfo[startvertices[vi]] = Dict{Profile{T}, T}()
        end
        distancesInfo[startvertices[vi]][startvertices[vj]] =
            scoreFunc(startvertices[vi], startvertices[vj])
      end
    end
    distancesInfo
  end

  function get{T}(distancesInfo :: DistanceHash{T}, k1 :: Profile{T}, k2 :: Profile{T})
    haskey(distancesInfo, k1) && haskey(distancesInfo[k1], k2) && return distancesInfo[k1][k2]
    haskey(distancesInfo, k2) && haskey(distancesInfo[k2], k1) && return distancesInfo[k2][k1]
    haskey(distancesInfo, k1) && error("invalid k2 key")
    haskey(distancesInfo, k2) && error("invalid k1 key")
    error("invalid keys in distance hash")
  end
end
