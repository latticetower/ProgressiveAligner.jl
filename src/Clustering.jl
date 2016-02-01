module Clustering
  export UPGMA, WPGMA, NeighbourJoining

  import ProfileAligner.Profile

  #typealias DistanceHash{T} Dict{Profile{T}, Dict{Profile{T}, T}}

  type DistanceMatrix{T}
    keys :: Dict{Profile{T}, Int64}
    hsh :: Array{T, 2}
    K :: Vector{T}
    size :: Int
    DistanceMatrix(n :: Int) = new(Dict{Profile{T}, Int64}(), zero(Array(T, 2*n, 2*n)), zero(Array(T, 2*n)), 0)
  end

  function haskey{T}(matrix :: DistanceMatrix{T}, p :: Profile{T})
    Base.haskey(matrix.keys, p)
  end

  function setValue{T}(
      distancesInfo :: DistanceMatrix{T},
      vi :: Profile{T},
      vj :: Profile{T},
      value :: T
      )
    if ! Base.haskey(distancesInfo.keys, vi)
      distancesInfo.keys[vi] = distancesInfo.size + 1
      distancesInfo.size += 1
    end
    if ! Base.haskey(distancesInfo.keys, vj)
      distancesInfo.keys[vj] = distancesInfo.size + 1
      distancesInfo.size += 1
    end
    distancesInfo.hsh[distancesInfo.keys[vi], distancesInfo.keys[vj]] = value
    distancesInfo.hsh[distancesInfo.keys[vj], distancesInfo.keys[vi]] = value
  end



  function getUPGMAWeight{T}(distancesInfo :: DistanceMatrix{T},
                          k :: Profile{T}, f :: Profile{T}, g :: Profile{T})
    (f.numberofstrings*get(distancesInfo, f, k) + g.numberofstrings*get(distancesInfo, g, k))/(g.numberofstrings + f.numberofstrings)
  end

  getCurrentPGMADistance{T}(
          distancesInfo :: DistanceMatrix{T},
          vertices :: Vector{Profile{T}},
          vertsize :: Int64,
          i :: Int,
          j :: Int
          ) = get(distancesInfo, vertices[i], vertices[j])

  function getCurrentNJDistance{T}(
          distancesInfo :: DistanceMatrix{T},
          vertices :: Vector{Profile{T}},
          vertsize :: Int64,
          i :: Int,
          j :: Int
          )
    (vertsize - 2)*get(distancesInfo, vertices[i], vertices[j]) - distancesInfo.K[distancesInfo.keys[vertices[i]]] - distancesInfo.K[distancesInfo.keys[vertices[j]]]
  end

  getWPGMAWeight{T}(
          distancesInfo :: DistanceMatrix{T},
          k :: Profile{T},
          f :: Profile{T},
          g :: Profile{T}
          ) = (get(distancesInfo, f, k) + get(distancesInfo, g, k)) * 0.5


  UPGMA{T}(startvertices :: Vector{Profile{T}}, scoreFunc, mergeFunc) =
          buildTree(startvertices,
              scoreFunc, mergeFunc, getUPGMAWeight)

  WPGMA{T}(startvertices :: Vector{Profile{T}}, scoreFunc, mergeFunc) = buildTree(
          startvertices, scoreFunc, mergeFunc, getWPGMAWeight)

  NeighbourJoining{T}(startvertices :: Vector{Profile{T}}, scoreFunc, mergeFunc) =
    buildTree(startvertices, scoreFunc, mergeFunc,
        getNJWeight, getNJParentDist,
        getCurrentNJDistance)


  function getNJWeight{T}(distancesInfo :: DistanceMatrix{T},
          k :: Profile{T},
          f :: Profile{T},
          g :: Profile{T})
    (
      get(distancesInfo, k, f) +
      get(distancesInfo, k, g) -
      get(distancesInfo, f, g)
    ) * 0.5
  end


  function getPGMAParentDist{T}(
          distancesInfo :: DistanceMatrix{T},
          vertices :: Vector{Profile{T}},
          vertsize :: Int64,
          i :: Int,
          j :: Int)
    get(distancesInfo, vertices[i], vertices[j]) * 0.5

  end

  function getNJParentDist{T}(
          distancesInfo :: DistanceMatrix{T},
          vertices :: Vector{Profile{T}},
          vertsize :: Int64,
          i :: Int,
          j :: Int)
    get(distancesInfo, vertices[i], vertices[j]) / 2.0
      + (distancesInfo.K[distancesInfo.keys[vertices[i]]]
      - distancesInfo.K[distancesInfo.keys[vertices[j]]])/(2.0*(vertsize - 2))
  end

  function buildTree{T}(startvertices :: Vector{Profile{T}},
          scoreFunc, mergeFunc, weightFunc,
          parentWeightFunc = getPGMAParentDist,
          currentDistFunc = getCurrentPGMADistance
          )
    distancesInfo = initDistanceMatrix(startvertices, scoreFunc)
    while(length(startvertices) > 1)
      simplifyTree!(
          distancesInfo,
          startvertices,
          mergeFunc,
          parentWeightFunc,
          weightFunc,
          currentDistFunc
          )
    end
    startvertices[1]
  end

  function simplifyTree!{T}(
          distancesInfo :: DistanceMatrix{T},
          vertices :: Vector{Profile{T}},
          mergeFunc,
          parentWeightFunc,
          weightFunc,
          currentDistFunc
          )
    n = length(vertices)
    #println("in simplify")
    #1. calculate Q

    (mini, minj) = (1, 2)
    Qmin = currentDistFunc(distancesInfo, vertices, n, 1, 2)

    for i = 1 : n - 1
      for j = i + 1 : n
        #println(vertices[i], vertices[j])
        Qcurrent = currentDistFunc(distancesInfo, vertices, n, i, j)
        if Qmin > Qcurrent
          (mini, minj) = (i, j)
          Qmin = Qcurrent
        end
      end
    end
    #
    newT = mergeFunc(vertices[mini], vertices[minj])
    computeDistances!(distancesInfo, vertices, n, mini, minj, newT, parentWeightFunc, weightFunc)

    #
    deleteat!(vertices, [mini, minj])
    push!(vertices, newT)
  end

  function computeDistances!{T}(
          distancesInfo :: DistanceMatrix{T},
          vertices :: Vector{Profile{T}},
          vertsize :: Int64,
          mini :: Int,
          minj :: Int,
          newT :: Profile{T},
          parentWeightFunc,
          weightFunc
          )
    #compute new distances from new to previous 2
    #distancesInfo[newT] = Dict{Profile{T}, T}()

    f = vertices[mini]
    g = vertices[minj]
    t2 = one(T) + one(T)
    setValue(distancesInfo, newT, f, parentWeightFunc(distancesInfo, vertices, vertsize, mini, minj))
    setValue(distancesInfo, newT, g, get(distancesInfo, f, g) - get(distancesInfo, newT, f))
    setValue(distancesInfo, newT, newT, zero(T)) #assume that profile aligns to profile with zero distance
    #compute to all others
    for i in 1 : length(vertices)
      i == mini && continue
      i == minj && continue
      k = vertices[i]
      setValue(distancesInfo, newT, k, weightFunc(distancesInfo, k, f, g))
      distancesInfo.K[distancesInfo.keys[k]] += (-get(distancesInfo, k, g) - get(distancesInfo, k, f) + get(distancesInfo, k, newT))
      distancesInfo.K[distancesInfo.keys[newT]] += get(distancesInfo, k, newT)
    end
  end

  function initDistanceMatrix{T}(startvertices :: Vector{Profile{T}}, scoreFunc)
    n = length(startvertices)
    distancesInfo = DistanceMatrix{T}(n)
    for i in 1 : n
      vi = startvertices[i]
      for j in i : n
        vj = startvertices[j]
        setValue(distancesInfo, vi,vj,
                scoreFunc(vi, vj))
        distancesInfo.K[distancesInfo.keys[vi]] += get(distancesInfo, vi, vj)
        if i != j
          distancesInfo.K[distancesInfo.keys[vj]] += get(distancesInfo, vi, vj)
        end
      end
    end
    distancesInfo
  end

  function get2{T}(distancesInfo :: DistanceMatrix{T}, k1 :: Profile{T}, k2 :: Profile{T})
    distancesInfo.hsh[distancesInfo.keys[k1], distancesInfo.keys[k2]]
  end

  function getKValue{T}(distancesInfo :: DistanceMatrix{T}, k1 :: Profile{T})
    distancesInfo.K[distancesInfo.keys[k1]]
  end


  function get{T}(distancesInfo :: DistanceMatrix{T}, k1 :: Profile{T}, k2 :: Profile{T})
    Base.haskey(distancesInfo.keys, k1) && Base.haskey(distancesInfo.keys, k2) && return get2(distancesInfo, k1, k2)
    Base.haskey(distancesInfo.keys, k1) && error("invalid k2 key in distance hash")
    Base.haskey(distancesInfo.keys, k2) && error("invalid k1 key in distance hash")
    error("invalid keys in distance hash")
  end
end
