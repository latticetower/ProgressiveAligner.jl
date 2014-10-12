module ProfileAligner

  export AlignmentMatrix, Profile, score, align, getstrings, scoreprofiles, measurequality

  debugprint(str :: Any) = 0
  #debugprintln(str :: Any) = println(str)

  immutable Profile{T}
    rawdata :: Array{Char, 2}
    data :: Array{Dict{Char, T}, 1}
    stringsize :: Int64
    numberofstrings :: Int64

    Profile(raw :: Array{Char, 2}) = getprofile(raw)
    Profile(str :: String) = Profile{T}( reshape([ letter for letter in str ], length(str), 1) )

    function getprofile(raw::Array{Char, 2})
      rawdata = copy(raw)
      data = Array(Dict{Char, T}, size(raw, 1))
      for row in 1 : size(raw, 1)
        data[row] = Dict{Char, T}()
        s = zero(T)
        for col in 1:size(raw, 2)
          letter = raw[row, col]
          data[row][letter] = get!(data[row], letter, zero(T)) + 1
          s += one(T)
        end
        if s > 0
          for (key, value) in data[row]
            data[row][key] /= s
          end
        end
      end
      new(rawdata, data, size(raw, 1), size(raw, 2))
    end
  end

  type AlignmentMatrix{T}
    matrix :: Array{T, 2}
    path :: Array{Char, 2}

    AlignmentMatrix(n :: Int64, m :: Int64) = build(n, m)

    function build(n :: Int64, m :: Int64)
      matrix = zero(Array(T, n, m))
      path = Array(Char, n, m)
      new(matrix, path)
    end
  end


input_file = open(ARGS[1], "r")

result = readdlm(input_file)
hsh = [
    getindex(result, outer_key, 1)[1] => [
      getindex(result, 1, inner_key)[1] =>
        getindex(result, outer_key, inner_key + 1)
          for inner_key = 1:24 ]
    for outer_key = 2:25
]

  #method returns pair of previous coordinates
  function getprev(direction :: Char, i :: Int64, j :: Int64)
    i < 2 && direction != 'R' && error("1st index is less than 2")
    j < 2 && direction != 'D' && error("2nd index is less than 2")
    direction == 'D' && return (i - 1, j)
    direction == 'R' && return (i, j - 1)
    direction == 'M' && return (i - 1, j - 1)
    error("unknown direction")
  end

  import Base.get
  get{T}(profile::Profile{T}, i::Int64, aa::Char) = get(profile.data[i], aa, zero(T))

  aminoacids = "ARNDCQEGHILKMFPSTWYVBZX"

  get_ri{T}(profile::Profile{T}, i::Int64) = length(keys(profile.data[i]))

  get_nij{T}(profile::Profile{T}, i::Int64, aa::Char) = get(profile.data[i], aa, 0.0)

  # 1. define weighting schemes
  function sequenceWeight1{T}(profile::Profile{T}, sequence::String)
    sum([
      1.0/(get_ri(profile, i) * get_nij(profile, i, sequence[i]))
      for i in 1:length(sequence)
    ])
  end

  # 2. define scores
  function score{T}(P::Profile{T}, Q::Profile{T}, i, j)
    sum([ get(P, i, aa) * get(Q, j, aa) for aa in aminoacids ] )
  end

  function score2{T}(P::Profile{T}, Q::Profile{T}, i, j)
    sum([
          get(P, i, aa1) * get(Q, j, aa2) * hsh[aa1][aa2]
            for aa1 in aminoacids, aa2 in aminoacids
        ] )
  end

  GAP_COST = -1

  function prepareAlignmentMatrix{T}(
                P :: Profile{T},
                Q :: Profile{T},
                matrixdata :: AlignmentMatrix{T}
                )
    n = P.stringsize
    m = Q.stringsize
    matrixdata.matrix[1, 1] = zero(T)
    matrixdata.path[1, 1] = 'U'
    for i = 1:n
      matrixdata.matrix[i + 1, 1] = matrixdata.matrix[i, 1] + GAP_COST
      matrixdata.path[i + 1, 1] = 'D'
    end
    for j = 1:m
      matrixdata.matrix[1, j + 1] = matrixdata.matrix[1, j] + GAP_COST
      matrixdata.path[1, j + 1] = 'R'
    end
  end

  function buildAlignmentMatrix{T}(
                  P::Profile{T},
                  Q :: Profile{T},
                  matrixdata :: AlignmentMatrix{T},
                  scoreFunc
                  )
    n = P.stringsize
    m = Q.stringsize
    for i = 1:n
      for j = 1:m
        # 1. Match
        matrixdata.matrix[i + 1, j + 1] = matrixdata.matrix[i, j] + scoreFunc(P, Q, i, j)
        matrixdata.path[i + 1, j + 1] = 'M'
        # 2. Down
        temp_value = matrixdata.matrix[i, j + 1] + GAP_COST
        if (matrixdata.matrix[i + 1, j + 1] < temp_value)
          matrixdata.matrix[i + 1, j + 1] = temp_value
          matrixdata.path[i + 1, j + 1] = 'D'
        end
        # 3. Right
        temp_value = matrixdata.matrix[i + 1, j] + GAP_COST
        if (matrixdata.matrix[i + 1, j + 1] < temp_value)
          matrixdata.matrix[i + 1, j + 1] = temp_value
          matrixdata.path[i + 1, j + 1] = 'R'
        end

      end
    end
  end

  #returns score from precomputed matrix and last column/row coordinates
  function getscore{T}(P::Profile{T}, Q :: Profile{T}, matrixdata :: AlignmentMatrix{T})
    n = P.stringsize
    m = Q.stringsize
    return (n + 1, m + 1, matrixdata.matrix[n + 1, m + 1])
  end

  function backtracing{T}(P :: Profile{T}, Q :: Profile{T},
                          matrixdata :: AlignmentMatrix{T},
                          lastrow :: Int64, lastcol :: Int64)
    debugprint(matrixdata.path)
    path = (Char, Int64, Int64)[]
    (direction, row, col) = (matrixdata.path[lastrow, lastcol], lastrow, lastcol)
    while (row > 1 || col > 1)
      push!(path, (direction, row - 1, col - 1))
      debugprint(matrixdata.path[row, col])
      (row, col) = getprev(direction, row, col)
      direction = matrixdata.path[row, col]
    end
    path
  end

  construct{T}(a1::Array{T}, a2::Array{T}) = hcat(a1, a2)

  function mixprofilecolumn{T}(P :: Profile{T}, Q :: Profile{T},
                                 direction :: Char, i :: Int64, j :: Int64)
    direction == 'M' && return construct(P.rawdata[i, 1:end], Q.rawdata[j, 1:end])
    direction == 'R' && return construct(reshape(['-' for k in 1 : P.numberofstrings], 1, P.numberofstrings),
              Q.rawdata[j, 1:end])
    direction == 'D' && return construct(
              P.rawdata[i, 1:end],
              reshape(['-' for k in 1 : Q.numberofstrings], 1, Q.numberofstrings)
              )
    error("unknown direction in mix profile column")
  end

  function mixprofiles{T}(P :: Profile{T}, Q :: Profile{T}, indices::Vector{(Char, Int64, Int64)})
    tempmatrix = [
      mixprofilecolumn(P, Q, indices[index][1], indices[index][2], indices[index][3])
      for index in 1 : length(indices)
    ]
    debugprint(tempmatrix)
    Profile{T}([
      tempmatrix[i][j]
      for i in length(indices):-1:1, j in length(tempmatrix[1]):-1:1
    ])
  end

  #
  # main alignment function for profiles.
  # returns consensus profile + alignment score
  #
  function align{T}(P :: Profile{T}, Q :: Profile{T})
    debugprint(P)
    debugprint(Q)
    n = P.stringsize
    m = Q.stringsize
    matrix = AlignmentMatrix{T}(n + 1, m + 1)
    prepareAlignmentMatrix(P, Q, matrix)
    buildAlignmentMatrix(P, Q, matrix, score)
    (lastrow, lastcol, sc) = getscore(P, Q, matrix)
    debugprint(sc)

    pairsofindices = backtracing(P, Q, matrix, lastrow, lastcol) # this should return new Profile{T} with resulting alignment and score
    debugprint(pairsofindices)
    res = mixprofiles(P, Q, pairsofindices)
    return res
  end

  function scoreprofiles{T}(P :: Profile{T}, Q :: Profile{T})
    debugprint(P)
    debugprint(Q)
    n = P.stringsize
    m = Q.stringsize
    matrix = AlignmentMatrix{T}(n + 1, m + 1)
    prepareAlignmentMatrix(P, Q, matrix)
    buildAlignmentMatrix(P, Q, matrix, score)
    (lastrow, lastcol, sc) = getscore(P, Q, matrix)
    debugprint(sc)
    return sc
  end

  #
  function measurequality{T}(P :: Profile{T})
    sum([
      sum([- get(P.data[i], aa, zero(T))*log2(get(P.data[i], aa, one(T)))
        for aa in setdiff(keys(P.data[i]), '-')])
      for i in 1 : P.stringsize])
  end

  #
  # method returns string array for given profile
  #
  getstrings{T}(P :: Profile{T}) =  [ ascii(P.rawdata[1:end, i]) for i in 1 : P.numberofstrings]

end
