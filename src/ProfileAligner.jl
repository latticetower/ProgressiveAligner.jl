module ProfileAligner

  export AlignmentMatrix, Profile,
         score, align, getstrings, scoreprofiles,
         measurequality, setScoringMatrix

  import DataReader.FastaRecord, DataReader.ScoreMatrix

  scoringMatrix = ScoreMatrix()

  setScoringMatrix(sm :: ScoreMatrix) = scoringMatrix = sm

  getScoringMatrixValue(matrix :: ScoreMatrix, i1 ::Int64, i2::Int64) = i1 > i2 ? matrix.hsh[i1][i2] : matrix.hsh[i2][i1]

  function getScoringMatrixValue(aa1 :: Char, aa2 :: Char)
    haskey(scoringMatrix.keys, aa1) && haskey(scoringMatrix.keys, aa2) &&
      return getScoringMatrixValue(scoringMatrix, scoringMatrix.keys[aa1], scoringMatrix.keys[aa2])
    aa1 == aa2 ? 1 : -1
  end

  debugprint(str :: Any) = 0
  #debugprintln(str :: Any) = println(str)

  immutable Profile{T}
    rawdata :: Array{Char, 2}
    data :: Array{Dict{Char, T}, 1}
    stringsize :: Int
    numberofstrings :: Int
    descriptions :: Array{ASCIIString, 1}

    #Profile(raw :: Array{Char, 2}, desc :: Array{ASCIIString, 1}) = getprofile{T}(raw, desc)
    #Profile(raw :: Array{Char, 2}, desc :: Array{String, 1}) = getprofile(raw, desc)
    Profile(str :: ASCIIString, desc :: ASCIIString = "") = Profile{T}( reshape([ letter for letter in str ], length(str), 1), [desc] )
    Profile(record :: FastaRecord) = Profile{T}(record.sequence, record.description)

    Profile(raw :: Array{Char, 2},
            d :: Array{Dict{Char, T}, 1},
            s1 :: Int,
            s2 :: Int,
            desc :: Array{ASCIIString, 1}) = new{T}(raw, d, s1, s2, desc)
    #function getprofile(raw :: Array{Char, 2}, descriptions :: Array{ASCIIString, 1})

    #function call{T}(::Type{SummedArray}, a::Vector{T})
    function  call{T}(::Type{Profile{T}}, raw :: Array{Char, 2}, descriptions :: Array{ASCIIString, 1} = [])
      rawdata = copy(raw)
      size_1 = size(raw, 1)
      size_2 = size(raw, 2)
      data = Array(Dict{Char, T}, size_1)
      for row in 1 : size_1
        data[row] = Dict{Char, T}()
        s = zero(T)
        for col in 1 : size_2
          letter = raw[row, col]
          data[row][letter] = get!(data[row], letter, zero(T)) + one(T)
          s += one(T)
        end
        if s > 0
          for (key, value) in data[row]
            data[row][key] /= s
          end
        end
      end
      Profile{T}(rawdata, data, size_1, size_2, descriptions)
    end
  end

  type AlignmentMatrix{T}
    matrix :: Array{T, 2}
    path :: Array{Char, 2}

    AlignmentMatrix(m :: Array{T, 2}, p :: Array{Char, 2}) = new(m, p)

    function call{T}(::Type{AlignmentMatrix{T}}, n :: Int64, m :: Int64)
      matrix = zero(Array(T, n, m))
      path = Array(Char, n, m)
      AlignmentMatrix{T}(matrix, path)
    end
  end


  #method returns pair of previous coordinates
  function getprev(direction :: Char, i :: Int64, j :: Int64)
    i < 2 && direction != 'R' && error("1st index is less than 2")
    j < 2 && direction != 'D' && error("2nd index is less than 2")
    direction == 'D' && return (i - 1, j)
    direction == 'R' && return (i, j - 1)
    direction == 'M' && return (i - 1, j - 1)
    error("unknown direction")
  end

  get{T}(profile::Profile{T}, i::Int64, aa :: Char) = Base.get(profile.data[i], aa, zero(T))

  aminoacids = "ARNDCQEGHILKMFPSTWYVBZX"

  get_ri{T}(profile::Profile{T}, i::Int64) = length(keys(profile.data[i]))

  get_nij{T}(profile::Profile{T}, i::Int64, aa::Char) = Base.get(profile.data[i], aa, 0.0)

  # 1. define weighting schemes
  function sequenceWeight1{T}(profile::Profile{T}, sequence::ASCIIString)
    sum([
      1.0/(get_ri(profile, i) * get_nij(profile, i, sequence[i]))
      for i in 1:length(sequence)
    ])
  end

  # 2. define scores
  function score2{T}(P::Profile{T}, Q::Profile{T}, i, j)
    sum([ get(P, i, aa) * get(Q, j, aa) for aa in aminoacids ] )
  end

  function score{T}(P::Profile{T}, Q::Profile{T}, i, j)
    sum([
          get(P, i, aa1) * get(Q, j, aa2) * getScoringMatrixValue(aa1, aa2)
            for aa1 in aminoacids, aa2 in aminoacids
        ] )
  end



  GAP_COST = -1.0

  function prepareAlignmentMatrix{T}(
                P :: Profile{T},
                Q :: Profile{T},
                matrixdata :: AlignmentMatrix{T}
                )
    n = P.stringsize
    m = Q.stringsize
    matrixdata.matrix[1, 1] = zero(T)
    matrixdata.path[1, 1] = 'U'
    for i = 1 : n
      matrixdata.matrix[i + 1, 1] = matrixdata.matrix[i, 1] + GAP_COST
      matrixdata.path[i + 1, 1] = 'D'
    end
    for j = 1 : m
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
    for j = 1 : m
      for i = 1 : n
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
    path = Tuple{Char, Int64, Int64}[]
    (direction, row, col) = (matrixdata.path[lastrow, lastcol], lastrow, lastcol)
    while (row > 1 || col > 1)
      push!(path, (direction, row - 1, col - 1))
      debugprint(matrixdata.path[row, col])
      (row, col) = getprev(direction, row, col)
      direction = matrixdata.path[row, col]
    end
    path
  end

  construct{T}(a1 :: Array{T}, a2 :: Array{T}) = hcat(a1, a2)

  function mixprofilecolumn{T}(P :: Profile{T}, Q :: Profile{T},
                                 direction :: Char, i :: Int64, j :: Int64)
    direction == 'M' && return construct(P.rawdata[i, 1:end], Q.rawdata[j, 1:end])
    direction == 'R' && return construct(reshape(['-' for k in 1 : P.numberofstrings], 1, P.numberofstrings),
              Q.rawdata[j, 1 : end])
    direction == 'D' && return construct(
              P.rawdata[i, 1 : end],
              reshape(['-' for k in 1 : Q.numberofstrings], 1, Q.numberofstrings)
              )
    error("unknown direction in mix profile column")
  end

  function mixprofiles{T}(P :: Profile{T}, Q :: Profile{T}, indices :: Vector{Tuple{Char, Int64, Int64}})

    newProfileSize = length(indices)
    tempMatrix = Array(Char, newProfileSize, P.numberofstrings + Q.numberofstrings)
    for index in  1 : newProfileSize
      newProfileColumn = mixprofilecolumn(P, Q,
          indices[index][1], indices[index][2], indices[index][3])
      for j in 1 : length(newProfileColumn)
        tempMatrix[newProfileSize - index + 1, j] = newProfileColumn[j]
      end
    end

    Profile{T}(tempMatrix, append!(P.descriptions, Q.descriptions))
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
      sum([- Base.get(P.data[i], aa, zero(T))*log2(Base.get(P.data[i], aa, one(T)))
        for aa in setdiff(keys(P.data[i]), '-')])
      for i in 1 : P.stringsize])
  end

  #
  # method should return FastaRecord array for given profile
  #
  getstrings{T}(P :: Profile{T}) =  [ FastaRecord(P.descriptions[i], ascii(P.rawdata[1:end, i])) for i in 1 : P.numberofstrings]

end
