module DataReader
  export readSequences, FastaRecord, readMatrix, ScoreMatrix

  type FastaRecord
    description :: ASCIIString
    sequence :: ASCIIString
  end

  type ScoreMatrix
    keys :: Dict{Char, Int}
    hsh :: Array{Array{Float64, 1}, 1}
  end
  ScoreMatrix() = ScoreMatrix(Dict{Char, Int}(), [])


  function readSequences(input_file_name :: String)
    records = FastaRecord[]
    temp_desc = ""
    temp_str = ""
    input_file = open(input_file_name,"r")
    while !eof(input_file)
      s = rstrip(readline(input_file), ['\r','\n'])
      if s[1] == '>'
        if length(temp_str) > 0
          push!(records, FastaRecord(temp_desc, temp_str))
          temp_str = ""
        end
        temp_desc = s[2:end]
      else
        temp_str = string(temp_str, s)
      end
    end

    if (temp_str != "")
      push!(records, FastaRecord(temp_desc, temp_str))
    end
    records
  end

  function readMatrix(input_file :: String)
    result = readdlm(input_file)
    k = [getindex(result, 1, inner_key)[1] => inner_key for inner_key = 1:24]
    h = [
        [
            getindex(result, outer_key, inner_key + 1) :: Float64
            for inner_key = 1:outer_key-1
        ]
        for outer_key = 2 : 25
    ]
    ScoreMatrix(k, h)
  end

end
