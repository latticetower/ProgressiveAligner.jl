module DataReader
  export readSequences, FastaRecord, readMatrix, ScoreMatrix

  import Base.==

  immutable FastaRecord
    description :: ASCIIString
    sequence :: ASCIIString

  end

  ==(a::FastaRecord, b::FastaRecord) = a.description == b.description && a.sequence == b.sequence

  immutable ScoreMatrix
    keys :: Dict{Char, Int}
    hsh  :: Array{Array{Float64, 1}, 1}
  end
  ScoreMatrix() = ScoreMatrix(Dict{Char, Int}(), [])


  function readSequences(input_file_name :: AbstractString)
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

  #helper function to avoid problems with readdlm from standard lib (empty lines skipping, \r\n, etc.)
  function my_readdlm(input_file :: AbstractString)
    removeEol = str-> strip(str, ['\r', '\n'])
    notCommentOrEmptyLine = x -> !startswith(x, "#") && !isempty(x)
    getNumbersAndChars = n -> ismatch(r"[^A-Z*]", n) ? float(n) : n[1]

    linesProcessor = lines -> map( x-> map(getNumbersAndChars, split(x)),
                              filter(notCommentOrEmptyLine, map(removeEol, lines)))

    f = open(input_file)
    lines = linesProcessor(readlines(f))
    close(f)
    lines
  end

  function readMatrix(input_file :: AbstractString)
    result = my_readdlm(input_file)
    #result = readdlm(input_file)
    # thus dirty fix needed to avoid error on empty lines in input file. in julia 0.3 readdlm doesn't skip them
    k = [ getindex(getindex(result, 1), inner_key)[1] => inner_key for inner_key = 1:24]
    h = [
        [
            getindex(getindex(result, outer_key), inner_key + 1) :: Float64
            for inner_key = 1 : outer_key - 1
        ]
        for outer_key = 2 : 25
    ]
    ScoreMatrix(k, h)
  end

end
