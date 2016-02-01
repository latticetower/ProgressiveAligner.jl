module DataWriter
  export writeSequences

  import DataReader.FastaRecord

  function writeSequences(output_file_name :: AbstractString, sequences :: Vector{FastaRecord})

    output_file = open(output_file_name, "w")
    for s in sequences
      println(output_file, ">", s.description)
      println(output_file, s.sequence)
    end
    close(output_file)
  end

end
