# ProgressiveAligner

[![Join the chat at https://gitter.im/latticetower/ProgressiveAligner.jl](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/latticetower/ProgressiveAligner.jl?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![Build Status](https://travis-ci.org/latticetower/ProgressiveAligner.jl.svg?branch=master)](https://travis-ci.org/latticetower/ProgressiveAligner.jl)

This package contains progressive alignment tool for protein sequences written in Julia language (http://julialang.org/).

Builds phylogenetic tree with neighbour joining, UPGMA or WPGMA algorithm, then aligns protein sequences by their profiles.

Usage examples currently can be found in `test` folder.

Typical usage pipeline
----------------------

1. call methods from `DataReader` submodule and read data from files.

  One way to get protein sequences - to read them from file:
  ```
  sequences = readSequences(dirname(@__FILE__()) * "/../data/input_test_sequences.faa")
  ```
  After this call, `sequences` is an Array of `FastaRecord` objects. `FastaRecord` object can also be created directly, from description and protein sequence string:
  ```
  fasta_record  = FastaRecord("test string1", "ACGT")
  ```

  Alignment algorithm uses alignment score matrix to score sequences; this matrix can be loaded from alignment matrix file, which can be loaded from [NCBI ftp](ftp://ftp.ncbi.nih.gov/blast/matrices/).
  ```
  matrix = readMatrix(dirname(@__FILE__()) * "/../data/blosum62.txt")
  ```

2. Convert and prepare data.

  Alignment algorithm don't use `FastaRecord` objects directly. It converts protein strings data to Profile objects, then merges these objects with different tree building methods to one `Profile`, which represents multiple alignment and can be converted to `Array` of `FastaRecord` objects with gaps.

  First step is to incorporate `FastaRecord`s or strings to Array of `Profile`s:
  ```
  strToProfiles(strings :: Vector{FastaRecord}) = [Profile{Float64}(record.sequence, record.description) :: Profile{Float64} for record in strings]

  profiles = strToProfiles(sequences)
  ```
  For collection of strings, profile array creation can be done in similar way:
  ```
  start_vertices2 = [Profile{Float64}(str)::Profile{Float64} for str in
      [
        "CAP",
        "CAPT",
        "APT",
        "PPT"
        ]]
  ```

  Second step is to set current scoring matrix. This can be done via call
  ```
  ProfileAligner.setScoringMatrix(score_matrix)
  ```

3. Define score and merge functions:

  ```
  function scoreFunc(p1 :: Profile{Float64}, p2 :: Profile{Float64})
    ProfileAligner.scoreprofiles(p1, p2)
  end

  function mergeFunc(p1 :: Profile{Float64}, p2 :: Profile{Float64})
    ProfileAligner.align(p1, p2)
  end
  ```
  Score function returns best alignment score for given pair of profiles. Merge function returns resulting profile object, which can be build by best-scored alignment. The main difference between these two methods - first one can be be computed faster and can consume less memory.

  Currently there is default implementation for both of these methods in `ProfileAligner` submodule (corresponding methods are shown in previous code example), which uses score matrix, set by `setScoringMatrix` call, to select best-scored profile alignment.

4. Select one of tree-based methods to perform multiple alignment.

  Currently 3 clustering algorithms are implemented - `NeighbourJoining`, `UPGMA`, `WPGMA`. They got similar signature and can be called like that:

  ```
  njResult = NeighbourJoining(profiles, scoreFunc, mergeFunc)
  wpgmaResult = WPGMA(profiles, scoreFunc, mergeFunc)
  upgmaResult = UPGMA(profiles, scoreFunc, mergeFunc)
  ```

  The result of each of these calls is a single Profile object, which represents multiple alignment.

5. Convert alignment result to readable way. Save results to file.

  `ProfileAligner` submodule provides utility method `getstrings`, which converts `Profile` representation back to array of `FastaRecord` objects (probably with gaps).
  ```
  getstrings(result)
  ```

  There is also utility submodule `DataWriter`, which can be used to save resulting multiple sequence alignment to .fasta-like file.
  ```
  writeSequences(output_file, getstrings(result))
  ```
  First parameter should contain file name to save these records.

  There is no differences in file format from typical .fasta, except one - in resulting sequence alignment, each string can contain gaps, represented by '-' symbol. Descriptions are kept to make it possible to find where aligned string came from (and we know that clustering algorithms mix input nodes and can change their order).
