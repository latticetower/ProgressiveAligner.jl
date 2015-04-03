importall ProfileAligner
importall Clustering
importall DataReader

start_vertices2 = [Profile{Float64}(str)::Profile{Float64} for str in
    [
      "CAP",
      "CAPT",
      "APT",
      "PPT"
      ]]

start_vertices = [Profile{Float64}(str)::Profile{Float64} for str in
    [
      "PEEKSAVTALWGKVNVDEVGG",
      "GEEKAAVLALWDKVNEEEVGG",
      "PADKTNVKAAWGKVGAHAGEYGA",
      "AADKTNVKAAWSKVGGHAGEYGA",
      "EHEWQLVLHVWAKVEADVAGHGQ"
      ]]

function scoreFunc(p1 :: Profile{Float64}, p2 :: Profile{Float64})
  scoreprofiles(p1, p2)
end

function mergeFunc(p1 :: Profile{Float64}, p2 :: Profile{Float64})
  align(p1, p2)
end

#TODO: add test to show difference between different tree construction methods
resulting_profile = UPGMA(start_vertices, scoreFunc, mergeFunc)
profile_strings1 = [
  FastaRecord("","GEEKAAVLALWDKV---NEEEVGG"),
  FastaRecord("","PADKTNVKAAWGKVGA-HAGEYGA"),
  FastaRecord("","PEEKSAVTALWGKVNVDEVG--G-"),
  FastaRecord("","AADKTNVKAAWSKV-GGHAGEYGA"),
  FastaRecord("","EHEWQLVLHVWAKVEADVAG-HGQ")
]

profile_strings_result1 = getstrings(resulting_profile)
for i in 1:5
  @test profile_strings_result1[i] == profile_strings1[i]
end

println(measurequality(resulting_profile))

resulting_profile2 = WPGMA(start_vertices, scoreFunc, mergeFunc)
profile_strings2 = [
  FastaRecord("","GEEKAAVLALWDKV---NEEEVGG"),
  FastaRecord("","PADKTNVKAAWGKVGA-HAGEYGA"),
  FastaRecord("","PEEKSAVTALWGKVNVDEVG--G-"),
  FastaRecord("","AADKTNVKAAWSKV-GGHAGEYGA"),
  FastaRecord("","EHEWQLVLHVWAKVEADVAG-HGQ")
]
profile_strings_result2 = getstrings(resulting_profile2)
for i in 1:5
  @test profile_strings_result2[i] == profile_strings2[i]
end
println(measurequality(resulting_profile2))

resulting_profile3 = NeighbourJoining(start_vertices, scoreFunc, mergeFunc)
profile_strings_result3 = getstrings(resulting_profile2)
profile_strings3 = [
  FastaRecord("","GEEKAAVLALWDKV---NEEEVGG"),
  FastaRecord("","PADKTNVKAAWGKVGA-HAGEYGA"),
  FastaRecord("","PEEKSAVTALWGKVNVDEVG--G-"),
  FastaRecord("","AADKTNVKAAWSKV-GGHAGEYGA"),
  FastaRecord("","EHEWQLVLHVWAKVEADVAG-HGQ")
]
for i in 1:5
  @test profile_strings_result3[i] == profile_strings3[i]
end
println(measurequality(resulting_profile3))
