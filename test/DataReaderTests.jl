importall DataReader

fasta_record  = FastaRecord("test string1", "ACGT")
fasta_record2 = FastaRecord("test string1", "ACGT")

println("testing FastaRecord...")
@test fasta_record.description == "test string1"
@test fasta_record.sequence == "ACGT"

@test fasta_record == fasta_record2

println("testing FastaRecords loading from file...")

sequences = readSequences(dirname(@__FILE__()) * "/../data/input_test_sequences.faa")
@test length(sequences) > 0
@test length(sequences) == 10
@test sequences[1].description == "1"
@test sequences[1].sequence == "LWYVRMHMYR"
@test sequences[10].description == "10"
@test sequences[10].sequence == "QNHHVGNGPA"

println("testing ScoreMatrix loading from sample BLOSUM62 matrix file from ncbi ftp...")

matrix = readMatrix(dirname(@__FILE__()) * "/../data/blosum62.txt")
@test length(matrix.keys) == length(matrix.hsh) == 24

for s = 1 : length(matrix.hsh)
  @test length(matrix.hsh[s]) == s
end
