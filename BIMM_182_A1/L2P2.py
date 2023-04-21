#use Proch::N50 qw(getStats getN50);
my filepath =
 
# Get N50 only: getN50(file) will return an integer
print "N50 only:\t", getN50($filepath), "\n";
 
# Full stats
my $seq_stats = getStats($filepath);
print Data::Dumper->Dump( [ $seq_stats ], [ qw(*FASTA_stats) ] );
# Will print:
# %FASTA_stats = (
#               'N50' => 65,
#               'N75' => 50,
#               'N90' => 4,
#               'min' => 4,
#               'max' => 65,
#               'dirname' => 'data',
#               'auN' => 45.02112,
#               'size' => 130,
#               'seqs' => 6,
#               'filename' => 'test.fa',
#               'status' => 1
#             );
 
# Get also a JSON object
my $seq_stats_with_JSON = getStats($filepath, 'JSON');
print $seq_stats_with_JSON->{json}, "\n";
# Will print:
# {
#    "status" : 1,
#    "seqs" : 6,
#    <...>
#    "filename" : "small_test.fa",
#    "N50" : 65,
# }
# Directly ask for the JSON object only:
my $json = jsonStats($filepath);
print $json;