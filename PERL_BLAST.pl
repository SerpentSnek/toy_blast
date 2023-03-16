#PERL_BLAST
open(IN,'perlblastdata.txt');
# $/ = "";
print "Give a k value:\n";
$k = <>;
chomp $k;

print "Give a threshold t value:\n";
$t = <>;
chomp $t;


# subroutine to find 4-mers and hash them
sub hashing {
	# my (%hash_var) = @_;
	$query = $_[0]; $kval = $_[1];
	%kmer = ();                      
	$i = 1;
	while (length($query) >= $kval) { # k is 4
		$query =~ m/(.{$kval})/; 
		# print "$1, $i \n";
		if (! defined $kmer{$1}) {
		   $kmer{$1} = $i;       
		}
		else { push (@{$kmer{$1}}, $i)} # push value into an ARRAY associated with the key
		$i++;
		$query = substr($query, 1, length($query) -1);
	}
	return %kmer;
}


# =begin
while($line = <IN>) {
	chomp $line;
    if ($line ne '') {
        $queryString = $line; # finding the query string, which should be the last nonempty line 
    }
} 
# print "$queryString\n";
# finding all 4-mers of the query string
%qKmer = hashing($queryString,$k);
# $qsize = keys %queryKmer;
%dbKmer = ();
%stringhash=();

close(IN);
@querystr =  split(//, $queryString); 
$nq = @querystr;
open(IN,'perlblastdata.txt');

# reading and hashing database strings from perlblastdata.txt
$s=1; # keep track of hashing index of stringhash (the hash table to prevent repeat reports)

while ($line = <IN>) {
	
	chomp $line;
	if($line =~ m/^\w{60}$/)  { # the db string is exaclty 60 chars long
		# hash the string
		# print "line: $line\n";
		%dbKmer = hashing($line,$k); # this is a hash table of S's 4-mers
		# push(@hasharray, %dbKmer);
		foreach my $qkey (keys(%qKmer)) {
			
			if(exists $dbKmer{$qkey}) {
				# determine the range for the substring				
				# convert the query string and S into character array
				
				@dbstr =  split(//, $line); 
				$mdb = @dbstr;
				# loop thru all places that have said kmer, iterate thru the hash
				foreach $dbkey (sort keys(%dbKmer)) {
					$L=$k; 
					# print "dbkey\t$dbkey\n";
					$dbstart=@{$dbKmer{$dbkey}}; # find where to start comparing in S
					$startsubstr=$qKmer{$qkey};

					for($i=$qKmer{$qkey}-1; $i < $nq && $dbstart < $mdb; $i++) {
						
						if($querystr[$i] eq $dbstr[$dbstart]) {
							# print "qcharR: $querystr[$i]\tdbcharR: $dbstr[$dbstart]\n";
							$dbstart++;
							$L++;
						} # end if 
						
						else {last;}
					} # end innermost for($i=$qkey; $i < $nq; $i++)
					
					# extend the search to the LEFT starting from qkey

					$dbstart=$dbKmer{$qkey}-1; # reset dbstart
					for($j=$qKmer{$qkey}-1; $j > 0 && $dbstart > 0; $j--) {	
						
						if($querystr[$j] eq $dbstr[$dbstart]) {
							# print "qcharL: $querystr[$j]\tdbcharL: $dbstr[$dbstart]\n";
							$dbstart--;
							$L++;
							$startsubstr--;
						} # end if 
						
						else {last;}
					} # end innermost for($i=$qkey; $i < $nq; $i++)
			
					if($L > $t) {
						# create a stringhash by hashing the substring
						# $L=$L-$k-1;
						$hspsubstr = substr($queryString, $startsubstr, $L-1);
						if(!defined $stringhash{$hspsubstr}) {
							print "L: $L\tsubstr:\t$hspsubstr\n"; 
							$stringhash{$hspsubstr} = $s;
							$s++;
						}
						
					}
				} # end for each key index
				
			} # endif exists $dbKmer{$qkey}		
			
		} # end foreach my $qkey (keys(%qKmer)
		
	} # end outermost if $line =~ m/^\w{60}$/
	

} # end while


# foreach $kmerkey (keys(%kmer)) {
	# print "The first occurrence of string $kmerkey is in position 
	# $kmer{$kmerkey}\n";
# }
close(IN);
# =cut
#end