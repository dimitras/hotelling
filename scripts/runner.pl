# USAGE: perl runner.pl ../data/spreadsheet_postnatal.csv 4 4 6 210

#use Modern::Perl;
use warnings;
use strict;
use File::Temp;

my ($ifile, $num_reps_cond1, $num_reps_cond2, $num_timepoints, $num_permutations) = @ARGV;

open(my $FH, '<', $ifile) or die "Cannot open file $ifile. $!\n";

my @results = ();
<$FH>;
while (my $l = <$FH>) {
	chomp $l;
	my @a = split(/,/, $l);
	my $id = shift(@a);
	my $name = pop(@a);

	my @colnames = map {"T".$_} 1..$num_timepoints;
	my @rownames = ((map {"C1R".$_} 1..$num_reps_cond1), (map {"C2R".$_} 1..$num_reps_cond2));

	my $x = columns_2d(\@a, $num_reps_cond1, $num_reps_cond2, $num_timepoints);

	if($x){
		my $tmp = File::Temp->new( DIR => "foo/", SUFFIX => '.txt' );
		print $tmp join("\t", "ID", @colnames)."\n";
		for (my $i=0; $i<@$x; $i++) {
			print $tmp join("\t", $rownames[$i], @{$x->[$i]})."\n";
		}
		print $tmp "";

		my $command = "perl hotelling.pl " . $tmp->filename . " $num_reps_cond1 $num_reps_cond2 $num_permutations";
		my $result = `$command`;
		my @split_res = split(/\n/, $result);
		$split_res[0] =~ s/Unpermuted Hotelling Statistic = //i;
		$split_res[-1] =~ s/num_less = //i;
		unshift @split_res, $id;
		push(@results, \@split_res);
	}
}
close $FH;

#print join(",", "id", "unpermuted_statistic", (map{"perm".$_} 1..$num_permutations), "num_less") . "\n"; # if I need the permutated values as well
print join(",", "id", "unpermuted_statistic", "num_less") . "\n";
foreach my $res (@results) {
	print join(",", @$res) . "\n";
}
#
# for(my $i=0; $i<@{$results[0]}; $i++) {
# 	my @row;
# 	for(my $j=0; $j<@results; $j++) {
# 		push @row, $results[$j][$i];
# 	}
# 	print join(",", @row) . "\n";
# }
#

sub columns_2d {
	my ($cols, $num_reps_cond1, $num_reps_cond2, $num_timepoints) = @_;

	my $start = 0;
	my $stop = ($num_reps_cond1 * $num_timepoints) - 1;
	my @a = @{$cols}[$start..$stop];
	my $first = condition_to_2d(\@a, $num_reps_cond1);

	$start = $num_reps_cond1 * $num_timepoints;
	$stop = $start + $num_reps_cond2 * $num_timepoints - 1;
	my @b = @{$cols}[$start..$stop];
	my $second = condition_to_2d(\@b, $num_reps_cond2);

	if($first && $second) {
		return [@{$first}, @{$second}];
	}
}

sub condition_to_2d {
	my ($a, $reps_cnt) = @_;

	my @x;
	my $zerocounts = 0;
	for (my $k=0; $k<@$a; $k++) {
		my $i = $k % $reps_cnt;
		my $j = int($k / $reps_cnt);

		if($a->[$k] == 0){
			$zerocounts++;
		}

		$x[$i][$j] = $a->[$k];
	}
	if($zerocounts<@$a-@$a*0.25){
		return \@x;
	}
	else{
		return;
	}
}
