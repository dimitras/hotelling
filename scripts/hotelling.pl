use Math::Matrix;
$num_reps_cond1 = $ARGV[1];
$num_reps_cond2 = $ARGV[2];
$num_perms = $ARGV[3];
open(INFILE, $ARGV[0]);
$line = <INFILE>;
for($row=0; $row<$num_reps_cond1; $row++) {
    $line = <INFILE>;
    chomp($line);
    @a = split(/\t/,$line);
    $numcols = @a;      # this count includes the id column
    for($i=1; $i<$numcols; $i++) {
	$x[$row][$i-1] = $a[$i];
    }
}
for($row=0; $row<$num_reps_cond2; $row++) {
    $line = <INFILE>;
    chomp($line);
    @a = split(/\t/,$line);
    $numcols = @a;      # this count includes the id column
    for($i=1; $i<$numcols; $i++) {
	$y[$row][$i-1] = $a[$i];
    }
}

#my @perm_array = ();

$unpermuted_hotelling_stat = Hotelling(\@x, \@y);
#push @perm_array, $unpermuted_hotelling_stat;
print "Unpermuted Hotelling Statistic = $unpermuted_hotelling_stat\n";

$num_less = 0;
for($perm=0; $perm<$num_perms; $perm++) {
    for($i=0; $i<1000; $i++) {
	$first = int(rand($num_reps_cond1));
	$second = int(rand($num_reps_cond2));
	$flip = int(rand(2));
	if($flip == 1) {
	    for($j=0; $j<$numcols; $j++) {
		$temp[$j] = $x[$first][$j];
	    }
	    for($j=0; $j<$numcols; $j++) {
		$x[$first][$j] = $y[$second][$j];
	    }
	    for($j=0; $j<$numcols; $j++) {
		$y[$second][$j] = $temp[$j];
	    }
	}
    }
    $permuted_hotelling_stat = Hotelling(\@x, \@y);
    if($permuted_hotelling_stat < $unpermuted_hotelling_stat) {
	$num_less++;
    }
	#push @perm_array, $permuted_hotelling_stat;
   # print "$permuted_hotelling_stat\n";
}

#push @perm_array, $num_less;
print "num_less = $num_less\n";

#for(my $row=0; $row<scalar(@perm_array); $row++){
#    print $perm_array[$row].", ";
#}
#print "";

sub Hotelling {
    ($M1_ref, $M2_ref) = @_;
    @{$M[0]} = @{$M1_ref};
    @{$M[1]} = @{$M2_ref};

    $num_vectors[0] = @{$M[0]};
    $num_vectors[1] = @{$M[1]};

    if($num_vectors[0] == 0 || $num_vectors[1] == 0) {
	print "Error: must be more than zero vectors\n";
	exit();
    }

    $vect_length1 = @{$M[0][0]};
    $vect_length2 = @{$M[1][0]};

    if($vect_length1 != $vect_length2) {
	print "Error: vector lengths not equal\n";
	exit();
    }
    if($vect_length1 == 0 || $vect_length2 == 0) {
	print "Error: vector lengths cannot be zero\n";
	exit();
    }
    $vect_length = $vect_length1;
    for($i=0; $i<$vect_length; $i++) {
	$zero_array[$i] = 0;
    }
    for($i=0; $i<$vect_length; $i++) {
	for($j=0; $j<$vect_length; $j++) {
	    $zero_2darray[$i][$j] = 0;
	    if($i == $j) {
		$identity_array[$i][$j] = 1;
	    }
	    else {
		$identity_array[$i][$j] = 0;
	    }
	}
    }

    for($group=0; $group<2; $group++) {
	for($i=0; $i<$num_vectors[$group]; $i++) {
	    $vectors[$group][$i] = new Math::Matrix($M[$group][$i]);
#	    print "v$i: ";
#	    $vectors[$group][$i]->print;
	}
    }
    for($group=0; $group<2; $group++) {
	$meanvector[$group] = new Math::Matrix(\@zero_array);
	$factor[$group] = new Math::Matrix(@zero_2darray);
	for($i=0; $i<$num_vectors[$group]; $i++) {
	    $meanvector[$group] = $meanvector[$group]->add($vectors[$group][$i]);
	}
	$f = 1/$num_vectors[$group];
	$meanvector[$group] = $meanvector[$group]->multiply_scalar($f);
#	print "mean group $group:\n";
#	$meanvector[$group]->print;
	for($i=0; $i<$num_vectors[$group]; $i++) {
	    $D1 = $vectors[$group][$i]->add($meanvector[$group]->multiply_scalar(-1));
	    $D2 = $D1->transpose;
	    $MM = $D2->multiply($D1);
	    $factor[$group] = $factor[$group]->add($MM);
	}
    }
    $F = $factor[0]->add($factor[1]);
    $W = $F->multiply_scalar(1/($num_vectors[0] + $num_vectors[1] - 2));
#    $W -> print("W:\n");
    $identity_matrix = new Math::Matrix(@identity_array);
#    $identity_matrix->print("Id\n");
    $C = $W->concat($identity_matrix);
    $W_inv = $C->solve;
#    $W_inv->print("W inverse\n");

#    $I = $W->multiply($W_inv);
#    $I->print("W * W_inv\n");

    $DR = $meanvector[0]->add($meanvector[1]->multiply_scalar(-1));
    $DT = $DR->transpose;

    $x = ($DR->multiply($W_inv))->multiply($DT);
#    $x->print("x\n");
    $coeff = ($num_vectors[0] * $num_vectors[1]) / ($num_vectors[0] + $num_vectors[1]);
    $z = $x->multiply_scalar($coeff);
#    $z->print("Hotelling: \n");
    
    $z =~ s/\s//gs;
    $z = $z + 0;

    return $z;
}

