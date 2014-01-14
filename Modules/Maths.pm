package Maths;

use base qw(Exporter);
our @EXPORT = qw(determineInfluentialPoints calculateCorrelation calculateSpearmanCorrelation); # access these using 'use Modules::CorrelationWithRapa'
our @EXPORT_OK = qw(determineInfluentialPoints calculateCorrelation calculateSpearmanCorrelation); 
our %EXPORT_TAGS = (	);

use PDL; # required for determineInfluentialPoints
use Statistics::Distributions;

# determineInfluentialPoints
# This function receives a set of points as input and uses Cooks distance to determine
# influential points in the set. Any point with influence greater than 4 / # of data points
# is considered an influential outlier
# input: 2 equal size arrays that serve as the x and y coordinates for a set of points
# returns: a hash containing influential outliers.
sub determineInfluentialPoints{
	my ($array1,$array2) = @_;

	my $x = &computeXmatrix([$array1]);
	my $y = transpose(pdl[$array2]);

	my $xTrans = transpose($x);
	my $nObservations=$x->dim(1);
	my $XtXinverse = inv($xTrans x $x);
	my $betas = $XtXinverse x ($xTrans x $y);
	my $yHat = $x x $betas;
	my $E = $y - $yHat;
	# $prms = 
	# $adev = 
	# $rms = deviation from the mean (also known as the root-mean-square deviation, or the square root of the variance)
	my ($mean,$prms,$median,$min,$max,$adev,$rms) = stats($E);
	#  convert $rms to sum of squares residual (or sum of squares error it is sometimes called)
	my $SSres = $rms*$rms*$nObservations;

	# number of dependent variables
	my $k = ($x->dim(0)-1);

	# degrees of freedom for the residuals
	my $DFres = $x->dim(1)-$k-1;
	my $MSres = $SSres / $DFres;
	my $MSresXTX = $XtXinverse x $MSres;
	my $se = $MSresXTX->diagonal(0,1);

	# $hat = H = X(XTX)-1XT
	my $hat = ($x x $XtXinverse) x $xTrans;

	my $leverageVector = $hat->diagonal(0,1);
	# from wikipedia...
	# There are different opinions regarding what cut-off values to use for spotting highly influential points. 
	# A simple operational guideline of D_i>1 has been suggested [2]. Others have indicated that D_i>4/n, 
	# where n is the number of observations, might be used [3].
	# 
	#my $influenceThreshold = 1 > 4/$nObservations ? 1 : 4/$nObservations;
	my $influenceThreshold = 4/$nObservations;
	my %influentialPoints = ();
	for (my $i = 0; $i < $nObservations; $i++) {
		my $predY = $yHat->at(0,$i);
		my $residual =  $y->at(0,$i)-$predY;
		#my $modMSE = ($MSres - $residual*$residual/( (1-$leverageVector->at($i,0)) * $DFres ) ) * $DFres / ($DFres-1) ;
		#my $rStudent = $residual/sqrt($modMSE*(1-$leverageVector->at($i,0)));
		
		#  excel function not yet converted to perl
		# my $tTest = TDIST(abs($rStudent),$DFres,2);

		my $cooksDistance = ($residual*$residual/$MSres)*($leverageVector->at($i,0) / ((1-$leverageVector->at($i,0))*(1-$leverageVector->at($i,0))) ) / ($k+1);
		# DFFITS is a diagnostic meant to show how influential a point is in a statistical regression. 
		#my $dffits = $rStudent * sqrt($leverageVector->at($i,0)/(1-$leverageVector->at($i,0)));
		
		$influentialPoints{$i} = [$x->at(1,$i), $y->at(0,$i), $cooksDistance->at(0,0)] if ($cooksDistance > $influenceThreshold);
	}
	return \%influentialPoints;
}

# computeXmatrix
# input: 1 or more arrays of equal length
# given a set of arrays this function will convert it into a matrix where the first row is
# all 1 and each subsequent row contains the values of the arrays passed to this function. 
# The matrix is then transposed and returned (this means that now the first COLUMN is 1s and
# the subsequent columns contain that array values)
# return value: a PDL matrix
sub computeXmatrix{
	my ($arrays) = shift;
	my %lengths = ();
	foreach my $array(@{$arrays}){
		$lengths{scalar(@{$array})} = 1;
	}
	my @lengths = keys %lengths;
	if(@lengths > 1){
		warn "Error! Unequal array lengths!";
		warn "Lengths = ". join(@lengths, ", ");
		exit;
	}
	my @ones = (1) x $lengths[0];
	unshift(@{$arrays},\@ones);
	return transpose(pdl[@{$arrays}]);
}


#	calculateCorrelation
#	input: 2 arrays
#  calculate Pearson correlation coefficient and p-value (based on the t-distribution) of 2 arrays
#  return values: (calculated correlation, p-value)
sub calculateCorrelation{
	my ($array1, $array2) = @_;
	my $length = scalar(@{$array1});
	if($length > 2){
		# sum of square values works, but does not inform direction of correlation
		# my $ssxx = &sumSqValues($array1,$array1,$length);
		# my $ssyy = &sumSqValues($array2,$array2,$length);
		# my $ssxy = &sumSqValues($array1,$array2,$length);
		# # calculate the correlation coefficient r (sometimes also denoted R) is then defined by:
		# my $corrCo = sqrt($ssxy*$ssxy/($ssxx*$ssyy));

		# therefore use the method found here:
		# http://www.mathsisfun.com/data/correlation.html instead
		my $meanX = &mean($array1,$length); # scalar
		my $meanY = &mean($array2,$length); # scalar
		my ($sqrDiffSumX, $sqrDiffSumY, $diffProduct) = (0,0,0);
		for (my $i = 0; $i < $length; $i++) {
			my $aDiff = $array1->[$i] - $meanX;
			my $bDiff = $array2->[$i] - $meanY;
			$sqrDiffSumX += $aDiff*$aDiff;
			$sqrDiffSumY += $bDiff*$bDiff;
			$diffProduct += $aDiff * $bDiff;
		}
		my $corrCo = $diffProduct / sqrt($sqrDiffSumX * $sqrDiffSumY);
		# The p-value is basically whether the correlation coefficient is significantly different from 0 or not, 
		# so this is a t test: t=r/sqrt[(1-r²)/(N-2)] with n-2 degrees of freedom.
		my $t = abs($corrCo / sqrt((1-($corrCo*$corrCo) ) / ($length-2)));
		return ($corrCo,2*Statistics::Distributions::tprob(($length-2),$t));
	}
	return (0,1);
}

# calculateSpearmanCorrelation
#	input: 2 arrays
#  modified from gene boggs:
#  http://search.cpan.org/~gene/Statistics-RankCorrelation-0.1203/lib/Statistics/RankCorrelation.pm#spearman
#  given 2 arrays of representing a dataset of x and y (or whatever) values this function will return the 
#  spearman correlation and the associated p-value. The Spearman correlation is different from calculateCorrelation
#  subroutine in that it is non-parametric
#  return values: (calculated correlation, p-value)
sub calculateSpearmanCorrelation{
	my ($array1, $array2) = @_;
	my $length = scalar(@{$array1});
	if($length < 3){ return (0,1); }
	($array1,$array2) = &co_sort($array1,$array2);
	my ($xRank,$xTies) = &rank($array1);
	my ($yRank,$yTies) = &rank($array2);

	# Algorithm contributed by Jon Schutz <Jon.Schutz@youramigo.com>:
	my($x_sum, $y_sum) = (0, 0);
	$x_sum += $_ for @{$xRank};
	$y_sum += $_ for @{$yRank};
	my $n = $length;
	my $x_mean = $x_sum / $n;
	my $y_mean = $y_sum / $n;
	# Compute the sum of the difference of the squared ranks.
	my($x_sum2, $y_sum2, $xy_sum) = (0, 0, 0);
	for( 0 .. $length - 1 ) {
	    $x_sum2 += ($xRank->[$_] - $x_mean) ** 2;
	    $y_sum2 += ($yRank->[$_] - $y_mean) ** 2;
	    $xy_sum += ($xRank->[$_] - $x_mean) * ($yRank->[$_] - $y_mean);
	}
	return (1,0) if $x_sum2 == 0 || $y_sum2 == 0;
	my $corr = $xy_sum / sqrt($x_sum2 * $y_sum2);
	my $t = abs($corr / sqrt((1-($corr*$corr) ) / ($length-2)));
	return ($corr,2*Statistics::Distributions::tprob(($length-2),$t));
}

# rank
# input: an array
# taken from http://cpansearch.perl.org/src/GENE/Statistics-RankCorrelation-0.1203/lib/Statistics/RankCorrelation.pm
# Calculates and returns an array containing the ranks of the values in the array that is passed to it.
# return values: the ranks arrayref in a scalar context and include ties if called in a list context.
sub rank {
  my $u = shift;
  # Make a list of tied ranks for each datum.
  my %ties;
  push @{ $ties{ $u->[$_] } }, $_ for 0 .. @$u - 1;

  my ($old, $cur) = (0, 0);

  # Set the averaged ranks.
  my @ranks;
  for my $x (sort { $a <=> $b } keys %ties) {
    # Get the number of ties.
    my $ties = @{ $ties{$x} };
    $cur += $ties;
    if ($ties > 1) {
      # Average the tied data.
      my $average = $old + ($ties + 1) / 2;
      $ranks[$_] = $average for @{ $ties{$x} };
    }
    else {
      # Add the single rank to the list of ranks.
      $ranks[ $ties{$x}[0] ] = $cur;
    }
    $old = $cur;
  }

  # Remove the non-tied ranks.
  while( my( $k, $v ) = each %ties ) {  delete $ties{$k} unless @$v > 1;  }

  # Return the ranks arrayref in a scalar context and include ties
  # if called in a list context.
  return wantarray ? (\@ranks, \%ties) : \@ranks;
}

# co-sort
# input: 2 arrays
# taken from http://cpansearch.perl.org/src/GENE/Statistics-RankCorrelation-0.1203/lib/Statistics/RankCorrelation.pm
# Sort the vectors as two dimensional data-point pairs with B<u> values sorted first.
# returns: the sorted pair of arrays
sub co_sort {
  my( $u, $v ) = @_;
  return unless @$u == @$v;
  # Ye olde Schwartzian Transforme:
  $v = [
    map { $_->[1] }
      sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] }
        map { [ $u->[$_], $v->[$_] ] }
          0 .. @$u - 1
  ];
  # Sort the independent vector last.
  $u = [ sort { $a <=> $b } @$u ];
  return ($u, $v);
}

#  mean
#  input: an array and it's length
#  calculates the mean of the values in the array
#  returns: the mean
sub mean {
	my ($a,$length)=@_;
	my $sum = 0;
	for (my $i=0;$i<$length;$i++){$sum=$sum+$a->[$i];}
	return ($sum/$length);
}

1;