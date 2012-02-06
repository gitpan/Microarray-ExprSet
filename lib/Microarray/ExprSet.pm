package Microarray::ExprSet;

# general class of different microarray data structure
# contains three elements
# data matrix
# feature (gene) names array
# sample names array

use List::Vectorize;
use strict;

our @ISA = qw();

our $VERSION = "0.10";

1;


sub new {

	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $self = { @_ };
	bless($self, $class);
	return $self;
	
}


# probe name
sub feature {

	my $self = shift;
	
	return $self->{feature};

}

sub set_feature {
	
	my $self = shift;
	my $feature = shift;
	
	$self->{feature} = $feature;
	
	return 1;

}

# sample name
sub phenotype {

	my $self = shift;
	
	return $self->{phenotype};

}

sub set_phenotype {
	
	my $self = shift;
	my $phenotype = shift;
	
	$self->{phenotype} = $phenotype;
	
	return 1;
	
}

sub set_matrix {

	my $self = shift;
	my $matrix = shift;
	
	$self->{matrix} = $matrix;
	
	return 1;
	
}

sub matrix {

	my $self = shift;
	
	return $self->{matrix};

}


sub remove_empty_features {

	my $self = shift;
	
	my $old_feature = $self->feature;
	my $old_matrix = $self->matrix;
	
	my $new_feature;
	my $new_matrix;
	
	for(my $i = 0; $i < len($old_feature); $i ++) {
	
		if($old_feature->[$i] !~ /^\s*$/) {
			push(@$new_feature, $old_feature->[$i]);
			push(@$new_matrix, $old_matrix->[$i]);
		}
	}
	
	$self->set_feature($new_feature);
	$self->set_matrix($new_matrix);
	undef($old_feature);
	undef($old_matrix);

	return 1;
}

sub n_feature {

	my $self = shift;
	
	return len($self->feature);

}

sub n_sample {

	my $self = shift;
	
	return len($self->phenotype);

}

# using mean or median
sub unique_feature {
	
	my $self = shift;
	my $method = shift || "mean";
	
	my $fun;
	if($method eq "mean") {
		$fun = \&mean;
	}
	elsif($method eq "median") {
		$fun = \&median;
	}
	
	my $fh;
	my $feature = $self->feature;
	for(my $i = 0; $i < len($feature); $i ++) {
		if($fh->{$feature->[$i]}) {
			push(@{$fh->{$feature->[$i]}}, $i);
		}
		else {
			$fh->{$feature->[$i]}->[0] = $i;
		}
	}

	my $new_feature;
	my $new_matrix;
	my $matrix = $self->matrix;
	foreach my $f (keys %$fh) {
	
		my $index = $fh->{$f};
		push(@$new_feature, $f);
		
		if(len($index) == 1) {
			push(@$new_matrix, $matrix->[$index->[0]]);
		}
		else {
			my $new_array;
			for(my $i = 0; $i < len($matrix->[0]); $i ++) {
				my $tmp_array;
				for(my $j = 0; $j < len($index); $j ++) {
					
					push(@$tmp_array, $matrix->[$index->[$j]]->[$i]);
				
				}
				push(@$new_array, &$fun($tmp_array));
			}
			push(@$new_matrix, $new_array);
		}
	}
	
	$self->set_feature($new_feature);
	$self->set_matrix($new_matrix);
	undef($feature);
	undef($matrix);

	return 1;
}

__END__

=pod

=head1 NAME

Microarray::ExprSet - General class of microarray data

=head1 SYNOPSIS

  use Microarray::ExprSet;
  
  my $mat = [[1, 2, 3, 4, 5, 6],
             [7, 8, 9, 10, 11, 12],
			 [13, 14, 15, 16, 17, 18],
			 [19, 20, 21, 22, 23, 24],
			 [25, 26, 27, 28, 29, 30],
			 [31, 32, 33, 34, 35, 36]];
  my $probe = ["gene1", "gene2", "gene2", "gene3", "", "gene4"];
  my $sample = ["c1", "c1", "c1", "c2", "c2", "c2"];
  
  my $expr = Microarray::ExprSet->new();
  $expr->set_matrix($mat);
  $expr->set_feature($probe);
  $expr->set_phenotype($sample);
  
  # do some preprocess
  $expr->remove_empty_features();
  # combine duplicated genes
  $expr->unique_feature("mean");  # you can use "median" too
  
  # now you can get content of the object
  my $new_mat = $expr->matrix;
  my $new_probe = $expr->feature;
  my $new_sample = $expr->phenotype;
  my $n_probe = $expr->n_feature;
  my $n_sample = $expr->n_phenotype;

=head1 DESCRIPTION

The C<Microarray::ExprSet> class object describes the data structure of microarray
data. It contains three elements: 1) data matrix
that stores the expression value; 2) array of features that are the probe names
or gene IDs; 3) array of phenotypes that are the settings of samples (e.g. control vs
treatment). Other information about the microarray experiment such as the protocal
or sample preparation is not included in the object.

Usually the C<Microarray::ExprSet> object is created by other modules such as
L<Microarray::GEO::SOFT>.

=head2 Subroutines

=over 4

=item C<new()>

Initial a C<Microarray::ExprSet> class object.

=item C<$expr-E<gt>set_matrix(MATRIX)>

Argument is the expression value matrix which is stored in an array reference of array
references.

=item C<$expr-E<gt>set_feature([ARRAY REF])>

Set the feature names. The length of features should be equal to the number of 
rows of the expression value matrix.

=item C<$expr-E<gt>set_phenotype([ARRAY REF])>

Set the phenotype names. The length of phenotypes should be equal to the number
of columns of the expression value matrix.

=item C<$expr-E<gt>matrix>

Get expression value matrix

=item C<$expr-E<gt>feature>

Get feature names

=item C<$expr-E<gt>phenotype>

Get phenotype names

=item C<$expr-E<gt>n_feature>

Get the number of features

=item C<$expr-E<gt>n_phenotype>

Get the number of phenotypes

=item C<$expr-E<gt>remove_empty_features>

Some features may not have names, so it is necessary to eliminate these features
without any names.

=item C<$expr-E<gt>unique_feature(mean | median)>

It is usually that features are measured repeatly, while some analysis procedures need
unified features. The second argument can be set to choose the method for multiple
feature merging.

=back

=head1 AUTHOR

Zuguang Gu E<lt>jokergoo@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright 2012 by Zuguang Gu

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.12.1 or,
at your option, any later version of Perl 5 you may have available.

=head1 SEE ALSO

L<Microarray::GEO::SOFT>

=cut
