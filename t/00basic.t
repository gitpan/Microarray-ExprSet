use strict;
use warnings;

use Test::More tests => 2;

BEGIN {
	use_ok('List::Vectorize') or BAIL_OUT "Unable to load List::Vectorize";
	use_ok('Microarray::ExprSet') or BAIL_OUT "Unable to load Microarray::ExprSet";
}
