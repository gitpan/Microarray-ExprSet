use ExtUtils::MakeMaker;

WriteMakefile(
	NAME          => 'Microarray::ExprSet',
	VERSION_FROM  => 'lib/Microarray/ExprSet.pm',
	PREREQ_PM     => {
		'List::Vectorize' => "",
		'Test::More' => "",
	},
	AUTHOR        => 'Zuguang Gu <jokergoo@gmail.com>',
	ABSTRACT      => 'General class of microarray data',
	LICENSE       => 'perl',
);
