use POSIX;
use strict;

my($file, $line, @item, @names, $i, $j, @temp);
$file = $ARGV[0];
#`perl -p -i -e 's/^/\t/g $file`; 
open(IN, "<$file");
while(!eof(IN)){
	$line = readline *IN;
#	print $line."\n";
	chomp $line;
	@temp = split/\t/, $line;
	$names[$i][0] = $temp[0];
	$names[$i][2] = $temp[1];
	$names[$i][1]= $temp[2];
#	print $names[$i][0]."\t".$names[$i][1]."\n";
	$i++;
}
open(OUT, ">Step_1.txt");
$i = 0;
$j = 1;
#print scalar(@names)."woof\n";
while ($i < scalar(@names)){
#	print $i."\n";
	$j = $names[$i][2];
	while($j < $names[$i][1]){
		print OUT $names[$i][0]."\t".$j."\n";
		$j++;
	}
	print OUT $names[$i][0]."\t".$j."\n";
	$i++;
}

