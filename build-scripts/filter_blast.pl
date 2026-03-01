#!/usr/bin/perl

$file = shift;
$id_cutoff = shift;

%hits = {};
%hitCount = {};
open (INPUT, $file);
while ($line = <INPUT>) {
    chomp($line);
    (@items) = split("\t", $line);
    $arg = @items[0];
    $plasmid = @items[1];
    $id = @items[2];
    if ($id >= $id_cutoff) {
        if (not(defined($hits{$arg}))) {
            $hits{$arg} = $plasmid;
            $hitCount{$arg} = 1;
        } else {
            $hits{$arg} = $hits{$arg} . "," . $plasmid;
            $hitCount{$arg}++;
        }
    }
}
close INPUT;

open (HITS, ">$file.hits.txt");
open (HITLIST, ">$file.hit_ids.txt");
open (HITCOUNT, ">$file.hit_count.txt");
foreach $arg (sort(keys(%hits))) {
    print HITS $arg . "\n";
    print HITLIST $arg . "\t" . $hits{$arg} . "\n";
    print HITCOUNT $arg . "\t" . $hitCount{$arg} . "\n";
}
close HITS;
close HITLIST;
close HITCOUNT;