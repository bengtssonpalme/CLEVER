#!/usr/bin/perl

$diamondout = shift;
open (DIAMONDOUT, $diamondout);
while ($line = <DIAMONDOUT>) {
    chomp($line);
    (@items) = split("\t", $line);
    $query = @items[0];
    if (defined($byQuery{$query})) {
        $byQuery{$query} = $byQuery{$query} . "\n" . $line;
    } else {
        $byQuery{$query} = $line;
    }
}
close DIAMONDOUT;

foreach $query (keys(%byQuery)) {
    (@lines) = split("\n", $byQuery{$query});
    $maxline = "";
    $maxid = 0;
    $maxcluster = "";
    foreach $line (@lines) {
        chomp($line);
        (@items) = split("\t", $line);
        $cluster = @items[1];
        $identity = @items[2];
        if ($identity > $maxid) {
            if ($cluster =~ m/C[0-9]*\|/) {
                $maxid = $identity;
                $maxline = $line;
                $maxcluster = $cluster;
            } else {
                if ($maxcluster =~ m/C[0-9]*\|/) {
                    $maxid = $maxid;
                    $maxline = $maxline;
                    $maxcluster = $maxcluster;
                } else {
                    $maxid = $identity;
                    $maxline = $line;
                    $maxcluster = $cluster;
                }
            }
            
        }
    }
    if ($maxline ne "") {
        print $maxline . "\n";
    }
}