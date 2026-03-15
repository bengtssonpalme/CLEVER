#!/usr/bin/perl

$genelist = shift;
$oldannotation = shift;
$manformat = shift;

%annotations = {};
%foundARGs = {};

open (LOG, ">$genelist.log.txt");

if ($manformat eq "1") {
    print STDERR "Basing annotation on intermediate format. Mobility and source information might be lost!\n";
    print LOG "Basing annotation on intermediate format. Mobility and source information might be lost!\n";
}

open (ANNOT, $oldannotation);
while ($line = <ANNOT>) {
    $blacklisted = 0;
    chomp($line);
    if ($manformat eq "1") {
        ($cleverID, $autogeneName, $geneName, $autoclass, $class, $est, $ver, $seqIDs) = split('\t', $line);
        (@seqIDs) = split(',', $seqIDs);
        if ($geneName =~ m/BLACKLISTED/) {
            $blacklisted = 1;
        }
        if ($geneName =~ m/EXCLUDED/) {
            $blacklisted = 1;
        }
        if ($autogeneName =~ m/BLACKLISTED/) {
            $blacklisted = 1;
        }
        if ($autogeneName =~ m/EXCLUDED/) {
            $blacklisted = 1;
        }
        if ($blacklisted == 1) {
            $annotations{$cleverID} = "EXCLUDED: $cleverID\t$class\t$est\t?\t$ver\t?\t$cleverID";
        } else {
            $annotations{$cleverID} = "$cleverID\t$class\t$est\t?\t$ver\t?\t$cleverID";
        }
        foreach $id (@seqIDs) {
            if ($blacklisted == 1) {
                $annotations{$id} = "EXCLUDED: $id\t$class\t$est\t?\t$ver\t?\t$id";
            } else {
                $annotations{$id} = "$id\t$class\t$est\t?\t$ver\t?\t$id";
            }
        }
    } else {
        ($cleverID, $class, $est, $mob, $ver, $source, $orgID) = split('\t', $line);
        $annotations{$cleverID} = $line;
    }
}
close ANNOT;

print "Original cluster name\tAutomated gene name\tCorrected gene name\tSuggested Class\tClass\tEstablished\tVerified\tSequences\n";
open (GENELIST, $genelist);
while ($line = <GENELIST>) {
    $blacklisted = 0;
    chomp($line);
    ($cleverID, $geneName, $seqs) = split('\t', $line);
    if ($geneName =~ m/BLACKLISTED/) {
        $blacklisted = 1;
    }
    if ($geneName =~ m/EXCLUDED/) {
        $blacklisted = 1;
    }
    if (defined($annotations{$cleverID})) {
        $annotation = $annotations{$cleverID};
    } else {
        (@entries) = split(',', $seqs);
        foreach $entry (@entries) {
            $annotation = $annotations{$entry};
            if ($annotation ne "") {
                last;
            }
        }
    }
    if ($annotation eq "") {
        if ($cleverID =~ m/^C[0-9]/) {
            ($version, $geneName, $cleverName, $attr, $source, $orgID) = split('\|', $cleverID);
            $est = substr($attr, 0, 1);
            $mob = substr($attr, 1, 1);
            $ver = substr($attr, 2, 1);
            ($class,$unused) = split('-', $cleverName);
            $annotation = "$cleverID\t$class\t$est\t$mob\t$ver\t$source\t$orgID";
        } else {
            print STDERR "Could not resolve annotation for $cleverID. You should check this entry manually.\n";
            print LOG "$cleverID : Could not resolve annotation for $cleverID. You should check this entry manually.\n";
        }
    }
    
    if ($class eq "") {
        print STDERR "Could not resolve class for $cleverID. You should check this entry manually.\n";
        print LOG "$cleverID : Could not resolve class for $cleverID. You should check this entry manually.\n";
    }

    ($cleverIDB, $class, $est, $mob, $ver, $source, $orgID) = split('\t', $annotation);
    if ($annotation =~ m/BLACKLISTED/) {
        $blacklisted = 1;
    }
    if ($annotation =~ m/EXCLUDED/) {
        $blacklisted = 1;
    }
    if ($blacklisted == 1) {
        print "$cleverID\tBLACKLISTED: $geneName\tEXCLUDED\t$class\t$class\t$est\t$ver\t$seqs\n";
        print LOG "$cleverID : BLACKLISTED. Gene name: $geneName\n";
    } else {
        print "$cleverID\t$geneName\t$geneName\t$class\t$class\t$est\t$ver\t$seqs\n";
    }
}
close GENELIST;
close LOG;

