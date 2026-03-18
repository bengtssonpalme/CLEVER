#!/usr/bin/perl

$originalfile = shift;
$updatedfile = shift;
$databaseFolder = shift;

push(@changeFiles, "$databaseFolder/CLEVER.families.faa");
push(@changeFiles, "$databaseFolder/CLEVER.lineages.faa");
push(@changeFiles, "$databaseFolder/CLEVER.variant-family-mapping.txt");
push(@changeFiles, "$databaseFolder/CLEVER.variants.faa");

push(@tsvFiles, "$databaseFolder/CLEVER.families.tsv");
push(@tsvFiles, "$databaseFolder/CLEVER.lineages.tsv");
push(@tsvFiles, "$databaseFolder/CLEVER.variants.tsv");

## This is for debugging
foreach $file (@changeFiles) {
    `cp $file $file.original`;
}

open (ORG, "$originalfile");
open (NEW, "$updatedfile");

open (VAR, ">$databaseFolder/CLEVER.variants.tsv.new");
open (FAM, ">$databaseFolder/CLEVER.families.tsv.new");
open (LIN, ">$databaseFolder/CLEVER.lineages.tsv.new");

while ($orgline = <ORG>) {
    chomp($orgline);
    $newline = <NEW>;
    chomp($newline);
    if ($orgline eq $newline) {
        ($nid, $ngene, $nest, $nmob, $nver, $nsource, $ndbid, $dbFile) = split('\t', $newline);
        if ($dbFile eq "V") {
            print VAR $newline;
        }
        if ($dbFile eq "F") {
            print FAM $newline;
        }
        if ($dbFile eq "L") {
            print LIN $newline;
        }
        next;
    } else {
        # C1|aac(6')-Ib|aac6p-25_4|ECV|ResFinder|M23634.1	aac6p	E	M	V	ResFinder	C1|aac(6')-Ib|aac6p-25_4|ECV|ResFinder|M23634.1
        ($oid, $ogene, $oest, $omob, $over, $osource, $odbid) = split('\t', $orgline);
        ($nid, $ngene, $nest, $nmob, $nver, $nsource, $ndbid, $dbFile) = split('\t', $newline);
        
        if ($dbFile eq "V") {
            print VAR $newline;
        }
        if ($dbFile eq "F") {
            print FAM $newline;
        }
        if ($dbFile eq "L") {
            print LIN $newline;
        }

        if ($oid ne $nid) {
            foreach $file (@changeFiles) {
                `mv $file $file.temp`;
                open (INFILE, "$file.temp");
                open (OUTFILE, ">$file");
                #print STDERR $oid . " => " . $nid . "\n";
                while ($line = <INFILE>) {
                    if ($line =~ m/\Q$oid\E/) {
                        if (substr($line, 0 , 1) eq ">") {
                            $repline = $newline;
                            $repline =~ s/\t/ /g;
                            print OUTFILE $repline . "\n";
                        } else {
                            $line =~ s/\Q$oid\E/$nid/g;
                            print OUTFILE $line;
                        }
                    } else {
                        print OUTFILE $line;
                    }
                }
                close OUTFILE;
                close INFILE;
            }
        }
        if ($ogene ne $ngene) {
            foreach $file (@changeFiles) {
                `mv $file $file.temp`;
                open (INFILE, "$file.temp");
                open (OUTFILE, ">$file");
                while ($line = <INFILE>) {
                    if ($line =~ m/\Q$nid\E/) {
                        $line =~ s/\Q$ogene\E/$ngene/g;
                    }
                    print OUTFILE $line;
                }
                close OUTFILE;
                close INFILE;
            }
        }
        if ($odbid ne $ndbid) {
            foreach $file (@changeFiles) {
                `mv $file $file.temp`;
                open (INFILE, "$file.temp");
                open (OUTFILE, ">$file");
                while ($line = <INFILE>) {
                    $line =~ s/\Q$odbid\E/$ndbid/g;
                    print OUTFILE $line;
                }
                close OUTFILE;
                close INFILE;
            }
        }
    }
}

foreach $file (@changeFiles) {
    `rm $file.temp`;
}

close ORG;
close NEW;
close LIN;
close FAM;
close VAR;