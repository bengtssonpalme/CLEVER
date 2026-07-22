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

open (VAR, ">$databaseFolder/CLEVER.variants.tsv");
open (FAM, ">$databaseFolder/CLEVER.families.tsv");
open (LIN, ">$databaseFolder/CLEVER.lineages.tsv");

print STDERR "Reading changes...\n";
while ($orgline = <ORG>) {
    chomp($orgline);
    $newline = <NEW>;
    chomp($newline);
    if ($orgline eq $newline) {
        ($nid, $ngene, $nest, $nmob, $nver, $nsource, $ndbid, $dbFile) = split('\t', $newline);
        if ($ndbid ne "") {
            $new_id{$ndbid} = $nid;
        }
        if ($dbFile eq "V") {
            print VAR $newline . "\n";
        }
        if ($dbFile eq "F") {
            print FAM $newline . "\n";
        }
        if ($dbFile eq "L") {
            print LIN $newline . "\n";
        }
        next;
    } else {
        # C1|aac(6')-Ib|aac6p-25_4|ECV|ResFinder|M23634.1	aac6p	E	M	V	ResFinder	C1|aac(6')-Ib|aac6p-25_4|ECV|ResFinder|M23634.1
        ($oid, $ogene, $oest, $omob, $over, $osource, $odbid) = split('\t', $orgline);
        ($nid, $ngene, $nest, $nmob, $nver, $nsource, $ndbid, $dbFile) = split('\t', $newline);
        if ($ndbid ne "") {
            $new_id{$ndbid} = $nid;
        }
        if ($odbid ne "") {
            $new_id{$odbid} = $nid;
        }
        if ($ngene ne "EXCLUDE") {
        
            if ($dbFile eq "V") {
                print VAR $newline . "\n";
            }
            if ($dbFile eq "F") {
                print FAM $newline . "\n";
            }
            if ($dbFile eq "L") {
                print LIN $newline . "\n";
            }
        }
    }
    if ($ogene ne $ngene) {
        $changes_gene{$oid} = "$ogene\n$ngene";
    }
    if ($ngene eq "EXCLUDE") {
        $changes_gene{$oid} = "$ogene\n$ngene";
    }
    if ($odbid ne $ndbid) {
        $changes_dbid{$oid} = "$odbid\n$ndbid";
    }
    if ($oid ne $nid) {
        $changes_id{$oid} = "$oid\n$nid";
    }
}

@changeDB_gene = keys(%changes_gene);
@changeDB_dbid = keys(%changes_dbid);
@changeDB_id = keys(%changes_id);
@updateMap = keys(%new_id);

foreach $file (@changeFiles) {
    print STDERR "\nProcessing file $file...\n";
    $lines = 0;
    $last_permille = 0;
    $last_percent = 0;
    
    $exclusion = 0;
    `mv $file $file.temp`;
    chomp($totalLines = `wc -l $file.temp`);
    open (INFILE, "$file.temp");
    open (OUTFILE, ">$file");
    while ($line = <INFILE>) {
        chomp($line);
        $lines++;
        $permille = int($lines * 100 / $totalLines);
        if ($permille > $last_permille) {
            print STDERR ".";
            $last_permille = $permille;
        }
        $percent = int($lines * 10 / $totalLines) * 10;
        if ($percent > $last_percent) {
            print STDERR " $percent% ";
            $last_percent = $percent;
        }
        foreach $gene_change (@changeDB_gene) {
            if ($line =~ m/\Q$gene_change\E/) {
                ($ogene,$ngene) = split('\n', $changes_gene{$gene_change});
                if ($ngene eq "EXCLUDE") {
                    if ($line =~ m/^>/) {
                        $exclusion = 1;
                    }
                } else {
                    $line =~ s/\Q$ogene\E/$ngene/g;
                }
            }
        }
        foreach $dbid_change (@changeDB_dbid) {
            if ($line =~ m/\Q$dbid_change\E/) {
                ($odbid,$ndbid) = split('\n', $changes_dbid{$dbid_change});
                $line =~ s/\Q$odbid\E/$ndbid/g;
            }
        }
        foreach $id_change (@changeDB_id) {
            if ($line =~ m/\Q$id_change\E/) {
                ($oid,$nid) = split('\n', $changes_id{$id_change});
                $line =~ s/\Q$oid\E/$nid/g;
            }
        }
        if ($line =~ m/^>/) {
            $exclusion = 0;
        }
        if ($exclusion == 0) {
            print OUTFILE $line . "\n";
        }
    }
    close OUTFILE;
    close INFILE;
}
foreach $file (@changeFiles) {
    `rm $file.temp`;
}

close ORG;
close NEW;
close LIN;
close FAM;
close VAR;

print STDERR "\nUpdating mapping IDs...\n";

`cp $databaseFolder/CLEVER.variant-family-mapping.txt $databaseFolder/CLEVER.variant-family-mapping-old.txt`;

open (MAP, "$databaseFolder/CLEVER.variant-family-mapping-old.txt");
open (NEWMAP, ">$databaseFolder/CLEVER.variant-family-mapping.txt");
$lines = 0;
$last_permille = 0;
$last_percent = 0;
chomp($totalLines = `wc -l $databaseFolder/CLEVER.variant-family-mapping-old.txt`);

while ($line = <MAP>) {
    $lines++;
    $permille = int($lines * 100 / $totalLines);
    if ($permille > $last_permille) {
        print STDERR ".";
        $last_permille = $permille;
    }
    $percent = int($lines * 10 / $totalLines) * 10;
    if ($percent > $last_percent) {
        print STDERR " $percent% ";
        $last_percent = $percent;
    }
    chomp($line);
    foreach $mapchange (@updateMap) {
        if ($mapchange ne "") {
            if ($line =~ m/\Q$mapchange\E/) {
                $replacement = $new_id{$mapchange};
                $line =~ s/\Q$mapchange\E/$replacement/g;
            }
        }
    }
    print NEWMAP $line . "\n";
}
close MAP;
close NEWMAP;
print STDERR "\nFinished!\n";

