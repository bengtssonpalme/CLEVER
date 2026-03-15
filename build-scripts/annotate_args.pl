#!/usr/bin/perl

$file = shift;
$blacklistfile = shift;
$cleverversion = shift;

%clusters = {};
%foundARGs = {};
@blacklist = ();

open (BLACKLIST, $blacklistfile);
while ($line = <BLACKLIST>) {
    chomp($line);
    push(@blacklist,$line);
}
close BLACKLIST;


open (INPUT, $file);
while ($line = <INPUT>) {
    chomp($line);
    (@items) = split("\t", $line);
    $arg = @items[0];
    if (not(defined($foundARGs{$arg}))) {
        $foundARGs{$arg} = 1;
        $cluster = @items[1];
        $id = @items[2];


        if (defined($clusters{$cluster})) {
            $clusters{$cluster} = $clusters{$cluster} . "\n" . $arg;
        } else {
            $clusters{$cluster} = $arg;
        }
    
    }
    
}
close INPUT;

open (GENELIST, ">$file.genelist.txt");
open (MAPPING, ">$file.mapping.txt");
open (LOG, ">$file.logfile.txt");
@families = sort(keys(%clusters));

foreach $family (@families) {
    #print $family . ":\n";
    #print $clusters{$family} . "\n";
    #print "--------------\n";

    (@args) = split('\n',$clusters{$family});
    undef(@argNames);
    undef(@accessions);
    #unshift(@args, $family);
    $addSelf = 1;
    foreach $arg (@args) {
        if ($arg eq $family) {
            $addSelf = 0;
        }
    }
    if ($addSelf == 1) {
        unshift(@args, $family);
    }
    foreach $arg (@args) {
        if (substr($arg, 0, 2) =~ m/C[0-9]/) {
            ## This is an existing CLEVER gene
            ($version, $argName, $cleverID, $attr, $source, $accession) = split('\|',$arg);
            push(@argNames, $argName);
            push(@accessions, $accession);
        } else {
            ($source, $rest) = split("-",$arg);
            $rest = $arg;
            $rest =~ s/^[^-]*-//;
            if ($source eq "ResFinder") {
                #ResFinder-qnrB90_1_MG182074_1
                ($name,$variant,$accession,$version) = split('_', $rest);
                $argName = $name;
                if ($argName =~ m/^[A-Z][A-Z]*-/) {
                    $argName = "bla" . $argName;
                } else {
                    if ($argName =~ m/^[A-Z][A-Z][A-Z]\(/) {
                        $argName = lc(substr($argName, 0, 3)) . substr($argName, 3);
                    } else {
                        if ($argName =~ m/^[A-Z]/) {
                            $argName = lc(substr($argName, 0, 1)) . substr($argName, 1);
                        }
                    }
                }
                $full_accession = $accession . "." . $version;
                push(@argNames, $argName);
                push(@accessions, $full_accession);
            }
            if ($source eq "CARD") {
                #CARD-gb|ADQ43424.1|ARO:3002746|QnrB31
                ($gb,$accession,$aro,$name) = split('\|', $rest);
                $argName = $name;
                if ($argName =~ m/^[A-Z][A-Z]*-/) {
                    $argName = "bla" . $argName;
                } else {
                    if ($argName =~ m/^[A-Z][A-Z][A-Z]\(/) {
                        $argName = lc(substr($argName, 0, 3)) . substr($argName, 3);
                    } else {
                        if ($argName =~ m/^[A-Z]/) {
                            $argName = lc(substr($argName, 0, 1)) . substr($argName, 1);
                        }
                    }
                }
                $full_accession = $accession;
                push(@argNames, $argName);
                push(@accessions, $full_accession);
            }
            if ($source eq "ResFinderFG") {
                #ResFinderFG-beta_lactamase|KU545081.1|feces|ATM
                ($class,$accession,$aro,$source,$antibiotic) = split('\|', $rest);
                $argName = "!" . $class;
                $full_accession = $accession;
                push(@argNames, $argName);
                push(@accessions, $full_accession);
            }
            if (($source eq "fARGene") || ($source eq "Inda-Diaz_2023")) {
                #fARGene-20k_concatenated-long-orfs_SFKQ01000019.1_seq1@@@methyltransferase_grp2_1
                ($seqinfo, $class)= split('@@@', $rest);
                ($junk1,$junk2,$accession,$variant) = split('_', $seqinfo);
                $argName = "!" . $class;
                $full_accession = $accession;
                push(@argNames, $argName);
                push(@accessions, $full_accession);
            }
        }
    }
## Wanted output:
## C[version]_[class]-[#]|[EL][MC][VP]|[~][gene]|[source]|[accession]
## 
## Ex:
## C1_van1|EMV|VanHDX|AY489045_3|ResFinder
## C1_blaC1-1|LCP|~C1_blaC1-1|Inda-Diaz|AAAAQY010000096.1

    %disagreement = {};

    if (scalar(@argNames) < 3) {
        $agreedName = @argNames[0];
        if ((@argNames[0] eq "") || (substr(@argNames[0],0,1) eq "!")) {
            if ((@argNames[1] ne "") || (substr(@argNames[1],0,1) eq "!")) {
                $agreedName = @argNames[1];
                print STDERR "Family $family lacked a gene name for the reference sequence. Settled on $agreedName for gene. You might want to check this entry manually.\n";
                print LOG "$family : Family $family lacked a gene name for the reference sequence. Settled on $agreedName for gene. You might want to check this entry manually.\n";
            } else {
                $agreedName = "!";
                if ($family !~ m/HASH/) {
                    print STDERR "Family $family lacked any valid gene names for any sequence. You might want to check this entry manually.\n";
                    print LOG "$family : Family $family lacked any valid gene names for any sequence. You might want to check this entry manually.\n";
                }
            }
        }
    } else {
        undef(@agreedArray);
        $problem = 0;
        for ($c = 0; $c < length(@argNames[0]); $c++) {
            if ($problem == 1) {
                last;
            }
            $mainName = @argNames[0];
            for ($i = 0; $i < scalar(@argNames); $i++) {
                if (defined($disagreement{"@argNames[$i]"})) {
                    $disagreement{"@argNames[$i]"} = $disagreement{"@argNames[$i]"} + 1; 
                } else {
                    $disagreement{"@argNames[$i]"} = 1;
                }
                if ($i > 0) {
                    if (@argNames[$i] =~ m/\Q$mainName\E/) {
                        next;
                    }
                    if (length(@argNames[$i]) < 2) {
                        next;
                    }
                }
                if (substr(@argNames[$i], 0 ,1) eq "!") {
                    next;
                } else {
                    if (substr($mainName, $c, 1) eq substr(@argNames[$i], $c, 1)) {
                        @agreedArray[$c] = substr($mainName, $c, 1);
                    } else {
                        @agreedArray[$c] = "";
                        $problem = 1;
                        last;
                    }
                }
            }
        }
        if (($mainName =~ m-/-) || ($problem == 1))  {
            $agreedName = "";
            @oldAgreedArray = @agreedArray;
            @agreedArray[0] = "?";
            @agreedArray[1] = "/";
            for ($c = 1; $c < scalar(@oldAgreedArray); $c++) {
                @agreedArray[$c+1] = @oldAgreedArray[$c-1];
                $agreedName = $agreedName . @agreedArray[$c+1];
            }
            if ($agreedName eq $mainName) {
                for ($c = 0; $c < length($mainName); $c++) {
                    @agreedArray[$c] = substr($mainName, $c, 1);
                }
            }
        }
        if (scalar(@agreedArray) > 2) {
            $agreedName = "";
            for ($c = 0; $c < scalar(@agreedArray); $c++) {
                $agreedName = $agreedName . @agreedArray[$c];
            }
            if (length($agreedName) < 2) {
                $agreedName = "!";
            }
        } else {
            $agreedName = "!";
        }
    }
    if (($agreedName =~ m/\([^)]*$/)) {
        $agreedName = $agreedName . ")";
    }
    if ($agreedName =~ m/^[?!]\//) {
        $testName = $agreedName;
        $testName =~ s/^[?!]\///;
        if ($mainName =~ m/\Q$testName\E/) {
            if (substr($testName, -1) eq "_") {
                $agreedName = substr($testName, 0, -1);
            } else {
                $agreedName = $mainName;
            }
        }
    }
    if (substr($agreedName, -1) =~ m/[-_.]/) {
        $agreedName = substr($agreedName, 0, -1);
    }
    $change = 0;
    if ($agreedName eq "!") {
        if (substr($mainName, 0 , 1) ne "!") {
            $agreedName = $mainName;
        } else {
            $agreedName = "";
        }

        foreach $name (reverse sort {$disagreement{$a} <=> $disagreement{$b}} keys %disagreement) {
            if ($name !~ m/HASH/) {
                if (substr($name, 0, 1) ne "!") {
                    if ($mainName =~ m/\Q$name\E/) {
                        #print "Do not add $name to $mainName\n";
                        $change = 1;
                    } else {
                        $agreedName = $agreedName . "/" . $name;
                        #print "Add $name to $mainName\n";
                        $change = 0;
                    }
                }
            }
        }
    }
    if ($change == 1) {
        $agreedName = $mainName;
    }
    (@argsincluster) = split('\n',$clusters{$family});
    $blacklisted = 0;
    foreach $gene (@argsincluster) {
        foreach $item (@blacklist) {
            if ($item ne "") {
                if ($gene =~ m/$item/) {
                    $blacklisted = 1;
                }
            }
        }
    }
    
    if (($agreedName =~ m/[?!]/) || ($agreedName eq "")) {
        if ($family !~ m/HASH/) {
            print LOG "$family : Family $family has inconistent gene name.\n";
            print STDERR "--------------\n";
            print STDERR "Inconsistent naming of family:\t". $family . "\n";
            print STDERR "Suggested name:             Original name:\n";
            print STDERR "==> " . $agreedName . "   " . $mainName . "\n";
            print STDERR "These are the genes in the cluster:\n";
            print STDERR $family . "\t";
            foreach $gene (@argsincluster) {
                print STDERR $gene . "\t";
            }
            print STDERR "\n";
            print STDERR "--------------\n";
            print STDERR "Please enter a solution to solve this issue. Enter:\n";
            print STDERR "   (1) new gene name\n";
            print STDERR "   (2) empty to accept suggested name above - $agreedName\n";
            print STDERR "   (3) 'O' or 'R' or 'ORIGINAL' or 'REVERT' to keep original name - $mainName\n";
            print STDERR "   (4) 'E' or 'X' or 'B' or 'EXCLUDE' or 'BLACKLIST' to exclude entry\n";
            $solution = <>;
            chomp($solution);
            if (length($solution) < 2) {
                $agreedName = $agreedName;
            } else {
                $agreedName = $solution;
            }
            if ((uc($solution) eq "B") || (uc($solution) eq "E") || (uc($solution) eq "X")) {
                $agreedName = "EXCLUDED: " . $agreedName;
                print STDERR "Excluded gene.\n";
                print LOG "$family : Solved by exluding gene.\n";
            }
            if ((uc($solution) eq "O") || (uc($solution) eq "R")) {
                $agreedName = $mainName;
                print STDERR "Reverted to original name $agreedName\n";
                print LOG "$family : Solved by reverting to original gene name: $agreedName\n";
            }
            
            if (($solution eq "BLACKLIST") || ($solution eq "EXCLUDE")) {
                $agreedName = "EXCLUDED: " . $agreedName;
                print STDERR "Excluded gene.\n";
                print LOG "$family : Solved by exluding gene.\n";
            }

            if (($solution eq "ORIGINAL") || ($solution eq "REVERT")) {
                $agreedName = $mainName;
                print STDERR "Reverted to original name $agreedName\n";
                print LOG "$family : Solved by reverting to original gene name: $agreedName\n";
            }

            if ($agreedName eq $solution) {
                print STDERR "Updated name to $agreedName\n";
                print LOG "$family : Solved by updating gene name to $agreedName.\n";
            }
            
            print STDERR "--------------\n";
        }
    }
    if ($blacklisted == 1) {
        print GENELIST $family . "\t" . "BLACKLISTED: " . $agreedName . "\t";
        print LOG "$family : Famuily $family was on a blacklist and has been excluded.\n";
    } else {
        print GENELIST $family . "\t" . $agreedName . "\t";
    }
    foreach $gene (@argsincluster) {
        print GENELIST $gene . ",";
        print MAPPING $gene . "\t" . $agreedName ."\n";
    }
    print GENELIST "\n";
}
close GENELIST;
close MAPPING;

sub getClass {
    $geneName = shift;
    $class = "";
    return $class;
}