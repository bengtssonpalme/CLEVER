#!/usr/bin/perl

$fastafile = shift;
$annotationfile = shift;
$mobilefile = shift;
$blacklistfile = shift;
$cleverversion = shift;

%repMap = {};
%argMap = {};
%classMap = {};
%familyNum = {};
%varNum = {};
%familyNumMap = {};
%varNumMap = {};
%estMap = {};
%mobMap = {};
%verMap = {};

open (ANNOT, $annotationfile);
while ($line = <ANNOT>) {
    chomp($line);
    ($clusterRep, $autoARG, $curatedARG, $autoclass, $class, $est, $ver, $clusterentries) = split('\t', $line);
    (@entries) = split(',', $clusterentries);
    if (defined($familyNum{$class})) {
        $familyNum{$class} = $familyNum{$class} + 1;
    } else {
        $familyNum{$class} = 1;
    }
    foreach $entry (@entries) {
        if (defined($varNum{"$class-$familyNum{$class}"})) {
            $varNum{"$class-$familyNum{$class}"} = $varNum{"$class-$familyNum{$class}"} + 1;
        } else {
            $varNum{"$class-$familyNum{$class}"} = 1;
        }
        $repMap{$entry} = $clusterRep;
        $argMap{$entry} = $curatedARG;
        $classMap{$entry} = $class;
        $estMap{$entry} = $est;
        $verMap{$entry} = $ver;
        $familyNumMap{$entry} = $familyNum{$class};
        $varNumMap{$entry} = $varNum{"$class-$familyNum{$class}"};
    }

}
close ANNOT;

open (MOBILELIST, $mobilefile);
while ($line = <MOBILELIST>) {
    chomp($line);
    $mobMap{$line} = "M";
}
close MOBILELIST;

@blacklist = ();

open (BLACKLIST, $blacklistfile);
while ($line = <BLACKLIST>) {
    chomp($line);
    push(@blacklist,$line);
}
close BLACKLIST;

$outputEntry = 0;
open (FASTA, $fastafile);
while ($line = <FASTA>) {
    chomp($line);
    if (substr($line, 0 ,1) eq ">") {
        $blacklisted = 0;
        foreach $item (@blacklist) {
            if ($item ne "") {
                if ($line =~ m/$item/) {
                    $blacklisted = 1;
                }
            }
        }
        if ($blacklisted == 1) {
            $outputEntry = 0;
            next;
        }
        ($arg) = split(' ', $line);
        $arg = substr($arg, 1);
        ($source, $rest) = split("-",$arg);
        $rest = $arg;
        $rest =~ s/^[^-]*-//;
        if (substr($arg, 0, 2) =~ m/C[0-9]/) {
            ## This is an existing CLEVER gene
            ($version, $argName, $cleverID, $attr, $source, $accession) = split('\|',$arg);
            push(@argNames, $argName);
            push(@accessions, $accession);
        } else {
            $cleverID = "";
            if ($source eq "ResFinder") {
                #ResFinder-qnrB90_1_MG182074_1
                ($name,$variant,$accession,$version) = split('_', $rest);
                $argName = $name;
                $full_accession = $accession . "." . $version;
            }
            if ($source eq "CARD") {
                #CARD-gb|ADQ43424.1|ARO:3002746|QnrB31
                ($gb,$accession,$aro,$name) = split('\|', $rest);
                $argName = $name;
                $full_accession = $accession;
            }
            if ($source eq "ResFinderFG") {
                #ResFinderFG-beta_lactamase|KU545081.1|feces|ATM
                ($class,$accession,$aro,$origin,$antibiotic) = split('\|', $rest);
                $argName = "~" . $class;
                $full_accession = $accession;
            }
            if (($source eq "fARGene") || ($source eq "Inda-Diaz_2023")) {
                #fARGene-20k_concatenated-long-orfs_SFKQ01000019.1_seq1@@@methyltransferase_grp2_1
                ($seqinfo, $class)= split('@@@', $rest);
                ($junk1,$junk2,$accession,$variant) = split('_', $seqinfo);
                $argName = "~" . $class;
                $full_accession = $accession;
            }
        }
        $rep = $repMap{"$arg"};
        $newname = $argMap{"$arg"};
        $class = $classMap{"$arg"};
        $established = $estMap{"$arg"};
        if (defined($mobMap{"$arg"})) {
            $mobile = "M";
        } else {
            $mobile = "C";
        }
        $verified = $verMap{"$arg"};

        if ($rep eq "") {
            ## Rep seems to be unused in output, skip this
        } else {
            if ($newname eq "") {
                $newname = $argMap{"$rep"};
            }
            if ($class eq "") {
                $class = $classMap{"$rep"};
            }
            if ($established eq "") {
                $established = $estMap{"$rep"};
            }
            if ($verified eq "") {
                $verified = $verMap{"$rep"};
            }
        }


        if ($newname eq "") {
            $newname = $argName;
        }
        if ($class eq "") {
            ($class,$unused) = split('-', $cleverID);
            if (substr($class, 0, 1) eq "~"){
                $class = substr($class, 1);
            }
        }
        if ($established eq "") {
            $established = substr($attr, 0, 1);
        }
        if (substr($attr, 1, 1) eq "M") {
            $mobile = "M";
        }
        if ($verified eq "") {
            $verified = substr($attr, 2, 1);
        }

        if ($class eq "") {
            $class = "unc";
        }

        $familyNumber = $familyNumMap{"$arg"};
        $variantNumber = $varNumMap{"$arg"};

        if (substr($newname, 0, 7) eq "CLEVER:") {
            $newname = substr($newname, 7);
            $newaccession = "C" . $cleverversion . "|" . "~" . $class . "-" . $familyNumber . "_" . $variantNumber . "|" . $class . "-" . $familyNumber . "_" . $variantNumber . "|" . $established . $mobile . $verified . "|" . $source . "|" . $full_accession;
        } else {
            if ($cleverID ne "") {
                $newaccession = $arg;
                $newaccession =~ s/C[0-9]*\|//;
                $newaccession = "C" . $cleverversion . "|" . $newaccession;
            } else {
                $newaccession = "C" . $cleverversion . "|" . $newname . "|" . $class . "-" . $familyNumber . "_" . $variantNumber . "|" . $established . $mobile . $verified . "|" . $source . "|" . $full_accession;
            }
        }
        
        if ($newname eq "EXCLUDE") {
            $outputEntry = 0;
        } else {
            $outputEntry = 1;
            print ">" . $newaccession . " " . $class  . " " . $established . " " . $mobile . " " . $verified . " " . $source . " " . $arg . "\n";
        }
    } else {
        if ($outputEntry == 1) {
            print $line . "\n";
        }
    }
}
close FASTA;