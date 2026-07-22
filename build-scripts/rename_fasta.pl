#!/usr/bin/perl

$fastafile = shift;
$annotationfile = shift;
$mappingfile = shift;
$mobilefile = shift;
$blacklistfile = shift;
$os_classfile = shift;
$cleverversion = shift;
$variantmap = shift;

%repMap = {};
%argMap = {};
%classMap = {};
%familyNum = {};
%maxFamilyNum = {};
%varNum = {};
%maxVarNum = {};
%familyNumMap = {};
%varNumMap = {};
%estMap = {};
%mobMap = {};
%verMap = {};
%lineageMap = {};
%familyMap = {};
%children = {};
%grandchildren = {};
%variantMap = {};

$level = "unknown";
if ($fastafile =~ m/variants/) {
    $level = "variant";
}
if ($fastafile =~ m/families/) {
    $level = "family";
}
if ($fastafile =~ m/lineages/) {
    $level = "lineage";
}

open (ANNOT, $annotationfile);
while ($line = <ANNOT>) {
    chomp($line);
    ($clusterRep, $autoARG, $curatedARG, $autoclass, $class, $lineage, $est, $ver, $clusterentries) = split('\t', $line);
    (@entries) = split(',', $clusterentries);
    
    if ($clusterRep =~ m/^C[0-9]/) {
        ($version, $geneName, $geneID) = split('\|', $clusterRep);
        ($classB, $family) = split("-", $geneID);
        ($famNumber, $variantNumber) = split("_", $family);
        if (defined($maxFamilyNum{$classB})) {
            if ($famNumber > $maxFamilyNum{$classB}) {
                $maxFamilyNum{$classB} = $famNumber;
            }
        } else {
            $maxFamilyNum{$classB} = $famNumber;
        }
        if (defined($maxVarNum{"$classB-$famNumber"})) {
            if ($variantNumber > $maxVarNum{"$classB-$famNumber"}) {
                $maxVarNum{"$classB-$famNumber"} = $variantNumber;
            }
        } else {
            $maxVarNum{"$classB-$famNumber"} = $variantNumber;
        }
        $classMap{$clusterRep} = $classB;
    }
    
    foreach $entry (@entries) {
        if ($entry =~ m/^C[0-9]/) {
            ($version, $geneName, $geneID) = split('\|', $entry);
            ($classB, $family) = split("-", $geneID);
            ($famNumber, $variantNumber) = split("_", $family);
            
            if (defined($maxVarNum{"$classB-$famNumber"})) {
                if ($variantNumber > $maxVarNum{"$classB-$famNumber"}) {
                    $maxVarNum{"$classB-$famNumber"} = $variantNumber;
                }
            } else {
                $maxVarNum{"$classB-$famNumber"} = $variantNumber;
            }
            $familyNumMap{$entry} = $famNumber;
            $varNumMap{$entry} = $variantNumber;
        } 

        $repMap{$entry} = $clusterRep;
        $argMap{$entry} = $curatedARG;
        $classMap{$entry} = $class;
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

#open (LINEAGE, $lineagefile);
#while ($line = <LINEAGE>) {
#    chomp($line);
#    ($gene, $lineage) = split('\t', $line);
#    ($version, $genename, $cleverID) = split('\|', $lineage);
#    $class = $cleverID;
#    $class =~ s/-.*//;
#    $lineageMap{$gene} = $class;
#}
#close LINEAGE;

open (MAPPING, $mappingfile);
while ($line = <MAPPING>) {
    chomp($line);
    ($variantID, $genename, $familyID, $lineageID, $class, $blacklisted) = split('\t', $line);
    ($linversion, $lingenename, $lincleverID) = split('\|', $lineageID);
    ($famversion, $famgenename, $famcleverID) = split('\|', $familyID);

    #if ($class ne $classMap{$variantID}) {
    #    if ($classMap{$variantID} ne "") {
    #        print STDERR "ERROR: Class of mapping and annotation sequence doesn't match for $variantID: $class vs. $classMap{$variantID}\n";
    #    }
    #}
    #$classMap{$variantID} = $class;
    $lineageMap{$variantID} = $lineageID;
    if (defined($children{$familyID})) {
        $children{$familyID} = $children{$familyID} . "\n" . $variantID;
    } else {
        $children{$familyID} = $variantID;
    }
    if (defined($grandchildren{$lineageID})) {
        $grandchildren{$lineageID} = $grandchildren{$lineageID} . "\n" . $variantID;
    } else {
        $grandchildren{$lineageID} = $variantID;
    }
    $familyMap{$variantID} = $familyID;
}
close MAPPING;

open (OSCLASSES, $os_classfile);
while ($line = <OSCLASSES>) {
    chomp($line);
    ($os_type, $clever_class) = split('\t', $line);
    $osMap{$os_type} = $clever_class;
    push(@os_terms, $os_type);
}
close OSCLASSES;

if (($level ne "variant") && ($variantmap ne "")) {
    open (VARIANTMAP, $variantmap);
    while ($line = <VARIANTMAP>) {
        chomp($line);
        ($newID, $oldID, $newGeneID) = split('\t', $line);
        $variantMap{$newID} = $newGeneID;
        $variantMap{$oldID} = $newGeneID;
    }
    close VARIANTMAP;
}

$outputEntry = 0;
if ($level eq "variant") {
    open (VARIANTMAP, ">$fastafile.variantmap.txt");
}
open (FASTA, $fastafile);
while ($line = <FASTA>) {
    chomp($line);
    if (substr($line, 0 ,1) eq ">") {
        $blacklisted = 0;
        foreach $item (@blacklist) {
            if ($item ne "") {
                if (($line =~ m/^$item[^A-Za-z]/) || ($line =~ m/[^A-Za-z]$item[^A-Za-z]/) || ($line =~ m/[^A-Za-z]$item$/)) {
                    $blacklisted = 1;
                }
            }
        }
        if ($blacklisted == 1) {
            $outputEntry = 0;
            next;
        }
        $attr = "";
        $class = "";
        $familyNumber = "";
        $variantNumber = "";
        $established = "";
        $mobile = "";
        $verified = "";
        $full_accession = "";
        $accession = "";
        $argName = "";
        $cleverID = "";
        $newaccession = "";
        $lineageID = "";
        $rep = "";
        $newname = "";

        ($arg) = split(' ', $line);
        $arg = substr($arg, 1);
        ($source, $rest) = split("-",$arg);
        $rest = $arg;
        $rest =~ s/^[^-]*-//;
        $established = "";
        $verified = "";
        $do_not_change_cleverID = 0;
        if (defined($variantMap{$arg})) {
            if ($variantMap{$arg} ne "") {
                $cleverID = $variantMap{$arg};
                $do_not_change_cleverID = 1;
                ($family, $variantNumber) = split("_", $cleverID);
                ($gclass, $familyNumber) = split("-", $family);
                $varNumMap{$arg} = $variantNumber;
                $familyNumMap{"$arg"} = $familyNumber;
            }
        }
        if (substr($arg, 0, 2) =~ m/C[0-9]/) {
            ## This is an existing CLEVER gene
            ($version, $argName, $cleverID, $attr, $source, $accession) = split('\|',$arg);
            ($gclass_fam, $variantNo) = split("_", $cleverID);
            if ($do_not_change_cleverID == 0) {
                $varNumMap{$arg} = $variantNo;
            }
            $established = substr($attr, 0, 1);
            $verified = substr($attr, 2, 1);
        } else {
            if ($do_not_change_cleverID == 0) {
                $cleverID = "";
            }
            if ($source eq "ResFinder") {
                #ResFinder-qnrB90_1_MG182074_1
                ($name,$variant,$accession,$version) = split('_', $rest);
                $argName = $name;
                $full_accession = $accession . "." . $version;
                $verified = "V";
            }
            if ($source eq "CARD") {
                #CARD-gb|ADQ43424.1|ARO:3002746|QnrB31
                ($gb,$accession,$aro,$name) = split('\|', $rest);
                $argName = $name;
                $full_accession = $accession;
                $verified = "V";
            }
            if ($source eq "ResFinderFG") {
                #ResFinderFG-beta_lactamase|KU545081.1|feces|ATM
                ($class,$accession,$aro,$origin,$antibiotic) = split('\|', $rest);
                $argName = "!" . $class;
                $full_accession = $accession;
                $verified = "V";
            }
            if (($source eq "fARGene") || ($source eq "Inda-Diaz_2023") || ($source eq "IndaDiaz_2023")) {
                #fARGene-20k_concatenated-long-orfs_SFKQ01000019.1_seq1@@@methyltransferase_grp2_1
                ($seqinfo, $classish)= split('@@@', $rest);
                foreach $term (@os_terms) {
                    if ($classish =~ /$term/) {
                        $class = $osMap{$term};
                        last;
                    }
                }
                ($junk1,$junk2,$accession,$variant) = split('_', $seqinfo);
                $argName = "!" . $class;
                $full_accession = $accession;
                $verified = "P";
            }
            if ($source eq "Victor_2025") {
                #HiAAMG_A_sequ1ence
                ($seqinfo, $junk) = split('_seq', $rest);
                foreach $term (@os_terms) {
                    if ($rest =~ /$term/) {
                        $class = $osMap{$term};
                        last;
                    }
                }
                $accession = $rest;
                $argName = "!" . $class;
                $full_accession = $accession;
                $verified = "P";
            }
            if ($source eq "Victor_2025b") {
                #Class_ANSF1A
                ($seqinfo, $junk) = split('NSF', $rest);
                $accession = $rest;
                foreach $term (@os_terms) {
                    if ($rest =~ /$term/) {
                        $class = $osMap{$term};
                        last;
                    }
                }
                $argName = "!" . $class;
                $full_accession = $accession;
                $verified = "P";
            }
            if ($source eq "Li_2025") {
                #class_a|cattle|k141_371442_seq1_1
                ($classish, $animal, $accession) = split('\|', $rest);
                $classish = $classish . "|";
                foreach $term (@os_terms) {
                    if ($classish =~ /$term/) {
                        $class = $osMap{$term};
                        last;
                    }
                }
                $argName = "!" . $class;
                $full_accession = $accession;
                $verified = "P";
            }
            if ($source eq "Sommerville_2026") {
                #>retrieved-translated_gene_catalog_Gene026555_seq1_1|class_A
                ($accession,$classish) = split('\|', $rest);
                foreach $term (@os_terms) {
                    if ($classish =~ /$term/) {
                        $class = $osMap{$term};
                        last;
                    }
                }
                $argName = "!" . $class;
                $full_accession = $accession;
                $verified = "P";
            }
            if (($source eq "Mustard") || ($source eq "Ruppe_2019")) {
                #1|311066|3|Human-Microbiome-3.9M-Gene-Catalog-(metahit-v2)|MC3.MG12.AS1.GP1.C65190.G5|aac2
                ($x, $y, $z, $refcatalog, $accession, $classish)= split('\|', $rest);
                foreach $term (@os_terms) {
                    if ($classish =~ /$term/) {
                        $class = $osMap{$term};
                        last;
                    }
                }
                $argName = "!" . $class;
                $full_accession = $accession;
                $verified = "S";
            }
            if ($source eq "Wang_2025") {
                #DfrA52_XPO54507.1
                ($geneName, $accession)= split('_', $rest);
                foreach $term (@os_terms) {
                    if ($rest =~ /$term/) {
                        $class = $osMap{$term};
                        last;
                    }
                }
                #$class = lc(substr($geneName, 0 , 3));
                $argName = "!" . $geneName;
                $full_accession = $accession;
                $verified = "S";
            }
        }
        $rep = $repMap{"$arg"};
        $newname = $argMap{"$arg"};
        $family = $familyMap{"$arg"};
        $lineage = $lineageMap{"$arg"};
        if (defined($classMap{"$arg"})) {
            if ($classMap{"$arg"} ne "") {
                $class = $classMap{"$arg"};
            }
        } else {

        }
        #$established = $estMap{"$arg"};
        if (defined($mobMap{"$arg"})) {
            $mobile = "M";
        } else {
            $mobile = "C";
        }
        if (defined($verMap{"$arg"})) {
            $verified = $verMap{"$arg"};
        } else {
            if ($verified eq "") {
                $verified = "P";
            }
        }

        if ($rep eq "") {
            ## Rep seems to be unused in output, skip this
        } else {
            if ($newname eq "") {
                $newname = $argMap{"$rep"};
            }
            
            if ($class eq "") {
                $class = $classMap{"$rep"};
                print STDERR "No class found for $arg Chosing class of rep ($rep): $classMap{$rep}\n";
            }
            if ($established eq "") {
                #$established = $estMap{"$rep"};
                $established = "L";
            }
            if ($verified eq "") {
                #if (defined($verMap{"$rep"})) {
                #    $verified = $verMap{"$rep"};
                #}
                $verified = "P";
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
            print STDERR "No class found for $arg Forced to base class on CLEVER ID: $class\n";
        }
        if ($class eq "") {
            if (defined($lineageMap{"$arg"})) {
                $lineageID = $lineageMap{"$arg"};
                ($linversion, $lingenename, $linclass) = split('\|', $lineageID);
                ($class, $no) = split('-', $linclass);
            }
            print STDERR "No class found for $arg Deriving class from lineage ($lineageID): $class\n";
        }
        if ($established eq "") {
            if ($attr ne "") {
                $established = substr($attr, 0, 1);
            } else {
                $established = "L";
            }
        }
        if ($attr ne "") {
            if (substr($attr, 1, 1) eq "M") {
                $mobile = "M";
            } else {
                $mobile = "C";
            }
        }
        if ($verified eq "") {
            if ($attr ne "") {
                $verified = substr($attr, 2, 1);
            } else {
                if (($source eq "Ruppe_2019") || ($source eq "Wang_2025")) {
                    $verified = "S";
                } else {
                    $verified = "P";
                }
            }
        }

        ## Look for children/grandchildren attributes
        if ($level eq "family") {
            $pest = substr($attr, 0, 1);
            $pmob = substr($attr, 1, 1);
            $pver = substr($attr, 2, 1);
            (@child_ids) = split('\n', $children{$arg});
            foreach $child (@child_ids) {
                if ($child =~ m/^C[0-9]/) {
                    ($cversion, $cargName, $ccleverID, $cattr, $csource, $caccession) = split('\|',$child);
                    $cest = substr($cattr, 0, 1);
                    $cmob = substr($cattr, 1, 1);
                    $cver = substr($cattr, 2, 1);
                    if ($cest eq "E") {
                        $established = "E";
                    }
                    if ($cmob eq "M") {
                        $mobile = "M";
                    }
                    if ($cver eq "V") {
                        $verified = "V";
                    }
                }
            }
        }

        if ($level eq "lineage") {
            $pest = substr($attr, 0, 1);
            $pmob = substr($attr, 1, 1);
            $pver = substr($attr, 2, 1);
            (@child_ids) = split('\n', $grandchildren{$arg});
            foreach $child (@child_ids) {
                if ($child =~ m/^C[0-9]/) {
                    ($cversion, $cargName, $ccleverID, $cattr, $csource, $caccession) = split('\|',$child);
                    $cest = substr($cattr, 0, 1);
                    $cmob = substr($cattr, 1, 1);
                    $cver = substr($cattr, 2, 1);
                    if ($cest eq "E") {
                        $established = "E";
                    }
                    if ($cmob eq "M") {
                        $mobile = "M";
                    }
                    if ($cver eq "V") {
                        $verified = "V";
                    }
                }
            }
        }

        if ($class eq "") {
            $class = "new";
            print STDERR "No class found for $arg Flagged as 'new' class\n";
        }

        if ($do_not_change_cleverID == 0) {
            $familyNumber = $familyNumMap{"$arg"};
            $variantNumber = $varNumMap{"$arg"};

            $correctFamily = "";
            if (($familyNumber eq "") && ($family ne "")) {
                if (defined($correctedFamilies{$family})) {
                    $correctFamily = $family;
                    $family = $correctedFamilies{$family};
                }
                if ($family =~ m/[A-Za-z]*-[0-9]*_[0-9]*/) {
                    ($xclassname, $famandclass) = split('-', $family);
                    ($familyNumber, $variantNumber) = split('_', $famandclass);
                    #$variantNumber =~ s/[^0-9].*//;
                    #$varNumMap{"$arg"} = $variantNumber;
                    $variantNumber = $maxVarNum{"$class-$familyNumber"} + 1;
                    $maxVarNum{"$class-$familyNumber"} = $variantNumber;
                    $varNumMap{$arg} = $maxVarNum{"$class-$familyNumber"};
                    $familyNumMap{"$arg"} = $familyNumber;
                    
                }
            }
            if (($familyNumber eq "") && ($newname ne "")) {
                if ($newname =~ m/[A-Za-z]*-[0-9]*_[0-9]*/) {
                    ($xclassname, $famandclass) = split('-', $newname);
                    ($familyNumber, $butnotvariantNumber) = split('_', $famandclass);
                    $familyNumMap{"$arg"} = $familyNumber;
                }
            }

            $addFam = 0;
            if ($familyNumber eq "") {
                $familyNumber = $maxFamilyNum{$class} + 1;
                $maxFamilyNum{$class} = $maxFamilyNum{$class} + 1;
                $maxVarNum{"$class-$familyNumber"} = 0;
                $addFam = 1;
                #print STDERR "New family: $class-$familyNumber Rep: $rep\n";
            }
            if ($variantNumber eq "") {
                $variantNumber = $maxVarNum{"$class-$familyNumber"} + 1;
                $maxVarNum{"$class-$familyNumber"} = $variantNumber;
                $varNumMap{$arg} = $maxVarNum{"$class-$familyNumber"};
                #print STDERR "$arg\t$class-$familyNumber\t$varNumMap{$arg}\t$variantNumber\n";
            }
            if ($correctFamily ne "") {
                $correctedFamilies{$correctFamily} = $class . "-" . $familyNumber . "_" . $variantNumber;
            }
        }

        if ($do_not_change_cleverID == 0) {
            if (substr($newname, 0, 1) eq "!") {
                $newname = substr($newname, 1);
                $newaccession = "C" . $cleverversion . "|" . "~" . $class . "-" . $familyNumber . "_" . $variantNumber . "|" . $class . "-" . $familyNumber . "_" . $variantNumber . "|" . $established . $mobile . $verified . "|" . $source . "|" . $full_accession;
            } else {
                if ($cleverID ne "") {
                    $newaccession = $arg;
                    ($cversion, $cargName, $ccleverID, $cattr, $csource, $caccession) = split('\|',$arg);
                    $newaccession = "C" . $cleverversion . "|" . $cargName . "|" . $ccleverID. "|" . $established . $mobile . $verified . "|" . $csource  . "|" . $caccession;
                } else {
                    if ($newname =~ m/[A-Za-z]*-[0-9]*_[0-9]*/) {
                        $newaccession = "C" . $cleverversion . "|" . "~" . $newname . "|" . $class . "-" . $familyNumber . "_" . $variantNumber . "|" . $established . $mobile . $verified . "|" . $source . "|" . $full_accession;
                    } else {
                        $newaccession = "C" . $cleverversion . "|" . $newname . "|" . $class . "-" . $familyNumber . "_" . $variantNumber . "|" . $established . $mobile . $verified . "|" . $source . "|" . $full_accession; 
                    }
                    
                }
            }
        } else {
            if (substr($newname, 0, 1) eq "!") {
                $newname = substr($newname, 1);
                $newaccession = "C" . $cleverversion . "|" . "~" . $class . "-" . $familyNumber . "_" . $variantNumber . "|" . $cleverID . "|" . $established . $mobile . $verified . "|" . $source . "|" . $full_accession;
            } else {
                if ($cleverID ne "") {
                    $newaccession = $arg;
                    ($cversion, $cargName, $ccleverID, $cattr, $csource, $caccession) = split('\|',$arg);
                    $newaccession = "C" . $cleverversion . "|" . $cargName . "|" . $cleverID . "|" . $established . $mobile . $verified . "|" . $csource  . "|" . $caccession;
                } else {
                    if ($newname =~ m/[A-Za-z]*-[0-9]*_[0-9]*/) {
                        $newaccession = "C" . $cleverversion . "|" . "~" . $newname . "|" . $cleverID . "|" . $established . $mobile . $verified . "|" . $source . "|" . $full_accession;
                    } else {
                        $newaccession = "C" . $cleverversion . "|" . $newname . "|" . $cleverID . "|" . $established . $mobile . $verified . "|" . $source . "|" . $full_accession;
                    }
                }
            }
        }

        if ($addFam == 1) {
            if ($do_not_change_cleverID == 0) {
                #print STDERR "$arg\t$newaccession\t$class-$familyNumber\_$variantNumber\t$class\t$familyNumber\t" . $varNum{"$class-$familyNumber"} . "\n";
                $repMap{$arg} = $newaccession;
                $argMap{$arg} = $class . "-" . $familyNumber . "_" . $variantNumber;
                $classMap{$arg} = $class;
                $familyNumMap{"$arg"} = $familyNumber;
                $varNumMap{"$arg"} = $variantNumber;
            } else {
                $repMap{$arg} = $newaccession;
                $argMap{$arg} = $cleverID;
                $classMap{$arg} = $class;
                $familyNumMap{"$arg"} = $familyNumber;
                $varNumMap{"$arg"} = $variantNumber;
            }
        }
        
        if ($newname eq "EXCLUDE") {
            $outputEntry = 0;
        } else {
            $outputEntry = 1;
            if ($level eq "variant") {
                print VARIANTMAP $arg . "\t" . $newaccession . "\t" . $class . "-" . $familyNumber . "_" . $variantNumber . "\n";
            }
            print ">" . $newaccession . " " . $class  . " " . $established . " " . $mobile . " " . $verified . " " . $source . " " . $arg . "\n";
        }
    } else {
        if ($outputEntry == 1) {
            print $line . "\n";
        }
    }
}
if ($level eq "variant") {
    close VARIANTMAP;
}
close FASTA;