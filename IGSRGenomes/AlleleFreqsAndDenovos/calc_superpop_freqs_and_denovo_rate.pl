#!/usr/bin/perl -w
#

use strict;

my $Usage = "Usage: calc_superpop_freqs_and_denovo_rate.pl <genotypes>\n";

$#ARGV==0
    or die $Usage;

my $geno_file = $ARGV[0];

my $samples_pops = '1000G_2504_high_coverage.sequence.index.txt';
my $samples_pops_related = '1000G_698_related_high_coverage.sequence.index.txt';
my $superpops = 'superpopulations.txt';
my $pedigree_file = '1kGP.3202_samples.pedigree_info.txt';

my $rh_sample_genos = read_sample_genos($geno_file);
my $rh_sample_pops = read_sample_pops($samples_pops, $samples_pops_related);
my $rh_superpops = read_superpops($superpops);
my $rh_parents = read_sample_parents($pedigree_file);

my %pop_samples = ();
my %superpop_samples = ();
foreach my $sample (keys %{$rh_sample_pops}) {
    my $pop = $rh_sample_pops->{$sample};
    my $superpop = $rh_superpops->{$pop};
    $pop_samples{$pop}++;
    $superpop_samples{$superpop}++;
}

# tally alleles by superpopulation:

my %total_alleles = ();
my %allele_counts = ();
my $total_child_alleles = 0;
my $no_children = 0;
my $no_consistent_children = 0;
my $no_denovo_children = 0;
my $total_alleles_seen = 0;
my %total_ac = ();
foreach my $sample (keys %{$rh_sample_genos}) {
    my $geno = $rh_sample_genos->{$sample}->{'geno'};
    my $ra_alleles = $rh_sample_genos->{$sample}->{'alleles'};
    my $samplepop = $rh_sample_pops->{$sample};
    if (!$samplepop) {
        die "No sample population for sample $sample!\n";
    }
    my $superpop = $rh_superpops->{$samplepop};
    $allele_counts{$superpop} = {} if (!$allele_counts{$superpop});

    if (!($rh_parents->{$sample}->{'father'}) && !($rh_parents->{$sample}->{'mother'})) {
        $total_alleles_seen += 2;
        $total_alleles{$superpop} = 0 if (!$total_alleles{$superpop});
        $total_alleles{$superpop} += 2;
        foreach my $allele (@{$ra_alleles}) {
            next if $allele eq 'WTYP';
            $allele_counts{$superpop}->{$allele} = 0 if (!$allele_counts{$superpop}->{$allele});
            $allele_counts{$superpop}->{$allele}++;
            $total_ac{$allele}++;
        }
    }

    if ($rh_parents->{$sample}) { # check if we have both parents
        my $father = $rh_parents->{$sample}->{'father'};
        my $mother = $rh_parents->{$sample}->{'mother'};
        if ($father && $mother && $rh_sample_genos->{$father} && $rh_sample_genos->{$mother}) { # have info for both parents
            my $father_geno = $rh_sample_genos->{$father}->{'geno'};
            my $mother_geno = $rh_sample_genos->{$mother}->{'geno'};
            my $ra_father_alleles = $rh_sample_genos->{$father}->{'alleles'};
            my $ra_mother_alleles = $rh_sample_genos->{$mother}->{'alleles'};
            $no_children++; 
            if ((grep {$_ eq $ra_alleles->[0]} @{$ra_father_alleles}) &&
                 (grep {$_ eq $ra_alleles->[1]} @{$ra_mother_alleles})) {
                #print "CONSISTENTCASE1\t$geno\t$father_geno\t$mother_geno\n";
                $no_consistent_children++;
            }
            elsif ((grep {$_ eq $ra_alleles->[0]} @{$ra_mother_alleles}) &&
                 (grep {$_ eq $ra_alleles->[1]} @{$ra_father_alleles})) {
                #print "CONSISTENTCASE2\t$geno\t$father_geno\t$mother_geno\n";
                $no_consistent_children++;
            }
            else {
                print "DENOVOCASE\t$geno\t$father_geno\t$mother_geno\t$superpop\n";
                $no_denovo_children++;
            }
        }
    }
}

print "$no_children children had two parents with genotypes: $no_consistent_children had consistent inheritance and $no_denovo_children had denovo mutation\n";

# report allele freqs for each superpopulation:

foreach my $superpop (keys %superpop_samples) {
    print "$superpop has $superpop_samples{$superpop} samples\n";
    my $no_alleles = $total_alleles{$superpop} || 0;
    foreach my $allele (sort keys %{$allele_counts{$superpop}}) {
        my $ac = $allele_counts{$superpop}->{$allele} || 0;
        my $af = $ac/$no_alleles;
        print "$superpop\t$allele\t$ac\t$no_alleles\t$af\n";
    }
}
foreach my $allele (sort keys %total_ac) {
    my $af = $total_ac{$allele}/$total_alleles_seen;
    print "Overall\t$allele\t$total_ac{$allele}\t$total_alleles_seen\t$af\n";
}


sub read_sample_genos {
    my $file = shift;

    open GENOS, $file
        or die "Couldn\'t open $file: $!\n";

    my %sample_genos = ();    
    while (<GENOS>) {
        chomp;
        next if (/^#/);
   
        my @fields = split /\s/, $_;
        my @alleles = split /_/, $fields[1];
        $sample_genos{$fields[0]} = {'geno' => $fields[1], 'alleles' => [@alleles]};
    }
    close GENOS;

    return {%sample_genos};
}

sub read_sample_pops {
    my @files = @_;

    my %samplepops = ();
    foreach my $file (@files) {
        open SAMPLES, $file
            or die "Couldn\'t open $file: $!\n";
    
        while (<SAMPLES>) {
            chomp;
            next if (/^#/);
    
            my @fields = split /\t/, $_;
            $fields[9] =~ s/\s+$//;
            $fields[10] =~ s/\s+$//;
            $samplepops{$fields[9]} = $fields[10];
        }
        close SAMPLES;
    }

    return {%samplepops};
}

sub read_superpops {
    my $file = shift;

    open SUPERS, $file
        or die "Couldn\'t open $file: $!\n";

    my %superpops = ();
    while (<SUPERS>) {
        chomp;
        next if (/^Population/);

        my @fields = split /\t/, $_;

        $superpops{$fields[0]} = $fields[6];
    }
    close SUPERS;

    return {%superpops};
}

sub read_sample_parents {
    my $file = shift;

    open PED, $file
        or die "Couldn\'t open $file: $!\n";

    my %sample_parents = ();    
    while (<PED>) {
        chomp;
        next if (/^sampleID/i);
   
        my ($sample, $fatherid, $motherid, $sex) = split /\s/, $_;
        if ($fatherid) {
            $sample_parents{$sample}->{'father'} = $fatherid;
        }
        if ($motherid) {
            $sample_parents{$sample}->{'mother'} = $motherid;
        }
    }
    close PED;
    return {%sample_parents};
}

