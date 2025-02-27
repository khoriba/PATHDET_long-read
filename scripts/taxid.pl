#!/usr/bin/env perl
use strict;
use warnings;
use FileHandle;
use DirHandle;
use File::Basename;

my $dir0=shift or die;
my $nodes=shift or die;
my @grps=qw/bacteria fungi viral/;

my $phr=&phr($nodes);
my $thr={};

if (1) {
  for my $grp (@grps) {
    my $dir=$dir0.'/'.$grp;
    -d $dir or next;
    my $dh=DirHandle->new($dir) or die;
    my @ac=grep {$_ !~ /^\.\.?$/} $dh->read;
    for my $ac(@ac) {
      my $acdir=$dir.'/'.$ac;
      my $name = &name($acdir);
      my $rep=$acdir.'/'.$ac.'_'.$name.'_assembly_report.txt';
      my $fa= $acdir.'/'.$ac.'_'.$name.'_genomic.fna.gz';
      -f $rep or die $rep;
      my $fh=FileHandle->new($rep) or die;
      while (<$fh>){
        s/\r\n$//;
        /^# Taxid:/ or next;
        /^# Taxid: +(\d+)$/ or die $_;
	my $id=$1;
	my $lin=&lineage($id);
	$thr->{$id}->{$ac}={lineage=>$lin, group=>$grp, fasta=>$fa};
      }
    }
  }
}
if (1) {
  my @buf=();
  for my $tid (sort {$a<=>$b} keys %$thr) {
    for my $ac ( keys %{$thr->{$tid}}) {
      push @buf, join "\t", $thr->{$tid}->{$ac}->{lineage}, $tid, $thr->{$tid}->{$ac}->{fasta};
    }
  }
  printf "%s\n", join "\n", @buf;
}

sub name {
  my $acdir=shift or die;
  my $ac=basename($acdir);
  my $dh=DirHandle->new($acdir) or die;
  my @rep=grep {$_ =~ /_assembly_report\.txt$/} $dh->read;
  $#rep == 0 or die $acdir;
  my $rep=$rep[0];
  $rep=~/${ac}_(.+)_assembly_report\.txt$/ or die $rep;
  $1;
}

sub phr{
  my $nodes=shift or die;
  my $hr={};
  my $fh=FileHandle->new($nodes) or die;
  for (<$fh>) {
    chomp;
    my @col = split /\t/, $_, -1;
    $col[0]=~/^\d+$/ or die;
    $col[1]=~/^|$/ or die;
    $col[2]=~/^\d+$/ or die;
    $hr->{$col[0]}=$col[2];
  }
  $hr;
}

sub lineage {
  my $id=shift or die;
  my @buf=($id);
  while(exists $phr->{$id}) {
    $id=$phr->{$id};
    push @buf, $id;
    last if $id == 1;
  }
  join ',', @buf;;
}
