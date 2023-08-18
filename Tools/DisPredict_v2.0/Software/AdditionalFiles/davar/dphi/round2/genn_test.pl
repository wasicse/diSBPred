#!/usr/bin/perl
#
#  A perl wrap for the NN fortran code
#  version: 1; build: 1 : original version
#  version: 1; build: 2 : added multi-outputs and output columns selection
#  version: 1; build: 3 : added option to scale inputs to a given interval
#			: test: test over all of db (database) file
# 
#  Eshel Faraggi 2009
if (($#ARGV >= 0))
  {
  print `head README`;
  die "\n\nSee README for more\nNo options necessary\n";
  }
else
  {
  $datestring = "echo \"Date: \"\`date\` >> nohup.out";
  `$datestring`;
  open(NOHUP,">>nohup.out");
  print NOHUP "Generating neural network ...\n";
  }
$mycompiler = $ENV{'fcompiler'};
$geninfl = "genn_test.in";
open (GENIN, "$geninfl") or die "Can't find file $geninfl\n";
@genin = <GENIN>;
close(GENIN);
$found = 0;
for (@genin)
  {
  @lv = split(/[\ \t\n]+/,$_);
  if ($lv[1] eq ":window:")  { $winsiz = $lv[0]; $found++; }
  if ($lv[1] eq ":level1:")  { $nhid1 = $lv[0]; $found++; }
  if ($lv[1] eq ":level2:")  { $nhid2 = $lv[0]; $found++; }
  if ($lv[1] eq ":dbfile:")  { $dbflnm = $lv[0]; $found++; }
  if ($lv[1] eq ":ninput:")  { $ninp = $lv[0]; $found++; }
  if ($lv[1] eq ":noutput:")  { $nout = $lv[0]; $found++; }
  if ($_ =~ ":outpos:")
    {
    $outpos = '';
    for ($i=0;$i<$nout-1;$i++)
         { $outpos .= "$lv[$i],"; }
    $outpos .= "$lv[$nout-1]";
    $found++;
    }
  }
if ($found<7)
  {
  die "Could not find all parameters ($found), (do no alter the white space structure of genn_test.in?)\nAborting\n";
  }
open(DBFL,"$dbflnm") or die "Can't find file $dbflnm\n";
@dbnm = <DBFL>;
close(DBFL);
$dbsize = @dbnm;
open(DBFL,"$dbnm[0]") or die "Can't find file $dbnm[0]\n";
@dbnm = <DBFL>;
close(DBFL);
@lv = split(/[\ \t\n]+/,$dbnm[0]);
$iflsz = @lv;
print NOHUP "Input: columns 1-$$ninp ; Outputs (prediction): column $iflsz\n";
print NOHUP "Window size: $winsiz ; Hidden: $nhid1, $nhid2\n";
print NOHUP "Database size: $dbsize\n";
`cp genn_test.in spiner.in`;
`cp genn_test.f spiner_test.f`;
#`replace ",nrx=21,nip=:filesize:,nw1=51,nw2=51,nrxo=1" ",nrx=$winsiz,nip=$ninp,nw1=$nhid1,nw2=$nhid2,nrxo=$nout" -- spiner_test.f`;
#`replace ":maxdbcol:)" "$iflsz)" -- spiner_test.f`;
#`replace ":listofoutputpos:" "$outpos" -- spiner_test.f`;
`cat spiner_test.f | sed "s/,nrx=21,nip=:filesize:,nw1=51,nw2=51,nrxo=1/,nrx=$winsiz,nip=$ninp,nw1=$nhid1,nw2=$nhid2,nrxo=$nout/g" > spiner_test1.f`;
`cat spiner_test1.f | sed "s/:maxdbcol:)/$iflsz)/g" > spiner_test2.f`;
`cat spiner_test2.f | sed "s/:listofoutputpos:/$outpos/g" > spiner_test3.f`;
`cp spiner_test3.f spiner_test.f`;
`rm spiner_test1.f spiner_test2.f spiner_test3.f`;
print NOHUP "Neural network initialized ...\n";
#$ifort = "/usr/local/intel/fc/9.1/bin/ifort";
#$ifort = "/data1/pub/intel/composerxe-2011.0.084/bin/intel64/ifort";
if ($mycompiler eq "ifort")
#if (-e "$ifort")
  {
  #print `$ifort -132 -shared-intel -mcmodel=medium -O2 -o spiner_test.e spiner_test.f`; 
  print `ifort -132 -shared-intel -mcmodel=medium -O2 -o spiner_test.e spiner_test.f`; 
  $complr = "ifort";
  }
else
  {
  #print NOHUP "Can't find $ifort\n";
  #print NOHUP "Trying gfortran compiler ...\n";
  $tt = system("gfortran -ffixed-line-length-none -O2 -o spiner_test.e spiner_test.f");
  if ($tt)
    {
    die "\n\nCan't compile with either gfortran or ifort, aborting\n\n";      
    }
  $complr = "gfortran";
  }
print NOHUP "Compiled with $complr ...\n";
print NOHUP "Predicting test set ...\n";
close(NOHUP);
print `cat nohup.out`;
print "\nOutput files spiner.valdav and nohup.out will be updated intermittently\n\n";
print "\Final prediction files are: Test-set: pred-testfold.out\n\n";
#exec('nohup time ./spiner_test.e &');
exec('./spiner_test.e >>nohup.out');
__END__
