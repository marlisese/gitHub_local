#!/usr/bin/perl -w
use File::Copy;
$path=shift; # location to create files
$pathBayes=shift;

for($i=1; $i<=100; $i++) {
        $fileMatrix="$path/$i-GFsim.txt";
        print  "DEBUG: $fileMatrix\n";
        system("$pathBayes/BayesTraits $path/tree_totest.nex $fileMatrix < $path/command_multiple.txt>>$path/results.txt");
        }   
exit 0;


# this script sends all simulations to BayesTraits to be fit under independent or dependent model
# uses one tree, but could use multiple if nexus file contains more 
# usage : ./runBayesTrait.pl path/to/data /path/to/program/
# specify names of simulations file , tree file and # of simulations inside
# command_multiple.txt is a file with the numbers of the analyses in BayesRates, i.e. 2 1 Run
