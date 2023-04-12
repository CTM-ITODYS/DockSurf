#!/usr/bin/perl
# *****************
# *  DockSurf.pl  *
# *****************
#
# usage: DockSurf.pl label pdbfile fine surfacetype sel clusconf
#
# Where:
# label is a nickname used for the generation of all output names
# pdbfile is a protein data bank file of a protein.
# fine is an integer of angle increment typical value is 10
# surfacetype is a name of surface type. For example Au111
# sel is an integer for the number of conformation/cluster to be selected. Typical value is 5
# cluscon is just a tag (clus or conf) to choose either $sel clusters or conformations
#
$| = 1;
use strict;
use Math::Trig;
#-----------------------------------------------------------------------------------------------------------------------------
# Get the input arguments and declarations
#
my $label=$ARGV[0];
my $pdb_input=$ARGV[1];
my $fine=$ARGV[2];
my $surf_type=$ARGV[3];
my $display=$ARGV[4];
my $clusconf=$ARGV[5];
my @type=();
my $i=0;
my $ind=0;
my $itype=0;
my $col_x=0;
my $col_y=0;
my $col_z=0;
my $x_ori=0;
my $y_ori=0;
my $z_ori=0;
my @tag=();
my @xb=();  #xb begin, xt after translation and x table for rotation
my @yb=();
my @zb=();
my @xt=();
my @yt=();
my @zt=();
my @x=();
my @y=();
my @z=();
my $XN=0;
my $YN=0;
my $ZN=0;
my @t=();
my @index=();
my @atom_type=();
my @res_type=();
my @res_numb=();
my $end=0;
my $endS=0;
my $prec=0;
my $nuta=0;
my $score_FF=0;
my $score_QM=0;
my $score_FFh=0;
my $score_QMh=0;
my $score_hydro=0;
my $zmin_value=0;
my $xmin_value=0;
my $ymin_value=0;
my $xMax_value=0;
my $yMax_value=0;
my $zmin_test=0;
my $zmin_testH=0;
my $dist=0;
my $rank=0;
my $res=0;
my $temp_FF=0;
my $temp_QM=0;
my $temp_hydro=0;
my $thresh=0;   # specific value for scoring function
my $threshH=0;  # hydrophobicity acts only at the surface contact
my $row_grid=0;
my $col_grid=0;
my $image="";
#my $image_Hmm="";
#my $image_Hqm="";
#my $image_Hsolv="";
#my $image_Gmm="";
#my $image_Gqm="";
my @tBA = ();
my @tCM = ();
my @tGA = ();
my @tCS = ();
my @tIM = ();
my @tAM = ();
my @tBE = ();
my @tME = ();
my @tCH = ();
my @prece=();
my @nutat=();
my @scoreHff=();
my @scoreHqm=();
my @scoreSolv=();
my @scoreGff=();
my @scoreGqm=();
my @zmin=();
my @xmin=();
my @ymin=();
my @xMax=();
my @yMax=();
my @table=();
my @indices=();
my @indicestries=();
my @name_scoring=();
my $iclus=0;
my $distP=0;
my $distN=0;
my $num_select=0;
my $tmp=0;
my $cluster_oui_non="";
#-----------------------------------------------------------------------------------------------------------------------------
# Declarations for the output data
# description
# A map file compiling the euler's angles, various scorings and geometrical parameters.
# A script for ploting 3D maps and generate various images
# A name for the generated pdb
my $mapfile="";
my $outpdb=""; #for the last part
my $scriptplotfile="";
my $score_Hff="";
my $score_Hqm="";
my $score_Solv="";
my $score_Gff="";
my $score_Gqm="";
my $espace="";
my $resid=0;
my $Ntag1="";
my $Ntag2="";
my $Ntag3="";
my $Ntag4="";
my $j=0;
my $mZ=0;
my $tag1="";
my $tag2="";
my $tag3="";
my $tag4="";
my $nb_Xunit=0;
my $nb_Yunit=0;
my $DistX=0;
my $DistY=0;
my $mX=0;
my $MX=0;
my $mY=0;
my $MY=0;
my $mZ=0;
my $dX=0;
my $dY=0;
my $cX=0;
my $cY=0;
my $ZminGlob=0;
my $Num_Unit_X=0;
my $Num_Unit_Y=0;
my $XNo=0;
my $YNo=0;
my $pdb_surface="";
my $temp_name="";
my @title_gnuplot=();
my $col=0;
#-----------------------------------------------------------------------------------------------------------------------------
# reading surface type, this part may be improved by adding more surfaces
#
if ($surf_type eq "Au111" )
{
  @type=('C','CA','CB','CD','CD1','CD2','CE','CE1','CE2','CE3','CG','CG1','CG2','CH2','CZ','CZ2','CZ3','N','ND1','ND2','NE','NE1','NE2','NH1','NH2','NZ','O','OD1','OD2','OE1','OE2','OG','OG1','OG2','OH','OXT','SD','SG');
  # Au{111}
  @tBA = (3.035,-4.14625,-4.08633,3758065.0,-7103.425,0.255010,2.320,-4.9565,-5.53098,89597.5,-1063.965,0.56486,-2.82);
  @tCM = (2.635,-16.15433,-15.51067,2667107.63333,-12099.23333,2.307,0.769,-4.42567,-4.94427,77316.66667,-942.7,0.56754,-30.40667);
  @tGA = (3.126,-8.749,-15.4325,9998178.25,-15054.2825,0.86918,2.286,-5.07925,-5.99153,77159.5,-989.425,0.68487,-15.6);
  @tCS = (3.172,-9.124,-4.9718,16937487.5,-23750.05,0.30467000,2.405,-7.4395,-8.54795,218325.5,-2066.65,0.96225,-1.53);
  @tIM = (3.31 ,-3.833,-3.21888,11943027.8,-12596.74,0.1913000,2.292,-4.4038,-5.32358,67309.2,-854.472,0.60499,-2.108);
  @tAM = (2.759,-39.646,-75.2894,11693610.0,-32938.07,4.627070,1.941,-26.501,-24.6274,58041.5,-1989.21,3.08458,-81.05);
  @tBE = (3.356,-3.17983,-2.6824,9952523.16667,-10497.28333,0.15942,2.223,-4.6575,-5.61167,49098.0,-742.23167,0.64519,-0.51167);
  @tME = (2.989,-4.1844,-4.17315,3914253.9,-7439.8,0.27999,2.431,-5.357,-5.9005,182161.5,-1624.76,0.66064,-4.175);
  @tCH = (3.494,-4.7807,-3.53935,22042382.0,-19579.7,0.21063,2.474,-8.6423,-8.87176,357731.0,-2947.04,0.99125,-0.05);
  $thresh=10.0;  # specific value for scoring function, value for gold
  $threshH=2.8;  # hydrophobicity acts only at the surface contact, 2x1.4 A means 2nd solvation
  $tag1="ATOM      1 AU1  SUR  ";
  $tag2="ATOM      2 AU2  SUR  ";
  $tag3="ATOM      3 AU3  SUR  ";
  $tag4="ATOM      4 AU4  SUR  ";
}
#-----------------------------------------------------------------------------------------------------------------------------
# Reading pdb and removing atoms which are not in the type table. 
#
open (FILE_INPUT, "< $pdb_input");
print "1/9 Read pdb and remove atoms which are not in the table...\t"; 
while (<FILE_INPUT>)
{
  @t = split(' ',$_);
  if ($t[0] eq 'ATOM')
    {
      for ($itype=0;$itype<38;$itype++)
      {
        if ($t[2] eq $type[$itype])
         {
             @tag[$i]=substr($_,0,30);
             @index[$i]=@t[1];
             @atom_type[$i]=@t[2];
             @res_type[$i]=@t[3];
             @res_numb[$i]=@t[4];
             if ($t[4] =~m/^\d+$/)  #pdb without chain ID (i.e. A, B..)
             {
               @xb[$i]=$t[5];
               @yb[$i]=$t[6];
               @zb[$i]=$t[7];
               $i=$i+1;
             }
             else{
               @xb[$i]=$t[6];        #pdb with a chain ID
               @yb[$i]=$t[7];
               @zb[$i]=$t[8];
               $i=$i+1;
             }
         }
      }
    }
}
close(FILE_INPUT);
$end=$i;
print "Done\n";
#-----------------------------------------------------------------------------------------------------------------------------
# Find the barycenter and make the translations
#
print "2/9 Find the barycenter and make translations...\t";
{
  $x_ori=origin($end,@xb);
  $y_ori=origin($end,@yb);
  $z_ori=origin($end,@zb);
  #printf ("PDB origin is at %.2f  %.2f  %.2f",$x_ori,$y_ori,$z_ori);
  for ($i=0;$i<$end;$i++)
  {
     ($xt[$i],$yt[$i],$zt[$i])=translate($xb[$i],$x_ori,$yb[$i],$y_ori,$zb[$i],$z_ori);
  }
}
print "Done\n";
#-----------------------------------------------------------------------------------------------------------------------------
# Euler proper rotations (precession then nutation) and scoring. A big output file is created: $label.mapfile
#
print "3/9 Euler proper rotations and score computations...\t";
$mapfile=$label."_score.map";
open (FMAP,"> $mapfile");
for ($prec=0;$prec<=180;$prec=$prec+$fine)
{
  for ($nuta=0;$nuta<=360;$nuta=$nuta+$fine)
  {
     for ($i=0;$i<$end;$i++)
     {
       ($XN,$YN,$ZN)=precession($xt[$i],$yt[$i],$zt[$i],$prec);
       ($x[$i],$y[$i],$z[$i])=nutation($XN,$YN,$ZN,$nuta);
     }
     #at this step, coordinates have been rotated by $prec and $nuta
     $xmin_value=minimum($end,@x);
     $ymin_value=minimum($end,@y);
     $xMax_value=maximum($end,@x);
     $yMax_value=maximum($end,@y);
     $zmin_value=minimum($end,@z);
     $zmin_test=$zmin_value+$thresh;
     $zmin_testH=$zmin_value+$threshH;
     # initializing the scores
     $score_FF=0.0;
     $score_QM=0.0;
     $score_hydro=0.0;
     # scoring
     for ($i=0;$i<$end;$i++)
     {
       if ($z[$i]<$zmin_test)
        {
           ($temp_FF,$temp_QM,$temp_hydro)=DLVO($z[$i],$threshH,$atom_type[$i],$res_type[$i],$res_numb[$i]);
           $score_FF=$score_FF+$temp_FF;
           $score_QM=$score_QM+$temp_QM;
           $score_hydro=$score_hydro+$temp_hydro;
        }
     }
     # write score in the map file
     $score_FFh=$score_FF-$score_hydro;
     $score_QMh=$score_QM-$score_hydro;
     printf FMAP ("%.2f    %.2f    %.2f    %.2f    %.2f     %.2f    %.2f    ",$prec,$nuta,$score_FF,$score_QM,$score_hydro,$score_FFh,$score_QMh);
     printf FMAP ("%.2f     %.2f     %.2f     %.2f     %.2f \n",$zmin_value,$xmin_value,$ymin_value,$xMax_value,$yMax_value);
  }
}
close(FMAP);
print "Done\n";
#------------------------------------------------------------------------------------------------------------------
# Reading the score.map file and make output scoring files
#
print "4/9 Reading the score.map file...\t";
{
$i=0;
open (Score_FILE, "<$mapfile");
while (<Score_FILE>)
 {
  @table = split(' ',$_);
  @prece[$i]=@table[0];
  @nutat[$i]=@table[1];
  @scoreHff[$i]=@table[2];
  @scoreHqm[$i]=@table[3];
  @scoreSolv[$i]=@table[4];
  @scoreGff[$i]=@table[5];
  @scoreGqm[$i]=@table[6];
  @zmin[$i]=@table[7];
  @xmin[$i]=@table[8];
  @ymin[$i]=@table[9];
  @xMax[$i]=@table[10];
  @yMax[$i]=@table[11];
  $i=$i+1;
 }
$endS=$i;
close(Score_FILE);
$score_Hff=$label."_Hff_rank.txt";
$score_Hqm=$label."_Hqm_rank.txt";
$score_Solv=$label."_Solv_rank.txt";
$score_Gff=$label."_Gff_rank.txt";
$score_Gqm=$label."_Gqm_rank.txt";
for ($i=0;$i<$endS;$i++)
 {
  $indices[$i]=$i;
 }
print "Done\n";
# Scoring conformations with Schwartzian Transformation
### conf option: provides the $display first energy conformations
print "5/9 Schwartzian transformation and selection...\t";
if ($clusconf eq "conf")
{
    # Scoring Hff
    open(FSELECT, ">$score_Hff");
    @indicestries = sort {$scoreHff[$a] <=> $scoreHff[$b]} (@indices) ;
    for ($i=0;$i<$display;$i++)
    {
      $tmp=$indicestries[$i];
      printf FSELECT ("%d      %.2f      %.2f      %.2f      %.2f \n",$i+1,$prece[$tmp],$nutat[$tmp],$zmin[$tmp],$scoreHff[$tmp]);
    }
    close(FSELECT);
    # Scoring Hqm
    open(FSELECT, ">$score_Hqm");
    @indicestries = sort {$scoreHqm[$a] <=> $scoreHqm[$b]} (@indices) ;
    for ($i=0;$i<$display;$i++)
    {
      $tmp=$indicestries[$i];
      printf FSELECT ("%d      %.2f      %.2f      %.2f      %.2f \n",$i+1,$prece[$tmp],$nutat[$tmp],$zmin[$tmp],$scoreHqm[$tmp]);
    }
    close(FSELECT);
    # Scoring Solv
    open(FSELECT, ">$score_Solv");
    @indicestries = sort {$scoreSolv[$a] <=> $scoreSolv[$b]} (@indices) ;
    for ($i=0;$i<$display;$i++)
    {
      $tmp=$indicestries[$i];
      printf FSELECT ("%d      %.2f      %.2f      %.2f      %.2f \n",$i+1,$prece[$tmp],$nutat[$tmp],$zmin[$tmp],$scoreSolv[$tmp]);
    }
    close(FSELECT);
    # Scoring Gff
    open(FSELECT, ">$score_Gff");
    @indicestries = sort {$scoreGff[$a] <=> $scoreGff[$b]} (@indices) ;
    for ($i=0;$i<$display;$i++)
    {
      $tmp=$indicestries[$i];
      printf FSELECT ("%d      %.2f      %.2f      %.2f      %.2f \n",$i+1,$prece[$tmp],$nutat[$tmp],$zmin[$tmp],$scoreGff[$tmp]);
    }
    close(FSELECT);
    # Scoring Gqm
    open(FSELECT, ">$score_Gqm");
    @indicestries = sort {$scoreGqm[$a] <=> $scoreGqm[$b]} (@indices) ;
    for ($i=0;$i<$display;$i++)
    {
      $tmp=$indicestries[$i];
      printf FSELECT ("%d      %.2f      %.2f      %.2f      %.2f \n",$i+1,$prece[$tmp],$nutat[$tmp],$zmin[$tmp],$scoreGqm[$tmp]);
    }
    close(FSELECT);
}
### clus option: provides the $display first clustered conformations
if ($clusconf eq "clus")
{
	# Scoring Hff
    open(FSELECT, ">$score_Hff");
    @indicestries = sort {$scoreHff[$a] <=> $scoreHff[$b]} (@indices) ;
 	$i=0;
    $tmp=$indicestries[$i];
    printf FSELECT ("%d      %.2f      %.2f      %.2f      %.2f \n",$i+1,$prece[$tmp],$nutat[$tmp],$zmin[$tmp],$scoreHff[$tmp]);
	$iclus=1;
	$j=0;
	while ($iclus < $display)
	{
        $cluster_oui_non="non";
		$j=$j+1;
		$distP=sqrt(($prece[$indicestries[$j]]-$prece[$indicestries[$i]])**2);
		$distN=sqrt(($nutat[$indicestries[$j]]-$nutat[$indicestries[$i]])**2);
		if ($distP >= 15 && $distP <= 165){$cluster_oui_non="oui";}
		if ($distN >= 15 && $distN <= 345){$cluster_oui_non="oui";}
		if ($cluster_oui_non eq "oui")
		{
			$tmp=$indicestries[$j];
			$iclus=$iclus+1;
			$i=0;
			$cluster_oui_non="non";
			printf FSELECT ("%d      %.2f      %.2f      %.2f      %.2f \n",$iclus,$prece[$tmp],$nutat[$tmp],$zmin[$tmp],$scoreHff[$tmp]);
			
        }
		if ($cluster_oui_non eq "non" )
		{
			
		}
	}
	close(FSELECT);
	# Scoring Hqm
    open(FSELECT, ">$score_Hqm");
    @indicestries = sort {$scoreHqm[$a] <=> $scoreHqm[$b]} (@indices) ;
 	$i=0;
    $tmp=$indicestries[$i];
    printf FSELECT ("%d      %.2f      %.2f      %.2f      %.2f \n",$i+1,$prece[$tmp],$nutat[$tmp],$zmin[$tmp],$scoreHqm[$tmp]);
	$iclus=1;
	$j=0;
	while ($iclus < $display)
	{
        $cluster_oui_non="non";
		$j=$j+1;
		$distP=sqrt(($prece[$indicestries[$j]]-$prece[$indicestries[$i]])**2);
		$distN=sqrt(($nutat[$indicestries[$j]]-$nutat[$indicestries[$i]])**2);
		if ($distP >= 15 && $distP <= 165){$cluster_oui_non="oui";}
		if ($distN >= 15 && $distN <= 345){$cluster_oui_non="oui";}
		if ($cluster_oui_non eq "oui")
		{
			$tmp=$indicestries[$j];
			$iclus=$iclus+1;
			$i=0;
			$cluster_oui_non="non";
			printf FSELECT ("%d      %.2f      %.2f      %.2f      %.2f \n",$iclus,$prece[$tmp],$nutat[$tmp],$zmin[$tmp],$scoreHqm[$tmp]);
			
        }
		if ($cluster_oui_non eq "non" )
		{
			
		}
	}
	close(FSELECT);
    # Scoring Solv
    open(FSELECT, ">$score_Solv");
    @indicestries = sort {$scoreSolv[$a] <=> $scoreSolv[$b]} (@indices) ;
 	$i=0;
    $tmp=$indicestries[$i];
    printf FSELECT ("%d      %.2f      %.2f      %.2f      %.2f \n",$i+1,$prece[$tmp],$nutat[$tmp],$zmin[$tmp],$scoreSolv[$tmp]);
	$iclus=1;
	$j=0;
	while ($iclus < $display)
	{
        $cluster_oui_non="non";
		$j=$j+1;
		$distP=sqrt(($prece[$indicestries[$j]]-$prece[$indicestries[$i]])**2);
		$distN=sqrt(($nutat[$indicestries[$j]]-$nutat[$indicestries[$i]])**2);
		if ($distP >= 15 && $distP <= 165){$cluster_oui_non="oui";}
		if ($distN >= 15 && $distN <= 345){$cluster_oui_non="oui";}
		if ($cluster_oui_non eq "oui")
		{
			$tmp=$indicestries[$j];
			$iclus=$iclus+1;
			$i=0;
			$cluster_oui_non="non";
			printf FSELECT ("%d      %.2f      %.2f      %.2f      %.2f \n",$iclus,$prece[$tmp],$nutat[$tmp],$zmin[$tmp],$scoreSolv[$tmp]);
			
        }
		if ($cluster_oui_non eq "non" )
		{
			
		}
	}
	close(FSELECT);    
	# Scoring Gff
    open(FSELECT, ">$score_Gff");
        @indicestries = sort {$scoreGff[$a] <=> $scoreGff[$b]} (@indices) ;
 	$i=0;
    $tmp=$indicestries[$i];
    printf FSELECT ("%d      %.2f      %.2f      %.2f      %.2f \n",$i+1,$prece[$tmp],$nutat[$tmp],$zmin[$tmp],$scoreGff[$tmp]);
	$iclus=1;
	$j=0;
	while ($iclus < $display)
	{
        $cluster_oui_non="non";
		$j=$j+1;
		$distP=sqrt(($prece[$indicestries[$j]]-$prece[$indicestries[$i]])**2);
		$distN=sqrt(($nutat[$indicestries[$j]]-$nutat[$indicestries[$i]])**2);
		if ($distP >= 15 && $distP <= 165){$cluster_oui_non="oui";}
		if ($distN >= 15 && $distN <= 345){$cluster_oui_non="oui";}
		if ($cluster_oui_non eq "oui")
		{
			$tmp=$indicestries[$j];
			$iclus=$iclus+1;
			$i=0;
			$cluster_oui_non="non";
			printf FSELECT ("%d      %.2f      %.2f      %.2f      %.2f \n",$iclus,$prece[$tmp],$nutat[$tmp],$zmin[$tmp],$scoreGff[$tmp]);
			
        }
		if ($cluster_oui_non eq "non" )
		{
			
		}
	}
	close(FSELECT);
	# Scoring Gqm
    open(FSELECT, ">$score_Gqm");
	    @indicestries = sort {$scoreGqm[$a] <=> $scoreGqm[$b]} (@indices) ;
 	$i=0;
    $tmp=$indicestries[$i];
    printf FSELECT ("%d      %.2f      %.2f      %.2f      %.2f \n",$i+1,$prece[$tmp],$nutat[$tmp],$zmin[$tmp],$scoreGqm[$tmp]);
	$iclus=1;
	$j=0;
	while ($iclus < $display)
	{
        $cluster_oui_non="non";
		$j=$j+1;
		$distP=sqrt(($prece[$indicestries[$j]]-$prece[$indicestries[$i]])**2);
		$distN=sqrt(($nutat[$indicestries[$j]]-$nutat[$indicestries[$i]])**2);
		if ($distP >= 15 && $distP <= 165){$cluster_oui_non="oui";}
		if ($distN >= 15 && $distN <= 345){$cluster_oui_non="oui";}
		if ($cluster_oui_non eq "oui")
		{
			$tmp=$indicestries[$j];
			$iclus=$iclus+1;
			$i=0;
			$cluster_oui_non="non";
			printf FSELECT ("%d      %.2f      %.2f      %.2f      %.2f \n",$iclus,$prece[$tmp],$nutat[$tmp],$zmin[$tmp],$scoreGqm[$tmp]);
			
        }
		if ($cluster_oui_non eq "non" )
		{
			
		}
	}
	close(FSELECT);
}
}
print "Done\n";
#------------------------------------------------------------------------------------------------------------------
# Generate a maximal surface for displaying all conformations
#
print "6/9 Generate a maximal $surf_type surface for displaying all conformations...\t";
if ($surf_type eq "Au111" )
{
#$mX=minimum($end,@x);
#$mY=minimum($end,@y);
#$mZ=minimum($end,@z);
#$MX=maximum($end,@x);
#$MY=maximum($end,@y);
$DistX=maximum($endS,@xMax)-minimum($endS,@xmin);
$DistY=maximum($endS,@yMax)-minimum($endS,@ymin);
$ZminGlob=minimum($endS,@zmin);
$Num_Unit_X=int(2*$DistX/(2*2.039))+1;
$Num_Unit_Y=int(2*$DistY/(2*2.039))+1;
#$XNo=$mX-0.25*$dX;
#$YNo=$mY-0.25*$dY;
#$ZN=$mZ-3.0; #could be tunable
$XNo=minimum($endS,@xmin);
$YNo=minimum($endS,@ymin);
$ZN=$ZminGlob - 3.0; #could be tunable
$pdb_surface=$label."-".$surf_type.".pdb";
open (FPDB,"> $pdb_surface");
for ($i=0;$i<$Num_Unit_X;$i++)
  {
    $XN=$XNo+$i*2.039*2;
    for($j=0;$j<$Num_Unit_Y;$j++)
    {
     $resid=$resid+1;
     if ($resid<1000){$espace="";}
     if ($resid<100){$espace=" ";}
     if ($resid<10){$espace="  ";}
     $Ntag1=$tag1.$espace.$resid." ";
     $Ntag2=$tag2.$espace.$resid." ";
     $Ntag3=$tag3.$espace.$resid." ";
     $Ntag4=$tag4.$espace.$resid." ";
     $YN=$YNo+$j*2.039*2;
     writepdb(($XN-2.039),($YN-2.039),($ZN-2.039),$Ntag1);
     writepdb(($XN-2.039),($YN),($ZN),$Ntag2);
     writepdb(($XN),($YN-2.039),($ZN),$Ntag3);
     writepdb(($XN),($YN),($ZN-2.039),$Ntag4);
    }
   }
close(FPDB);
}
print "Done\n";
#-----------------------------------------------------------------------------------------------------------------------------
# Write a gnuplot input file and generate the various plots
#
print "7/9 Generate various map plots...\t";
{
	$row_grid=2*(int(180/$fine)+1);
    $col_grid=2*(int(360/$fine)+1);
    $scriptplotfile=$label."_plot.plt";
    #$image_Hmm=$label."-Hmm.png";
    #$image_Hqm=$label."-Hqm.png";
    #$image_Hsolv=$label."-Hsolv.png";
    #$image_Gmm=$label."-Gmm.png";
    #$image_Gqm=$label."-Gqm.png";
	open (FPLOT, "> $scriptplotfile");
    print FPLOT "set dgrid3d $row_grid,$col_grid splines\n";
    print FPLOT "set xlabel 'Precession'\n";
    print FPLOT "set ylabel 'Nutation'\n";
    print FPLOT "set xrange [0:180]\n";
    print FPLOT "set yrange [0:360.1]\n";
    print FPLOT "set xtics 30 \n";
    print FPLOT "set ytics 30 \n";
    print FPLOT "set pm3d at b\n";
    print FPLOT "unset surface\n";
    print FPLOT "unset key\n";
    print FPLOT "unset clabel\n";
    print FPLOT "set size square\n";
    print FPLOT "set view map\n";
    print FPLOT "set contour base\n";
    print FPLOT "set palette rgbformulae 33,13,10\n";
    print FPLOT "set term png enhanced size 1000,1000 font \"arial,20\"\n";
    @title_gnuplot=("Hmm","Hqm","Solvation","Gmm","Gqm");
	@name_scoring=($score_Hff,$score_Hqm,$score_Solv,$score_Gff,$score_Gqm);
	$col=3;
    for($ind=0;$ind<5;$ind++)
    {
      $col=3+$ind;
	  $image=$label.$title_gnuplot[$ind].".png";
	  print FPLOT "set title \"$label $title_gnuplot[$ind] contribution\"\n";
	  open (FILE_INPUT, "< $name_scoring[$ind]");
      while (<FILE_INPUT>)
      {
        @t = split(' ',$_);
		print FPLOT "set label $t[0] \"$t[0]\" at $t[1],$t[2],1 font \"arial bold,20\" textcolor \"black\" front \n";  
      }
	  print FPLOT "set output \"$image\"\n";
	  print FPLOT "splot \"$mapfile\" u 1:2:$col w li lt -1 lw 2\n";
	  
    }
    print FPLOT "quit\n";
    close(FPLOT);
    system ("gnuplot < $scriptplotfile");
}
print "Done\n";
#-----------------------------------------------------------------------------------------------------------------------------
# Write all final pdb
#
print "8/9 Write all final pdb...\t";
@name_scoring=($score_Hff,$score_Hqm,$score_Solv,$score_Gff,$score_Gqm);
for($ind=0;$ind<5;$ind++)
 {
   open (FILE_INPUT, "< $name_scoring[$ind]");
   #print "Opening $name_scoring[$ind]";
   while (<FILE_INPUT>)
   {
     @t = split(' ',$_);
	 $j=$t[0];
     $outpdb=$name_scoring[$ind];
	 $outpdb=~ s/_rank.txt//;
	 $outpdb=$outpdb."-".$j."-P".$t[1]."N".$t[2].".pdb";
	 open (FPDB,"> $outpdb");
     for ($i=0;$i<$end;$i++)
     {
        ($XN,$YN,$ZN)=precession($xt[$i],$yt[$i],$zt[$i],$t[1]);
        ($x[$i],$y[$i],$z[$i])=nutation($XN,$YN,$ZN,$t[2]);
     }
	 for ($i=0;$i<$end;$i++)
     { 
       writepdb(($x[$i]+$DistX/2.0),($y[$i]+$DistY/2.0),($z[$i]-($t[3]-$ZminGlob)),$tag[$i]);
     }
     print FPDB "END\n";
     close(FPDB);
   }
 }
 print "Done\n";
 print "9/9 END. Please cite .....\n";
#-----------------------------------------------------------------------------------------------------------------------------
#############################
#                           #
#  FUNCTIONS DECLARATTIONS  #
#                           #
#############################
#
#
######################## Function origin  ###########################
sub origin {                                                        #
    my ($fin,@table)=@_;                                            #
    my $val=$table[0];                                              #
    my $ind=1;                                                      #
    for ($ind=1;$ind<$fin;$ind++)                                   #
    {                                                               #
      $val=$val+$table[$ind];                                       #
    }                                                               #
    $val=$val/$fin;                                                 #
    return ($val);                                                  #
    }                                                               #
#####################################################################
######################## Function writepdb ##########################
sub writepdb {                                                      #
    my ($a,$b,$c,$t)=@_;                                            #
    my $text="";                                                    #
    my $pad_len=30;                                                 #
    $text = sprintf ("%-${pad_len}s%8.2f%8.2f%8.2f \n",$t,$a,$b,$c);#
    print FPDB "$text";                                             #
    return();                                                       #
    }                                                               #
#####################################################################
######################## Fonction translate #########################
sub translate {                                                     #
    my ($a,$at,$b,$bt,$c,$ct)=@_;                                   #
    my $new_x=$a-$at;                                               #
    my $new_y=$b-$bt;                                               #
    my $new_z=$c-$ct;                                               #
    return ($new_x,$new_y,$new_z);                                  #
    }                                                               #
#####################################################################
######################## Function nutation ##########################
sub nutation {                                                      #
    my ($nut_a,$nut_b,$nut_c,$T)=@_;                                #
    my $nut_X=0;                                                    #
    my $nut_Y=0;                                                    #
    my $nut_Z=0;                                                    #
    $T=($T*3.141592654/180.0);                                      #
    $nut_X=$nut_a;                                                  #
    $nut_Y=$nut_b*cos($T)-$nut_c*sin($T);                           #
    $nut_Z=$nut_b*sin($T)+$nut_c*cos($T);                           #
    return ($nut_X,$nut_Y,$nut_Z);                                  #
    }                                                               #
#####################################################################
######################## Function precession ########################
sub precession {                                                    #
    my ($prec_a,$prec_b,$prec_c,$T)=@_;                             #
    my $prec_X=0;                                                   #
    my $prec_Y=0;                                                   #
    my $prec_Z=0;                                                   #
    $T=($T*3.141592654/180.0);                                      #
    $prec_X=$prec_a*cos($T)-$prec_c*sin($T);                        #
    $prec_Y=$prec_b;                                                #
    $prec_Z=$prec_a*sin($T)+$prec_c*cos($T);                        #
    return ($prec_X,$prec_Y,$prec_Z);                               #
    }                                                               #
#####################################################################
######################## Function minimum ###########################
sub minimum {                                                       #
    my ($fin,@table)=@_;                                            #
    my $min=$table[0];                                              #
    my $ind=1;                                                      #
    for ($ind=1;$ind<$fin;$ind++)                                   #
    {                                                               #
      if ($table[$ind]<$min){$min=$table[$ind];}                    #
    }                                                               #
    return ($min);                                                  #
    }                                                               #
#####################################################################
######################## Function maximum ###########################
sub maximum {                                                       #
    my ($fin,@table)=@_;                                            #
    my $max=$table[0];                                              #
    my $ind=1;                                                      #
    for ($ind=1;$ind<$fin;$ind++)                                   #
    {                                                               #
      if ($table[$ind]>$max){$max=$table[$ind];}                    #
    }                                                               #
    return ($max);                                                  #
    }                                                               #
#####################################################################
######################## Function DLVO ##############################
sub DLVO {                                                          #
    my ($d,$h,$Atype,$Restype,$numres)=@_;                          #
    my $FF=0;                                                       #
    my $QM=0;                                                       #
    my $solv=0;                                                     #
    my @tr=();                                                      #
    if ($Atype eq "CA"){@tr=@tBA;}                                  #
    if ($Atype eq "N"){@tr=@tBA;}                                   #
    if ($numres == 1){if ($Atype eq "N"){@tr=@tAM;}}                #
    if ($Atype eq "C"){@tr=@tBA;}                                   #
    if ($Atype eq "O"){@tr=@tBA;}                                   #
    if ($Atype eq "CB"){@tr=@tCH;} #by default                      #
    if ($Restype eq "ASN")                                          #
      {                                                             #
      if ($Atype eq "CB"||"CG"||"ND2"||"OD1")                       #
      {@tr=@tBA;}                                                   #
      }                                                             #
    if ($Restype eq "GLN")                                          #
      {                                                             #
      if ($Atype eq "CD"||"CG"||"NE2"||"OE1")                       #
      {@tr=@tBA;}                                                   #
      }                                                             #
    if ($Restype eq "ASP")                                          #
      {                                                             #
      if ($Atype eq "CG"||"OD1"||"OD2")                             #
      {@tr=@tCM;}                                                   #
      }                                                             #
    if ($Restype eq "GLU")                                          #
      {                                                             #
      if ($Atype eq "CD"||"OE1"||"OE2")                             #
      {@tr=@tCM;}                                                   #
      }                                                             #
    if ($Restype eq "ARG")                                          #
      {                                                             #
      if ($Atype eq "NE"||"CZ"||"NH2"||"NH1")                       #
      {@tr=@tGA;}                                                   #
      }                                                             #
    if ($Restype eq "CYS"||"CYX")                                   #
      {                                                             #
      if ($Atype eq "SG"||"CB")                                     #
      {@tr=@tCS;}                                                   #
      }                                                             #
    if ($Restype eq "MET")                                          #
      {                                                             #
      if ($Atype eq "CG"||"SD"||"CE")                               #
      {@tr=@tCS;}                                                   #
      }                                                             #
    if ($Restype eq "HIS"||"HIE"||"HID"||"HIP")                     #
      {                                                             #
      if ($Atype eq "CG"||"CD2"||"ND1"||"CE1"||"NE2")               #
      {@tr=@tIM;}                                                   #
      }                                                             #
    if ($Restype eq "TRP")                                          #
      {                                                             #
      if ($Atype eq "CG"||"CD1"||"NE1")                             #
      {@tr=@tIM;}                                                   #
      }                                                             #
    if ($Restype eq "PHE"||"TYR")                                   #
      {                                                             #
      if ($Atype eq "CG"||"CD1"||"CD2"||"CE1"||"CE2"||"CZ")         #
      {@tr=@tBE;}                                                   #
      }                                                             #
     if ($Restype eq "TRP")                                         #
      {                                                             #
      if ($Atype eq "CD2"||"CE2"||"CE3"||"CZ2"||"CZ3"||"CH2")       #
      {@tr=@tBE;}                                                   #
      }                                                             #
     if ($Restype eq "SER")                                         #
      {                                                             #
      if ($Atype eq "CB"||"OG")                                     #
      {@tr=@tME;}                                                   #
      }                                                             #
     if ($Restype eq "THR")                                         #
      {                                                             #
      if ($Atype eq "CB"||"OG1"||"CG2")                             #
      {@tr=@tME;}                                                   #
      }                                                             #
     if ($Restype eq "TYR")                                         #
      {                                                             #
      if ($Atype eq "OH")                                           #
      {@tr=@tME;}                                                   #
      }                                                             #
     if ($Restype eq "ARG"||"GLN"||"GLU"||"LEU"||"LYS"||"PRO")      #
      {                                                             #
      if ($Atype eq "CG")                                           #
      {@tr=@tCH;}                                                   #
      }                                                             #
     if ($Restype eq "ILE"||"VAL")                                  #
      {                                                             #
      if ($Atype eq "CG1"||"CG2")                                   #
      {@tr=@tCH;}                                                   #
      }                                                             #
     if ($Restype eq "ARG"||"LYS"||"PRO")                           #
      {                                                             #
      if ($Atype eq "CD")                                           #
      {@tr=@tCH;}                                                   #
      }                                                             #
     if ($Restype eq "ILE"||"LEU")                                  #
      {                                                             #
      if ($Atype eq "CD1"||"CD2")                                   #
      {@tr=@tCH;}                                                   #
      }                                                             #
    if ($Atype eq "OXT"){@tr=@tCM;}                                 #
    if ($Atype eq "NZ"){@tr=@tAM;}                                  #
    if ($d < $tr[0]){$QM=$tr[1];}                                   #
    if ($d >= $tr[0])                                               #
    {                                                               #
     $QM=($tr[2]/$d+($tr[3]/$d**12)+($tr[4]/$d**6)+$tr[5]);         #
     if($QM>0){$QM=0}                                               #
    }                                                               #
    if ($d < $tr[6]){$FF=$tr[7];}                                   #
    if ($d >= $tr[6])                                               #
    {                                                               #
    $FF=($tr[8]/$d+($tr[9]/$d**12)+($tr[10]/$d**6)+$tr[11])/$tr[0]; #
    if($FF>0){$FF=0}                                                #
    }                                                               #
    if($d<$h){$solv=$tr[12];}                                       #
    return ($FF,$QM,$solv);                                         #
    }                                                               #
#####################################################################
######################## Function radian ############################
sub radian {                                                        #
    my ($a)=@_;                                                     #
    my $new=$a/180.0*3.141592654;                                   #
    return ($new);                                                  #
    }                                                               # 
#####################################################################







