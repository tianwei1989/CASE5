#!MC 1200
# Created by Tecplot 360 build 12.0.0.3116
$!VarSet |MFBD| = '.'
$!READDATASET  '.\result.plt '

$!GLOBALTHREED SLICE{NORMAL{X = 1}}
$!GLOBALTHREED SLICE{NORMAL{Z = 0}}
$!GLOBALTHREED SLICE{ORIGIN{X = 1.22}}
$!CREATESLICEZONEFROMPLANE 
  SLICESOURCE = VOLUMEZONES
  FORCEEXTRACTIONTOSINGLEZONE = YES
  COPYCELLCENTEREDVALUES = NO
$!PLOTTYPE = CARTESIAN2D
$!TWODAXIS XDETAIL{VARNUM = 2}
$!TWODAXIS YDETAIL{VARNUM = 3}
$!GLOBALTWODVECTOR UVAR = 8
$!GLOBALTWODVECTOR VVAR = 9
$!RESETVECTORLENGTH 
$!FIELDLAYERS SHOWVECTOR = YES
$!RemoveVar |MFBD|