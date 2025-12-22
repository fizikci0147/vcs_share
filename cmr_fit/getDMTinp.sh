#!/bin/bash

Q2=$1
W=$2
THETA=$3
TH=`echo "180-$THETA" | bc`

rm -f amp0a.htm amp1a.htm amp0b.htm amp1b.htm cglns.htm

# l=0, (p(1/2), n(1/2), 3/2)
wget -q -O amp0a.htm "http://maid.kph.uni-mainz.de/cgi-bin/maid1?switch=104&param0=1&value1=0&param1=1&param3=1&param4=1&param5=1&param6=1&param11=1&param12=1&param13=1&param14=1&param15=1&param16=1&param17=1&param18=1&param19=1&param20=1&param21=1&value51=1.0&value52=1.0&value53=1.0&value54=1.0&value55=0.0&value56=1.0&value57=1.0&value58=0.0&value59=1.0&value60=1.0&value61=1.0&value62=1.0&value63=0.0&value64=1.0&value65=0.0&value66=1.0&value67=0.0&value68=1.0&value69=1.0&value70=0.0&value71=1.0&value72=1.0&value73=1.0&value74=1.0&value75=0.0&value76=1.0&value77=0.0&value78=1.0&value79=0.0&value80=1.0&value81=1.0&value82=0.0&param50=2&value35=${Q2}&value36=${W}&value41=1000&value42=1300"
# l=1, (p(1/2), n(1/2), 3/2)
wget -q -O amp1a.htm "http://maid.kph.uni-mainz.de/cgi-bin/maid1?switch=104&param0=1&value1=1&param1=1&param3=1&param4=1&param5=1&param6=1&param11=1&param12=1&param13=1&param14=1&param15=1&param16=1&param17=1&param18=1&param19=1&param20=1&param21=1&value51=1.0&value52=1.0&value53=1.0&value54=1.0&value55=0.0&value56=1.0&value57=1.0&value58=0.0&value59=1.0&value60=1.0&value61=1.0&value62=1.0&value63=0.0&value64=1.0&value65=0.0&value66=1.0&value67=0.0&value68=1.0&value69=1.0&value70=0.0&value71=1.0&value72=1.0&value73=1.0&value74=1.0&value75=0.0&value76=1.0&value77=0.0&value78=1.0&value79=0.0&value80=1.0&value81=1.0&value82=0.0&param50=2&value35=${Q2}&value36=${W}&value41=1000&value42=1300"
#l=0 charge channels
wget -q -O amp0b.htm "http://maid.kph.uni-mainz.de/cgi-bin/maid1?switch=104&param0=3&value1=0&param1=1&param3=1&param4=1&param5=1&param6=1&param11=1&param12=1&param13=1&param14=1&param15=1&param16=1&param17=1&param18=1&param19=1&param20=1&param21=1&value51=1.0&value52=1.0&value53=1.0&value54=1.0&value55=0.0&value56=1.0&value57=1.0&value58=0.0&value59=1.0&value60=1.0&value61=1.0&value62=1.0&value63=0.0&value64=1.0&value65=0.0&value66=1.0&value67=0.0&value68=1.0&value69=1.0&value70=0.0&value71=1.0&value72=1.0&value73=1.0&value74=1.0&value75=0.0&value76=1.0&value77=0.0&value78=1.0&value79=0.0&value80=1.0&value81=1.0&value82=0.0&param50=2&value35=${Q2}&value36=${W}&value41=1000&value42=1300"
#l=1 charge channels
wget -q -O amp1b.htm "http://maid.kph.uni-mainz.de/cgi-bin/maid1?switch=104&param0=3&value1=1&param1=1&param3=1&param4=1&param5=1&param6=1&param11=1&param12=1&param13=1&param14=1&param15=1&param16=1&param17=1&param18=1&param19=1&param20=1&param21=1&value51=1.0&value52=1.0&value53=1.0&value54=1.0&value55=0.0&value56=1.0&value57=1.0&value58=0.0&value59=1.0&value60=1.0&value61=1.0&value62=1.0&value63=0.0&value64=1.0&value65=0.0&value66=1.0&value67=0.0&value68=1.0&value69=1.0&value70=0.0&value71=1.0&value72=1.0&value73=1.0&value74=1.0&value75=0.0&value76=1.0&value77=0.0&value78=1.0&value79=0.0&value80=1.0&value81=1.0&value82=0.0&param50=2&value35=${Q2}&value36=${W}&value41=1000&value42=1300"
# CGLN Fs
wget -q -O cglns.htm "http://maid.kph.uni-mainz.de/cgi-bin/maid1?switch=103&param0=1&param1=2&param11=1&param12=1&param13=1&param14=1&param15=1&param16=1&param17=1&param18=1&param19=1&param20=1&param21=1&value51=1.0&value52=1.0&value53=1.0&value54=1.0&value55=0.0&value56=1.0&value57=1.0&value58=0.0&value59=1.0&value60=1.0&value61=1.0&value62=1.0&value63=0.0&value64=1.0&value65=0.0&value66=1.0&value67=0.0&value68=1.0&value69=1.0&value70=0.0&param50=3&value35=${Q2}&value36=${W}&value37=${TH}&value41=1000&value42=180"



e0pr=`awk 'NR==56 {print $2;exit}' amp0b.htm`
e0pi=`awk 'NR==56 {print $3;exit}' amp0b.htm`
l0pr=`awk 'NR==68 {print $2;exit}' amp0b.htm`
l0pi=`awk 'NR==68 {print $3;exit}' amp0b.htm`

e1pr=`awk 'NR==56 {print $2;exit}' amp1b.htm`
e1pi=`awk 'NR==56 {print $3;exit}' amp1b.htm`
m1pr=`awk 'NR==60 {print $2;exit}' amp1b.htm`
m1pi=`awk 'NR==60 {print $3;exit}' amp1b.htm`
m1mr=`awk 'NR==64 {print $2;exit}' amp1b.htm`
m1mi=`awk 'NR==64 {print $3;exit}' amp1b.htm`
l1pr=`awk 'NR==68 {print $2;exit}' amp1b.htm`
l1pi=`awk 'NR==68 {print $3;exit}' amp1b.htm`
l1mr=`awk 'NR==72 {print $2;exit}' amp1b.htm`
l1mi=`awk 'NR==72 {print $3;exit}' amp1b.htm`

e0p1r=`awk 'NR==56 {print $2;exit}' amp0a.htm`
e0p1i=`awk 'NR==56 {print $3;exit}' amp0a.htm`
e0p3r=`awk 'NR==56 {print $6;exit}' amp0a.htm`
e0p3i=`awk 'NR==56 {print $7;exit}' amp0a.htm`
l0p1r=`awk 'NR==68 {print $2;exit}' amp0a.htm`
l0p1i=`awk 'NR==68 {print $3;exit}' amp0a.htm`
l0p3r=`awk 'NR==68 {print $6;exit}' amp0a.htm`
l0p3i=`awk 'NR==68 {print $7;exit}' amp0a.htm`

e1p1r=`awk 'NR==56 {print $2;exit}' amp1a.htm`
e1p1i=`awk 'NR==56 {print $3;exit}' amp1a.htm`
e1p3r=`awk 'NR==56 {print $6;exit}' amp1a.htm`
e1p3i=`awk 'NR==56 {print $7;exit}' amp1a.htm`
m1p1r=`awk 'NR==60 {print $2;exit}' amp1a.htm`
m1p1i=`awk 'NR==60 {print $3;exit}' amp1a.htm`
m1p3r=`awk 'NR==60 {print $6;exit}' amp1a.htm`
m1p3i=`awk 'NR==60 {print $7;exit}' amp1a.htm`
m1m1r=`awk 'NR==64 {print $2;exit}' amp1a.htm`
m1m1i=`awk 'NR==64 {print $3;exit}' amp1a.htm`
m1m3r=`awk 'NR==64 {print $6;exit}' amp1a.htm`
m1m3i=`awk 'NR==64 {print $7;exit}' amp1a.htm`
l1p1r=`awk 'NR==68 {print $2;exit}' amp1a.htm`
l1p1i=`awk 'NR==68 {print $3;exit}' amp1a.htm`
l1p3r=`awk 'NR==68 {print $6;exit}' amp1a.htm`
l1p3i=`awk 'NR==68 {print $7;exit}' amp1a.htm`
l1m1r=`awk 'NR==72 {print $2;exit}' amp1a.htm`
l1m1i=`awk 'NR==72 {print $3;exit}' amp1a.htm`
l1m3r=`awk 'NR==72 {print $6;exit}' amp1a.htm`
l1m3i=`awk 'NR==72 {print $7;exit}' amp1a.htm`

F1r=`awk 'NR==53 {print $2;exit}' cglns.htm`
F1i=`awk 'NR==53 {print $3;exit}' cglns.htm`
F2r=`awk 'NR==53 {print $4;exit}' cglns.htm`
F2i=`awk 'NR==53 {print $5;exit}' cglns.htm`
F3r=`awk 'NR==53 {print $6;exit}' cglns.htm`
F3i=`awk 'NR==53 {print $7;exit}' cglns.htm`
F4r=`awk 'NR==57 {print $2;exit}' cglns.htm`
F4i=`awk 'NR==57 {print $3;exit}' cglns.htm`
F5r=`awk 'NR==57 {print $4;exit}' cglns.htm`
F5i=`awk 'NR==57 {print $5;exit}' cglns.htm`
F6r=`awk 'NR==57 {print $6;exit}' cglns.htm`
F6i=`awk 'NR==57 {print $7;exit}' cglns.htm`


echo "99  $Q2  $W  $THETA  phi  eps  sig  err   0   0   $e0pr $e0pi  $l0pr $l0pi  $e1p3r $e1p3i  $m1p3r $m1p3i  $l1p3r $l1p3i  $m1mr $m1mi  $l1mr $l1mi   $F1r $F1i  $F2r $F2i  $F3r $F3i  $F4r $F4i  $F5r $F5i  $F6r $F6i  $e1pr $e1pi  $m1pr $m1pi  $l1pr $l1pi  $e0p3r $e0p3i  $l0p3r  $l0p3i  $m1m3r $m1m3i  $l1m3r $l1m3i  $e0p1r $e0p1i  $l0p1r $l0p1i  $e1p1r $e1p1i  $m1p1r $m1p1i  $l1p1r $l1p1i  $m1m1r $m1m1i  $l1m1r $l1m1i"

rm amp0a.htm amp1a.htm amp0b.htm amp1b.htm cglns.htm
