//#######################
// parameters
//##############
// span of one wing
halfSpan = 1.;
// aspect ratio
aspectRatio = 3.1;
// angle of attack (deg)
alpha = 40;
// position of centre of rotation
centreX = 0.0;
centreY = 0.0;
centreZ = 0.0;
// relative thickness
thicknessToHalfSpan = 0.0107;
//##############
// sizes around planform
outerSide = 2.0  * halfSpan;
innerSide = 0.05 * halfSpan;
top       = 2.5  * halfSpan;
bottom    = 2.5  * halfSpan;

chord     = halfSpan / aspectRatio;
thickness = halfSpan * thicknessToHalfSpan;

//##############
// MESH D - approx 16.7M nodes - for full circle
numEl_front     = 50;
numEl_chord     = 30;
numEl_back      = 70;
numEl_quadrant  = 60;
numEl_outerSide = 30;
numEl_innerSide = 20;
numEl_span      = 30;
numEl_thickness = 5;
numEl_bottom    = 50;
numEl_top       = 50;
// progressions
progression_front  = 0.98; 
progression_chord  = 1.; 
progression_back   = 1.01; 
progression_quadrant = 1.; 
progression_outer  = 0.93; 
progression_inner  = 0.2; 
progression_span   = 1.00; 
progression_top    = 1.07; 
progression_bottom = 0.93; 

// number of blocks
npointi = 5; // tangential direction
npointj = 4; // radial direction
npointk = 4; // axial direction

stepj = npointi;
stepk = 1000;

// layer of inclusion of plate
k = 2;

//Enclose the circle?
EncloseCircle = 1;
//#######################

//#######################
// auxiliary functions 
Function ProjectPoint 

  // In the following commands we use the reserved variable name
  // `newp', which automatically selects a new point number. This
  // number is chosen as the highest current point number, plus
  // one. (Note that, analogously to `newp', the variables `newc',
  // `news', `newv' and `newreg' select the highest number amongst
  // currently defined curves, surfaces, volumes and `any entities
  // other than points', respectively.)

  // expects: 
  // norm1, norm2, norm3 - plane normal components
  // dir1, dir2, dir3 - direction of projection
  // xp, yp, zp - a point in the plane to which we project
  // oldnumber - old number of the point
  // newnumber - new number of the point

  x[] = Point{oldnumber};
  c = (norm1 * (xp - x[0]) + norm2 * (yp - x[1]) + norm3 * (zp - x[2])) / (norm1*dir1 + norm2*dir2 + norm3*dir3);

  Point(newnumber) = {x[0] + c*dir1,  x[1] + c*dir2,  x[2] + c*dir3} ;

Return

Function RotatePoint
Rotate {{0, 0, 1}, {centreX, centreY, centreZ}, -alpha/180*Pi} {
   Point{pointIndex};
}
Return

Function ComputeEllipsePoint
   // on input: 
   //   a_ellipse
   //   b_ellipse
   //   centreX
   //   centreZ
   //   x_coord
   // on output:
   //   z_coord

   z_coord = centreZ - Sqrt( b_ellipse^2 * ( 1 - (x_coord - centreX)^2 / a_ellipse^2 ));
Return
//#######################

// planform in 2D
//#######################
OriginX = 0;
OriginY = 0;
OriginZ = 1.0667;

 Point(         (k-1)*stepk+101 ) = {  -9.23075999999999931E-003 ,    0.0000000000000000      ,   3.07691999999999948E-003  };
 Point(         (k-1)*stepk+1*stepj+2 ) = {  -7.69230000000000053E-002 ,    0.0000000000000000      ,   4.92307199999999917E-002  };
 Point(         (k-1)*stepk+103 ) = {  -0.11384604000000000      ,    0.0000000000000000      ,   9.53845199999999865E-002  };
 Point(         (k-1)*stepk+104 ) = {  -0.13230755999999999      ,    0.0000000000000000      ,   0.14153832000000000       };
 Point(         (k-1)*stepk+105 ) = {  -0.14461524000000001      ,    0.0000000000000000      ,   0.18769212000000002       };
 Point(         (k-1)*stepk+106 ) = {  -0.15076908000000000      ,    0.0000000000000000      ,   0.23384591999999998       };
 Point(         (k-1)*stepk+107 ) = {  -0.15692292000000002      ,    0.0000000000000000      ,   0.27999972000000001       };
 Point(         (k-1)*stepk+108 ) = {  -0.15692292000000002      ,    0.0000000000000000      ,   0.32615352000000003       };
 Point(         (k-1)*stepk+109 ) = {  -0.15692292000000002      ,    0.0000000000000000      ,   0.37230732000000000       };
 Point(         (k-1)*stepk+110 ) = {  -0.15076908000000000      ,    0.0000000000000000      ,   0.41846112000000002       };
 Point(         (k-1)*stepk+111 ) = {  -0.14461524000000003      ,    0.0000000000000000      ,   0.46461491999999999       };
 Point(         (k-1)*stepk+112 ) = {  -0.13846140000000004      ,    0.0000000000000000      ,   0.51076871999999995       };
 Point(         (k-1)*stepk+113 ) = {  -0.12923064000000004      ,    0.0000000000000000      ,   0.55692251999999998       };
 Point(         (k-1)*stepk+114 ) = {  -0.11692296000000003      ,    0.0000000000000000      ,   0.60307631999999989       };
 Point(         (k-1)*stepk+115 ) = {  -0.10461528000000005      ,    0.0000000000000000      ,   0.64923012000000002       };
 Point(         (k-1)*stepk+116 ) = {  -8.92306800000000344E-002 ,    0.0000000000000000      ,   0.69538392000000004       };
 Point(         (k-1)*stepk+117 ) = {  -8.30768400000000407E-002 ,    0.0000000000000000      ,   0.74153771999999996       };
 Point(         (k-1)*stepk+118 ) = {  -7.38460800000000361E-002 ,    0.0000000000000000      ,   0.78769151999999998       };
 Point(         (k-1)*stepk+119 ) = {  -6.76922400000000563E-002 ,    0.0000000000000000      ,   0.83384532000000000       };
 Point(         (k-1)*stepk+120 ) = {  -6.15384000000000556E-002 ,    0.0000000000000000      ,   0.87999911999999991       };
 Point(         (k-1)*stepk+121 ) = {  -5.53845600000000549E-002 ,    0.0000000000000000      ,   0.92615291999999994       };
 Point(         (k-1)*stepk+2*stepj+2 ) = {  -5.23076400000000649E-002 ,    0.0000000000000000      ,   0.97230672000000007       };
 Point(         (k-1)*stepk+123 ) = {   4.30768799999999355E-002 ,    0.0000000000000000      ,   0.98461440000000000       };
 Point(         (k-1)*stepk+2*stepj+3 ) = {   6.46153199999999345E-002 ,    0.0000000000000000      ,   0.93846059999999998       };
 Point(         (k-1)*stepk+125 ) = {   8.61537599999999404E-002 ,    0.0000000000000000      ,   0.89230680000000007       };
 Point(         (k-1)*stepk+126 ) = {   0.14153831999999994      ,    0.0000000000000000      ,   0.84615300000000004       };
 Point(         (k-1)*stepk+127 ) = {   0.21538439999999995      ,    0.0000000000000000      ,   0.79999920000000002       };
 Point(         (k-1)*stepk+128 ) = {   0.24923051999999996      ,    0.0000000000000000      ,   0.75384540000000000       };
 Point(         (k-1)*stepk+129 ) = {   0.27384587999999999      ,    0.0000000000000000      ,   0.70769159999999998       };
 Point(         (k-1)*stepk+130 ) = {   0.28307663999999993      ,    0.0000000000000000      ,   0.66153779999999995       };
 Point(         (k-1)*stepk+131 ) = {   0.29230740000000000      ,    0.0000000000000000      ,   0.61538400000000004       };
 Point(         (k-1)*stepk+132 ) = {   0.29538431999999992      ,    0.0000000000000000      ,   0.56923020000000002       };
 Point(         (k-1)*stepk+133 ) = {   0.29846123999999996      ,    0.0000000000000000      ,   0.52307640000000000       };
 Point(         (k-1)*stepk+134 ) = {   0.29538431999999992      ,    0.0000000000000000      ,   0.47692259999999997       };
 Point(         (k-1)*stepk+135 ) = {   0.28923047999999996      ,    0.0000000000000000      ,   0.43076880000000001       };
 Point(         (k-1)*stepk+136 ) = {   0.27692279999999997      ,    0.0000000000000000      ,   0.38461499999999998       };
 Point(         (k-1)*stepk+137 ) = {   0.26153819999999994      ,    0.0000000000000000      ,   0.33846120000000002       };
 Point(         (k-1)*stepk+138 ) = {   0.25538435999999998      ,    0.0000000000000000      ,   0.29230740000000000       };
 Point(         (k-1)*stepk+139 ) = {   0.23999976000000001      ,    0.0000000000000000      ,   0.24615360000000003       };
 Point(         (k-1)*stepk+140 ) = {   0.21846131999999999      ,    0.0000000000000000      ,   0.19999980000000001       };
 Point(         (k-1)*stepk+124 ) = {   0.19076904000000000      ,    0.0000000000000000      ,   0.15384600000000001       };
 Point(         (k-1)*stepk+122 ) = {   0.15999984000000000      ,    0.0000000000000000      ,   0.10769220000000002       };
 Point(         (k-1)*stepk+1*stepj+3 ) = {   0.11692295999999999      ,    0.0000000000000000      ,   6.15384000000000070E-002  };
 Point(         (k-1)*stepk+102 ) = {   4.92307199999999986E-002 ,    0.0000000000000000      ,   1.53846000000000035E-002  };

// move centre origin
Translate {-OriginX, -OriginY, -OriginZ} {
  Point{(k-1)*stepk+101, (k-1)*stepk+1*stepj+2,   (k-1)*stepk+103, (k-1)*stepk+104, (k-1)*stepk+105, (k-1)*stepk+106, (k-1)*stepk+107, (k-1)*stepk+108, 
        (k-1)*stepk+109, (k-1)*stepk+110, (k-1)*stepk+111, (k-1)*stepk+112, (k-1)*stepk+113, (k-1)*stepk+114, (k-1)*stepk+115, (k-1)*stepk+116, 
        (k-1)*stepk+117, (k-1)*stepk+118, (k-1)*stepk+119, (k-1)*stepk+120, (k-1)*stepk+121, (k-1)*stepk+2*stepj+2  , (k-1)*stepk+123, (k-1)*stepk+2*stepj+3  , 
        (k-1)*stepk+125, (k-1)*stepk+126, (k-1)*stepk+127, (k-1)*stepk+128, (k-1)*stepk+129, (k-1)*stepk+130, (k-1)*stepk+131, (k-1)*stepk+132, 
        (k-1)*stepk+133, (k-1)*stepk+134, (k-1)*stepk+135, (k-1)*stepk+136, (k-1)*stepk+137, (k-1)*stepk+138, (k-1)*stepk+139, (k-1)*stepk+140, 
        (k-1)*stepk+124, (k-1)*stepk+122, (k-1)*stepk+1*stepj+3,   (k-1)*stepk+102 };
}

//Spline(1) = {1, 123, 2};
//Line(1) = {1, 2};
//Spline(2) = {4, 101, 102, 3};
////Line(2) = {4, 3};
//Spline(3) = {4, 104, 106, 108, 110, 112, 114, 116, 118, 120, 121, 1};
//Spline(4) = {3, 124, 139, 137, 135, 133, 131, 129, 128, 127, 126, 125, 2};

//#######################

// fluid domain
x1[] = Point{(k-1)*stepk+2*stepj+2};
x2[] = Point{(k-1)*stepk+2*stepj+3};
x3[] = Point{(k-1)*stepk+1*stepj+2};
x4[] = Point{(k-1)*stepk+1*stepj+3};

r1 = Sqrt( (centreZ - x1[2])^2 + (centreX - x1[0])^2 );
r2 = Sqrt( (centreZ - x2[2])^2 + (centreX - x2[0])^2 );
r3 = Sqrt( (centreZ - x3[2])^2 + (centreX - x3[0])^2 );
r4 = Sqrt( (centreZ - x4[2])^2 + (centreX - x4[0])^2 );

phi1 = Fabs(Atan( (centreX - x1[0]) / (centreZ - x1[2]) ));
phi2 = Fabs(Atan( (centreX - x2[0]) / (centreZ - x2[2]) ));
phi3 = Fabs(Atan( (centreX - x3[0]) / (centreZ - x3[2]) ));
phi4 = Fabs(Atan( (centreX - x4[0]) / (centreZ - x4[2]) ));

r0 = r1 - innerSide; // internal circle
r5 = (r3 + r4)/2 + outerSide; // external circle

// points on outer ellipse
j = 1;
//Point((k-1)*stepk+(j-1)*stepj+1) = {centreX - r5/Cos(alpha/180*Pi),   0, centreZ};
Point((k-1)*stepk+(j-1)*stepj+1) = {centreX - r5,   0, centreZ};
// ellipse major and secondary axes
//a_ellipse = r5/Cos(alpha/180*Pi);
a_ellipse = r5;
b_ellipse = r5;
// compute point on ellipse
x_coord = centreX - r5*Sin(2*phi4);
Call ComputeEllipsePoint;
Point((k-1)*stepk+(j-1)*stepj+2) = {x_coord,    0, z_coord };
//Point((k-1)*stepk+(j-1)*stepj+2) = {centreX - r5*Sin(phi4),    0, centreZ - r5*Cos(phi4) };
// compute point on ellipse
x_coord = centreX + r5*Sin(2*phi3);
Call ComputeEllipsePoint;
Point((k-1)*stepk+(j-1)*stepj+3) = {x_coord,    0, z_coord };
//Point((k-1)*stepk+(j-1)*stepj+3) = {centreX + r5*Sin(phi3),    0, centreZ - r5*Cos(phi3) };
//Point((k-1)*stepk+(j-1)*stepj+4) = {centreX,   0, centreZ + r5};
//Point((k-1)*stepk+(j-1)*stepj+4) = {centreX + r5/Cos(alpha/180*Pi),   0, centreZ};
Point((k-1)*stepk+(j-1)*stepj+4) = {centreX + r5,   0, centreZ};
Point((k-1)*stepk+(j-1)*stepj+5) = {centreX, 0, centreZ + r5};

// points on outer edge 
j = 2;
Point((k-1)*stepk+(j-1)*stepj+1) = {centreX - r4/Cos(alpha/180*Pi),   0, centreZ};
//Point((k-1)*stepk+(j-1)*stepj+1) = {centreX - r4,   0, centreZ};
//Point((k-1)*stepk+(j-1)*stepj+4) = {centreX,   0, centreZ + r3};
Point((k-1)*stepk+(j-1)*stepj+4) = {centreX + r3/Cos(alpha/180*Pi),   0, centreZ};
//Point((k-1)*stepk+(j-1)*stepj+4) = {centreX + r3,   0, centreZ};
Point((k-1)*stepk+(j-1)*stepj+5) = {centreX,     0, centreZ + r3};

// points on inner edge 
j = 3;
Point((k-1)*stepk+(j-1)*stepj+1) = {centreX - r1/Cos(alpha/180*Pi),   0, centreZ};
//Point((k-1)*stepk+(j-1)*stepj+4) = {centreX,   0, centreZ + r2};
Point((k-1)*stepk+(j-1)*stepj+4) = {centreX + r2/Cos(alpha/180*Pi),   0, centreZ};
Point((k-1)*stepk+(j-1)*stepj+5) = {centreX,     0, centreZ + r2};

// points on inner ellipse
j = 4;
// ellipse major and secondary axes
a_ellipse = r0/Cos(alpha/180*Pi);
b_ellipse = r0;
Point((k-1)*stepk+(j-1)*stepj+1) = {centreX - r0/Cos(alpha/180*Pi),   0, centreZ};
// compute point on ellipse
x_coord = centreX - r0*Sin(phi1);
Call ComputeEllipsePoint;
Point((k-1)*stepk+(j-1)*stepj+2) = {x_coord,    0, z_coord };
//Point((k-1)*stepk+(j-1)*stepj+2) = {centreX - r0*Sin(phi1),    0, centreZ - r0*Cos(phi1) };
// compute point on ellipse
x_coord = centreX + r0*Sin(phi2);
Call ComputeEllipsePoint;
Point((k-1)*stepk+(j-1)*stepj+3) = { x_coord,    0, z_coord };
//Point((k-1)*stepk+(j-1)*stepj+3) = {centreX + r0*Sin(phi2),    0, centreZ - r0*Cos(phi2) };
//Point((k-1)*stepk+(j-1)*stepj+4) = {centreX,   0, centreZ + r0};
Point((k-1)*stepk+(j-1)*stepj+4) = {centreX + r0/Cos(alpha/180*Pi),   0, centreZ};
Point((k-1)*stepk+(j-1)*stepj+5) = {centreX,     0, centreZ + r0};

// centre of rotation
Point((k-1)*stepk+100) = {centreX, 0, centreZ};

// Rotate points
k = 2;
//For j In {1:npointj}
//   For i In {1:npointi}
//      pointIndex = (k-1)*stepk + (j-1)*stepj + i;
//      Call RotatePoint;
//   EndFor
//EndFor
For j In {2:npointj}
   For i In {1:npointi-1}
      pointIndex = (k-1)*stepk + (j-1)*stepj + i;
      Call RotatePoint;
   EndFor
EndFor
For p In {100:140}
   pointIndex = (k-1)*stepk + p;
   Call RotatePoint;
EndFor

norm1 = Sin(alpha/180*Pi); norm2 = Cos(alpha/180*Pi); norm3 = 0;
dir1  = norm1; dir2 = norm2; dir3 = norm3;
xp = thickness*norm1; yp = thickness*norm2; zp = 0;

// project points on planform orthogonally
k = 2;
For p In {101:140}
    oldnumber = (k-1)*stepk + p;
    newnumber = oldnumber + stepk;
    Call ProjectPoint;
EndFor
k = 2;
For j In {2:3}
   For i In {2:3}
      oldnumber = (k-1)*stepk + (j-1)*stepj + i;
      newnumber = oldnumber + stepk;
      Call ProjectPoint;
   EndFor
EndFor
dir1  = 0; dir2 = 1; dir3 = 0;
For i In {1:1}
   For j In {2:3}
      oldnumber = (k-1)*stepk + (j-1)*stepj + i;
      newnumber = oldnumber + stepk;
      Call ProjectPoint;
   EndFor
EndFor
For i In {4:npointi}
   For j In {2:3}
      oldnumber = (k-1)*stepk + (j-1)*stepj + i;
      newnumber = oldnumber + stepk;
      Call ProjectPoint;
   EndFor
EndFor

For j In {npointj:npointj}
   For i In {1:npointi}
      oldnumber = (k-1)*stepk + (j-1)*stepj + i;
      newnumber = oldnumber + stepk;
      Call ProjectPoint;
   EndFor
EndFor
// project centre in skewed direction
p = 100;
oldnumber = (k-1)*stepk + p;
newnumber = oldnumber + stepk;
Call ProjectPoint;

// project outer points in skewed direction
dir1  = 0; dir2 = 1; dir3 = 0;
norm1 = 0; norm2 = 1; norm3 = 0;
k = 2;
For j In {1:1}
   For i In {1:npointi}
      oldnumber = (k-1)*stepk + (j-1)*stepj + i;
      newnumber = oldnumber + stepk;
      Call ProjectPoint;
   EndFor
EndFor

// upper surface from layer 3
k = 3;
norm1 = 0; norm2 = 1; norm3 = 0;
dir1  = norm1; dir2 = norm2; dir3 = norm3;
xp = 0; yp = top; zp = 0;

For p In {101:140}
    oldnumber = (k-1)*stepk + p;
    newnumber = oldnumber + stepk;
    Call ProjectPoint;
EndFor
For i In {1:npointi}
   For j In {1:npointj}
      oldnumber = (k-1)*stepk + (j-1)*stepj + i;
      newnumber = oldnumber + stepk;
      Call ProjectPoint;
   EndFor
EndFor
// project centre in skewed direction
p = 100;
oldnumber = (k-1)*stepk + p;
newnumber = oldnumber + stepk;
Call ProjectPoint;

// lower surface from layer 2
k = 2;
norm1 = 0; norm2 = 1; norm3 = 0;
dir1  = -norm1; dir2 = -norm2; dir3 = -norm3;
xp = 0; yp = -bottom; zp = 0;

For p In {101:140}
    oldnumber = (k-1)*stepk + p;
    newnumber = oldnumber - stepk;
    Call ProjectPoint;
EndFor
For i In {1:npointi}
   For j In {1:npointj}
      oldnumber = (k-1)*stepk + (j-1)*stepj + i;
      newnumber = oldnumber - stepk;
      Call ProjectPoint;
   EndFor
EndFor
// project centre in skewed direction
p = 100;
oldnumber = (k-1)*stepk + p;
newnumber = oldnumber - stepk;
Call ProjectPoint;

lineShiftj = 50;
lineShiftk = 100;

// Generate lines
For k In {1:4}

   // lines front to back
   For i In {1:npointi-1}
      For j In {1:npointj}
         // in x direction
         If ( !(i == 2 && j == 2 ) && !(i == 2 && j == 3) )
            Ellipse((k-1)*stepk + (j-1)*stepj + i)  = { (k-1)*stepk + (j-1)*stepj + i,  (k-1)*stepk + 100, (k-1)*stepk + 1, (k-1)*stepk + (j-1)*stepj + i+1 };
         EndIf
      EndFor
   EndFor
   If (EncloseCircle == 1) 
      // enclose circle
      i = npointi;
      For j In {1:npointj}
         // in x direction
         Ellipse((k-1)*stepk + (j-1)*stepj + i)  = { (k-1)*stepk + (j-1)*stepj + i,  (k-1)*stepk + 100, (k-1)*stepk + 1, (k-1)*stepk + (j-1)*stepj + 1 };
      EndFor
   EndIf

   i = 2; j = 2;
   Spline((k-1)*stepk + (j-1)*stepj + i) = {(k-1)*stepk + (j-1)*stepj + i, (k-1)*stepk + 101, (k-1)*stepk + 102, (k-1)*stepk + (j-1)*stepj + i+1};
   i = 2; j = 3;
   Line((k-1)*stepk + (j-1)*stepj + i) = {(k-1)*stepk + (j-1)*stepj + i, (k-1)*stepk + (j-1)*stepj + i+1};

   // lines right to left
   For i In {1:npointi}
      For j In {1:npointj-1}
         // in z direction
         If ( !(i == 2 && j == 2 ) && !(i == 3 && j == 2) )
            Line(lineShiftj + (k-1)*stepk + (j-1)*stepj + i)  = { (k-1)*stepk + (j-1)*stepj + i,  (k-1)*stepk + (j)*stepj + i };
         EndIf
      EndFor
   EndFor

   i = 2; j = 2;
   Spline(lineShiftj + (k-1)*stepk + (j-1)*stepj + i) = {(k-1)*stepk + (j-1)*stepj + i, (k-1)*stepk + 104, (k-1)*stepk + 106, (k-1)*stepk + 108, (k-1)*stepk + 110, (k-1)*stepk + 112, (k-1)*stepk + 114, (k-1)*stepk + 116, (k-1)*stepk + 118, (k-1)*stepk + 120, (k-1)*stepk + 121, (k-1)*stepk + (j)*stepj + i};
   i = 3; j = 2;
   Spline(lineShiftj + (k-1)*stepk + (j-1)*stepj + i) = {(k-1)*stepk + (j-1)*stepj + i, (k-1)*stepk + 124, (k-1)*stepk + 139, (k-1)*stepk + 137, (k-1)*stepk + 135, (k-1)*stepk + 133, (k-1)*stepk + 131, (k-1)*stepk + 129, (k-1)*stepk + 128, (k-1)*stepk + 127, (k-1)*stepk + 126, (k-1)*stepk + 125,  (k-1)*stepk + (j)*stepj + i};

   // lines bottom to up
   For i In {1:npointi}
      For j In {1:npointj}
         // in y direction
         If ( k != 4 )
            Line(lineShiftk + (k-1)*stepk + (j-1)*stepj + i)  = { (k-1)*stepk + (j-1)*stepj + i,  (k)*stepk + (j-1)*stepj + i };
         EndIf
      EndFor
   EndFor

   // division front to back
   For i In {1:1}
      For j In {1:npointj}
         Transfinite Line {(k-1)*stepk + (j-1)*stepj + i} = numEl_front Using Progression progression_front;
      EndFor
   EndFor
   For i In {2:2}
      For j In {1:npointj}
         Transfinite Line {(k-1)*stepk + (j-1)*stepj + i} = numEl_chord Using Bump progression_chord;
      EndFor
   EndFor
   For i In {3:3}
      For j In {1:npointj}
         Transfinite Line {(k-1)*stepk + (j-1)*stepj + i} = numEl_back Using Progression progression_back;
      EndFor
   EndFor
   For i In {4:4}
      For j In {1:npointj}
         Transfinite Line {(k-1)*stepk + (j-1)*stepj + i} = numEl_quadrant Using Progression progression_quadrant;
      EndFor
   EndFor
   If (EncloseCircle == 1) 
      For i In {5:5}
         For j In {1:npointj}
            Transfinite Line {(k-1)*stepk + (j-1)*stepj + i} = numEl_quadrant;
         EndFor
      EndFor
   EndIf

   // division right to left
   For i In {1:npointi}
      For j In {1:1}
         Transfinite Line {lineShiftj + (k-1)*stepk + (j-1)*stepj + i} = numEl_outerSide Using Progression progression_outer;
      EndFor
   EndFor
   For i In {1:npointi}
      For j In {2:2}
         Transfinite Line {lineShiftj + (k-1)*stepk + (j-1)*stepj + i} = numEl_span Using Progression progression_span;
      EndFor
   EndFor
   For i In {1:npointi}
      For j In {3:3}
         Transfinite Line {lineShiftj + (k-1)*stepk + (j-1)*stepj + i} = numEl_innerSide Using Bump progression_inner;
      EndFor
   EndFor
EndFor

// division botom to up
k = 1;
For i In {1:npointi}
   For j In {1:npointj}
      Transfinite Line {lineShiftk + (k-1)*stepk + (j-1)*stepj + i} = numEl_bottom Using Progression progression_bottom;
   EndFor
EndFor
k = 2;
For i In {1:npointi}
   For j In {1:npointj}
      Transfinite Line {lineShiftk + (k-1)*stepk + (j-1)*stepj + i} = numEl_thickness;
   EndFor
EndFor
k = 3;
For i In {1:npointi}
   For j In {1:npointj}
      Transfinite Line {lineShiftk + (k-1)*stepk + (j-1)*stepj + i} = numEl_top Using Progression progression_top;
   EndFor
EndFor

lineLoopShift = 10000;

// Generate horizontal surfaces
For k In {1:npointk}
   For j In {1:npointj-1}
      For i In {1:npointi-1}
         lineLoopIndex = lineLoopShift + (k-1)*stepk + (j-1)*stepj + i;
         Line Loop(lineLoopIndex) = { lineShiftj + (k-1)*stepk + (j-1)*stepj + i, 
                                                   (k-1)*stepk + (j  )*stepj + i,
                                    -(lineShiftj + (k-1)*stepk + (j-1)*stepj + i+1),
                                                  -((k-1)*stepk + (j-1)*stepj + i) };
         newSurfIndex = (k-1)*stepk + (j-1)*stepj + i;
         Ruled Surface(newSurfIndex) = {lineLoopIndex};
         Transfinite Surface {newSurfIndex} = { (k-1)*stepk + (j  )*stepj + i, 
                                                (k-1)*stepk + (j  )*stepj + i+1,
                                                (k-1)*stepk + (j-1)*stepj + i+1,
                                                (k-1)*stepk + (j-1)*stepj + i };
         Recombine Surface {newSurfIndex};
      EndFor
      If (EncloseCircle == 1) 
         // enclosing surface
         i = npointi;
         lineLoopIndex = lineLoopShift + (k-1)*stepk + (j-1)*stepj + i;
         Line Loop(lineLoopIndex) = { lineShiftj + (k-1)*stepk + (j-1)*stepj + i, 
                                                   (k-1)*stepk + (j  )*stepj + i,
                                    -(lineShiftj + (k-1)*stepk + (j-1)*stepj + 1),
                                                  -((k-1)*stepk + (j-1)*stepj + i) };
         newSurfIndex = (k-1)*stepk + (j-1)*stepj + i;
         Ruled Surface(newSurfIndex) = {lineLoopIndex};
         Transfinite Surface {newSurfIndex} = { (k-1)*stepk + (j  )*stepj + i, 
                                                (k-1)*stepk + (j  )*stepj + 1,
                                                (k-1)*stepk + (j-1)*stepj + 1,
                                                (k-1)*stepk + (j-1)*stepj + i };
         Recombine Surface {newSurfIndex};
      EndIf
   EndFor
EndFor

// Generate vertical surfaces in radial direction
For k In {1:npointk-1}
   For j In {1:npointj-1}
      For i In {1:npointi}
         lineLoopIndex = lineLoopShift + lineShiftj + (k-1)*stepk + (j-1)*stepj + i;
         Line Loop(lineLoopIndex) = { lineShiftj + (k-1)*stepk + (j-1)*stepj + i, 
                                      lineShiftk + (k-1)*stepk + (j  )*stepj + i,
                                    -(lineShiftj + (k  )*stepk + (j-1)*stepj + i),
                                    -(lineShiftk + (k-1)*stepk + (j-1)*stepj + i) };
         newSurfIndex = lineShiftj + (k-1)*stepk + (j-1)*stepj + i;
         Ruled Surface(newSurfIndex) = {lineLoopIndex};
         Transfinite Surface {newSurfIndex} = { (k-1)*stepk + (j-1)*stepj + i, 
                                                (k-1)*stepk + (j  )*stepj + i,
                                                (k  )*stepk + (j  )*stepj + i,
                                                (k  )*stepk + (j-1)*stepj + i };
         Recombine Surface {newSurfIndex};
      EndFor
   EndFor
EndFor

// Generate vertical surfaces in tangential direction
For k In {1:npointk-1}
   For j In {1:npointj}
      For i In {1:npointi-1}
         lineLoopIndex = lineLoopShift + lineShiftk + (k-1)*stepk + (j-1)*stepj + i;
         Line Loop(lineLoopIndex) = {              (k-1)*stepk + (j-1)*stepj + i, 
                                      lineShiftk + (k-1)*stepk + (j-1)*stepj + i+1,
                                    -(             (k  )*stepk + (j-1)*stepj + i),
                                    -(lineShiftk + (k-1)*stepk + (j-1)*stepj + i) };
         newSurfIndex = lineShiftk + (k-1)*stepk + (j-1)*stepj + i;
         Ruled Surface(newSurfIndex) = {lineLoopIndex};
         Transfinite Surface {newSurfIndex} = { (k-1)*stepk + (j-1)*stepj + i, 
                                                (k-1)*stepk + (j-1)*stepj + i+1,
                                                (k  )*stepk + (j-1)*stepj + i+1,
                                                (k  )*stepk + (j-1)*stepj + i };
         Recombine Surface {newSurfIndex};
      EndFor
      
      If (EncloseCircle == 1) 
         // enclose circle
         i = npointi;
         lineLoopIndex = lineLoopShift + lineShiftk + (k-1)*stepk + (j-1)*stepj + i;
         Line Loop(lineLoopIndex) = {              (k-1)*stepk + (j-1)*stepj + i, 
                                      lineShiftk + (k-1)*stepk + (j-1)*stepj + 1,
                                    -(             (k  )*stepk + (j-1)*stepj + i),
                                    -(lineShiftk + (k-1)*stepk + (j-1)*stepj + i) };
         newSurfIndex = lineShiftk + (k-1)*stepk + (j-1)*stepj + i;
         Ruled Surface(newSurfIndex) = {lineLoopIndex};
         Transfinite Surface {newSurfIndex} = { (k-1)*stepk + (j-1)*stepj + i, 
                                                (k-1)*stepk + (j-1)*stepj + 1,
                                                (k  )*stepk + (j-1)*stepj + 1,
                                                (k  )*stepk + (j-1)*stepj + i };
         Recombine Surface {newSurfIndex};
      EndIf
   EndFor
EndFor

// Generate volumes
surfaceLoopShift = lineLoopShift;
For k In {1:npointk-1}
   For j In {1:npointj-1}
      For i In {1:npointi-1}
         // skip plate itself
         If ( !(i == 2 && j == 2 && k == 2) )
            surfaceLoopIndex = surfaceLoopShift + (k-1)*stepk + (j-1)*stepj + i;
            Surface Loop(surfaceLoopIndex) = { (k-1)*stepk + (j-1)*stepj + i, 
                                              -((k  )*stepk + (j-1)*stepj + i),
                                                lineShiftk + (k-1)*stepk + (j-1)*stepj + i,
                                              -(lineShiftk + (k-1)*stepk + (j  )*stepj + i),
                                                lineShiftj + (k-1)*stepk + (j-1)*stepj + i,
                                              -(lineShiftj + (k-1)*stepk + (j-1)*stepj + i+1) };
            newVolumeIndex = (k-1)*stepk + (j-1)*stepj + i;
            Volume(newVolumeIndex) = {surfaceLoopIndex};
            Transfinite Volume{newVolumeIndex};
         EndIf
      EndFor
      If (EncloseCircle == 1) 
         // enclose circle
         i = npointi;
         surfaceLoopIndex = surfaceLoopShift + (k-1)*stepk + (j-1)*stepj + i;
         Surface Loop(surfaceLoopIndex) = { (k-1)*stepk + (j-1)*stepj + i, 
                                          -((k  )*stepk + (j-1)*stepj + i),
                                            lineShiftk + (k-1)*stepk + (j-1)*stepj + i,
                                          -(lineShiftk + (k-1)*stepk + (j  )*stepj + i),
                                            lineShiftj + (k-1)*stepk + (j-1)*stepj + i,
                                          -(lineShiftj + (k-1)*stepk + (j-1)*stepj + 1) };
         newVolumeIndex = (k-1)*stepk + (j-1)*stepj + i;
         Volume(newVolumeIndex) = {surfaceLoopIndex};
         Transfinite Volume{newVolumeIndex};
      EndIf
   EndFor
EndFor

// outer boundary
inflowList[] = {};
// lower surface
k = 1;
For j In {1:npointj-1}
   For i In {1:npointi-1}
      inflowList[] = {inflowList[], (k-1)*stepk + (j-1)*stepj + i };
   EndFor
   If (EncloseCircle == 1) 
      i = npointi;
      inflowList[] = {inflowList[], (k-1)*stepk + (j-1)*stepj + i };
   EndIf
EndFor

// upper surface
k = 4;
For j In {1:npointj-1}
   For i In {1:npointi-1}
      inflowList[] = {inflowList[], (k-1)*stepk + (j-1)*stepj + i };
   EndFor
   If (EncloseCircle == 1) 
      i = npointi;
      inflowList[] = {inflowList[], (k-1)*stepk + (j-1)*stepj + i };
   EndIf
EndFor
For k In {1:npointk-1}
   // outer circle
   j = 1;
   For i In {1:npointi-1}
      inflowList[] = {inflowList[], lineShiftk + (k-1)*stepk + (j-1)*stepj + i };
   EndFor
   If (EncloseCircle == 1) 
      i = npointi;
      inflowList[] = {inflowList[], lineShiftk + (k-1)*stepk + (j-1)*stepj + i };
   EndIf
   // inner circle
   j = 4;
   For i In {1:npointi-1}
      inflowList[] = {inflowList[], lineShiftk + (k-1)*stepk + (j-1)*stepj + i };
   EndFor
   If (EncloseCircle == 1) 
      i = npointi;
      inflowList[] = {inflowList[], lineShiftk + (k-1)*stepk + (j-1)*stepj + i };
   EndIf
   // Add inflow and outflow cuts
   If (EncloseCircle == 0) 
      i = 1;
      For j In {1:npointj-1}
         inflowList[] = {inflowList[], lineShiftj + (k-1)*stepk + (j-1)*stepj + i };
      EndFor
      i = npointi;
      For j In {1:npointj-1}
         inflowList[] = {inflowList[], lineShiftj + (k-1)*stepk + (j-1)*stepj + i };
      EndFor
   EndIf
EndFor
Physical Surface("boundary") = { inflowList[] };
//Printf("New members '%g','%g'", inflowList[1],inflowList[20] ); 

// plate
plateList[] = {};
//horizontal faces
i = 2; j = 2; 
For k In {2:3}
   plateList[] = { plateList[], (k-1)*stepk + (j-1)*stepj + i };
EndFor
//vertical faces
k = 2; 
j = 2;
For i In {2:3}
   plateList[] = { plateList[], lineShiftj + (k-1)*stepk + (j-1)*stepj + i };
EndFor
k = 2; 
i = 2;
For j In {2:3}
   plateList[] = { plateList[], lineShiftk + (k-1)*stepk + (j-1)*stepj + i };
EndFor

Physical Surface("plate") = { plateList[] };

// domain
domainList[] = {};
For k In {1:npointk-1}
   For j In {1:npointj-1}
      For i In {1:npointi-1}
         // skip plate itself
         If ( !(i == 2 && j == 2 && k == 2) )
            volumeIndex = (k-1)*stepk + (j-1)*stepj + i;
            domainList[] = { domainList[], volumeIndex };
         EndIf
      EndFor
      If (EncloseCircle == 1) 
         i = npointi;
         // skip plate itself
         If ( !(i == 2 && j == 2 && k == 2) )
            volumeIndex = (k-1)*stepk + (j-1)*stepj + i;
            domainList[] = { domainList[], volumeIndex };
         EndIf
      EndIf
   EndFor
EndFor
Physical Volume("domain") = {domainList[]};
//Physical Volume("domain") = {20, 22, 25, 19, 26, 24, 18, 21, 23, 1, 4, 6, 2, 7, 3, 5, 8, 9, 12, 14, 10, 17, 15, 11, 13, 16};
