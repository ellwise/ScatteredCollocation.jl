// -*- C++ -*- 

c = 1;
t = 0.01875 * c;
w = 2*c;

numE_t = 10;
numE_c = 100;
numE_w = 80;

// parameters in xy plane 
l1 = 10*c;
l2 = 20*c;
l5 = 10*c;
l7 = 3*c;
l8 = l7;
l6 = l5;

progression_l1 = 1.10; 
progression_l2 = 1.09; 
progression_l7 = 1.1; 
progression_l6_l8 = 1.2; 
progression_l8 = progression_l7; 
progression_l5_l7 = progression_l6_l8;

numDiv_c = numE_c;
numDiv_t = numE_t;
numDiv_l7 = 40;
numDiv_l2 = 60;
numDiv_l1 = 50;
numDiv_l6_l8 = 10;
numDiv_l8 = numDiv_l7;
numDiv_l5_l7 = numDiv_l6_l8;

// geometry in xy plane
Point(32) = {0,0,0};
Point(1) = {c,0,0};
Point(2) = {c,t,0};
Point(3) = {0,t,0};

Point(4) = {-l1,-l5,0};
Point(6) = {0,-l5,0};
Point(7) = {c,-l5,0};
Point(9) = {c+l2,-l5,0};

Point(10) = {-l1,-l7,0};
Point(12) = {0,-l7,0};
Point(13) = {c,-l7,0};
Point(15) = {c+l2,-l7,0};

Point(33) = {-l1,0,0};
Point(36) = {c+l2,0,0};

Point(16) = {-l1,t,0};
Point(19) = {c+l2,t,0};

Point(20) = {-l1,t+l8,0};
Point(22) = {0,t+l8,0};
Point(23) = {c,t+l8,0};
Point(25) = {c+l2,t+l8,0};

Point(26) = {-l1,t+l6,0};
Point(28) = {0,t+l6,0};
Point(29) = {c,t+l6,0};
Point(31) = {c+l2,t+l6,0};

Line(2) = {4, 6};
Line(4) = {6, 7};
Line(6) = {7, 9};
Line(7) = {10, 12};
Line(9) = {12, 13};
Line(11) = {13, 15};
Line(12) = {33, 32};
Line(16) = {1, 36};
Line(17) = {16, 3};
Line(21) = {2, 19};
Line(22) = {26, 28};
Line(24) = {28, 29};
Line(26) = {29, 31};
Line(27) = {20, 22};
Line(29) = {22, 23};
Line(31) = {23, 25};
Line(33) = {32, 1};
Line(35) = {3, 2};

Line(37) = {4, 10};
Line(38) = {10, 33};
Line(39) = {33, 16};
Line(40) = {16, 20};
Line(41) = {20, 26};
Line(47) = {6, 12};
Line(48) = {12, 32};
Line(49) = {32, 3};
Line(50) = {3, 22};
Line(51) = {22, 28};
Line(52) = {7, 13};
Line(53) = {13, 1};
Line(54) = {1, 2};
Line(55) = {2, 23};
Line(56) = {23, 29};
Line(62) = {9, 15};
Line(63) = {15, 36};
Line(64) = {36, 19};
Line(65) = {19, 25};
Line(66) = {25, 31};

Transfinite Line {39,49,54,64} = numDiv_t;
Transfinite Line {4,9,33,35,29,24} = numDiv_c;

Transfinite Line {40,50,55,65} = numDiv_l8 Using Progression progression_l8;
Transfinite Line {38,48,53,63} = numDiv_l7 Using Progression 1/progression_l7;

Transfinite Line {2,7,12,17,27,22} = numDiv_l1 Using Progression 1/progression_l1;
Transfinite Line {6, 11, 16, 21, 31, 26} = numDiv_l2 Using Progression progression_l2;

Transfinite Line {41,51,56,66} = numDiv_l6_l8 Using Progression progression_l6_l8;
Transfinite Line {37,47,52,62} = numDiv_l5_l7 Using Progression 1/progression_l5_l7;



Line Loop(67) = {2,47,-7,-37};
Ruled Surface(1) = {67};
Transfinite Surface {1} = {4,6,12,10};

Line Loop(69) = {4,52,-9,-47};
Ruled Surface(3) = {69};
Transfinite Surface {3} = {6,7,13,12};

Line Loop(71) = {6,62,-11,-52};
Ruled Surface(5) = {71};
Transfinite Surface {5} = {7,9,15,13};

Line Loop(72) = {7,48,-12,-38};
Ruled Surface(6) = {72};
Transfinite Surface {6} = {10,12,32,33};

Line Loop(74) = {9,53,-33,-48};
Ruled Surface(8) = {74};
Transfinite Surface {8} = {12,13,1,32};

Line Loop(76) = {11,63,-16,-53};
Ruled Surface(10) = {76};
Transfinite Surface {10} = {13,15,36,1};

Line Loop(77) = {12,49,-17,-39};
Ruled Surface(11) = {77};
Transfinite Surface {11} = {33,32,3,16};

Line Loop(79) = {33,54,-35,-49};
Ruled Surface(13) = {79};
Transfinite Surface {13} = {32,1,2,3};

Line Loop(81) = {16,64,-21,-54};
Ruled Surface(15) = {81};
Transfinite Surface {15} = {1,36,19,2};

Line Loop(82) = {17,50,-27,-40};
Ruled Surface(16) = {82};
Transfinite Surface {16} = {16,3,22,20};

Line Loop(84) = {35,55,-29,-50};
Ruled Surface(18) = {84};
Transfinite Surface {18} = {3,2,23,22};

Line Loop(86) = {21,65,-31,-55};
Ruled Surface(20) = {86};
Transfinite Surface {20} = {2,19,25,23};

Line Loop(87) = {27,51,-22,-41};
Ruled Surface(21) = {87};
Transfinite Surface {21} = {20,22,28,26};

Line Loop(89) = {29,56,-24,-51};
Ruled Surface(23) = {89};
Transfinite Surface {23} = {22,23,29,28};

Line Loop(91) = {31,66,-26,-56};
Ruled Surface(25) = {91};
Transfinite Surface {25} = {23,25,31,29};

Recombine Surface {1,3,5,6,8,10,11,13,15,16,18,20,21,23,25};

//parameters in yz plane
l9 =  5*c;
l11 = 2*c;

progression_l9_l11 = 1.5;
progression_l11 = 1.2;

numDiv_l11 = 15;
numDiv_l9_l11 = 5;

// geometry in yz plane
Extrude{0,0,w}{
    Surface{1,3,5,6,8,10,11,13,15,16,18,20,21,23,25}; 
    Layers{numE_w}; 
    Recombine;
}

// find layer locations in l11
r1 = progression_l11;
n1 = numDiv_l11;

a1 = (r1 - 1) / (r1^n1 - 1);
one1[0] = 1;
layer1[0] = a1;

For i In {1:n1-1}
one1[i] = 1;
layer1[i] = layer1[i-1] + a1 * r1^i;
EndFor

Extrude{0,0,-l11}{
    Surface{1,3,5,6,8,10,11,13,15,16,18,20,21,23,25}; 
    Layers{one1[],layer1[]}; 
    Recombine;
}

Extrude{0,0,l11}{
    Surface{113, 135, 157, 179, 201, 223, 245, 267, 289, 311, 333, 355, 377, 399, 421}; 
    Layers{one1[],layer1[]}; 
    Recombine;
}


// find layer locations in l9_l11
r2 = progression_l9_l11;
n2 = numDiv_l9_l11;

a2 = (r2 - 1) / (r2^n2 - 1);
one2[0] = 1;
layer2[0] = a2;

For i In {1:n2-1}
one2[i] = 1;
layer2[i] = layer2[i-1] + a2 * r2^i;
EndFor

Extrude{0,0,-(l9-l11)}{
    Surface{443, 465, 487, 509, 531, 553, 575, 597, 619, 641, 663, 685, 707, 729, 751}; 
    Layers{one2[],layer2[]}; 
    Recombine;
}

Extrude{0,0,(l9-l11)}{
    Surface{773, 795, 817, 839, 861, 883, 905, 927, 949, 971, 993, 1015, 1037, 1059, 1081}; 
    Layers{one2[],layer2[]}; 
    Recombine;
}


Physical Volume("Domain") = {1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75};

// Physical Line("surface") = {33,54,35,49};
// Physical Line("boundary") = {2,4,6,62,63,64,65,66,23,24,26,37,38,39,40,41};
// Physical Surface("Domain")= {1,3,5,6,8,10,11,15,16,18,20,21,23,25};
