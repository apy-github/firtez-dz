!
MODULE DAMPING
  !
  ! J M Borrero
  ! July 22, 2013
  ! KIS, Freiburg
  !
  !================================================
  !
  USE CONS_PARAM, ONLY: SP, DP, EVOLT, HPLA, LIGHT &
      , KBOL, DPI, MAMU, RBOHR
  USE ATOM_DATABASE, ONLY: XI, XII, MATOM, ABUND, REFRAX
  USE DERIVVAR
  USE CODE_MODES, ONLY: MRESPFUNCT
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  !  
  !  This is the s-p data from Anstee and O'Mara 1995, MNRAS 276,859
  !  
  REAL(SP), PARAMETER :: CS_SP(18,21) = RESHAPE(SHAPE= (/18,21 /), SOURCE= &
       (/ 126,   140,   165,  202,  247,  299,  346,  383,  435,  491,  553,  617,  685,  769,  838,  925, 1011, 1082,&
       140,   150,   162,  183,  218,  273,  327,  385,  440,  501,  557,  620,  701,  764,  838,  923, 1025, 1085,&
       154,   167,   175,  192,  216,  251,  299,  357,  423,  487,  549,  617,  684,  759,  834,  910, 1014, 1064,&
       166,   180,   192,  206,  226,  253,  291,  339,  397,  459,  532,  600,  676,  755,  832,  896, 1002, 1055,&
       208,   194,   207,  223,  242,  265,  296,  335,  384,  445,  511,  583,  656,  726,  817,  889,  988, 1044,&
       262,   254,   220,  239,  261,  283,  310,  344,  388,  442,  496,  568,  635,  725,  791,  890,  970, 1036,&
       311,   306,   299,  251,  280,  304,  330,  361,  396,  443,  500,  563,  630,  704,  796,  880,  951, 1033,&
       358,   359,   350,  338,  293,  323,  352,  381,  416,  455,  511,  566,  635,  706,  780,  859,  946, 1039,&
       411,   409,   405,  392,  370,  340,  375,  406,  439,  478,  525,  580,  644,  714,  790,  873,  961, 1050,&
       462,   463,   459,  450,  443,  400,  394,  432,  467,  501,  546,  595,  650,  711,  786,  873,  963, 1050,&
       522,   525,   529,  524,  516,  518,  438,  454,  495,  532,  565,  621,  671,  741,  813,  874,  951, 1034,&
       589,   593,   590,  583,  579,  568,  565,  483,  517,  560,  600,  644,  691,  752,  821,  904,  978, 1048,&
       658,   655,   666,  657,  649,  653,  649,  587,  549,  592,  674,  674,  728,  782,  833,  902,  992, 1084,&
       738,   742,   747,  725,  721,  729,  699,  730,  626,  622,  668,  721,  765,  809,  887,  938, 1001, 1109,&
       838,   838,   810,  809,  790,  800,  769,  815,  757,  679,  704,  755,  806,  854,  901,  974, 1034, 1105,&
       942,   946,   925,  901,  918,  895,  919,  897,  933,  890,  785,  797,  859,  908,  976, 1020, 1115, 1173,&
       1059,  1061,  1056, 1061, 1074, 1031, 1036, 1036,  993, 1038,  932,  852,  878,  943, 1003, 1074, 1131, 1200,&
       1069,  1076,  1083, 1095, 1102, 1091, 1126, 1156, 1103, 1149, 1157, 1036,  972, 1007, 1064, 1124, 1209, 1283,&
       1338,  1350,  1356, 1354, 1324, 1301, 1312, 1318, 1257, 1239, 1297, 1233, 1089, 1059, 1106, 1180, 1218, 1317,&
       1409,  1398,  1367, 1336, 1313, 1313, 1409, 1354, 1317, 1287, 1353, 1386, 1279, 1158, 1141, 1188, 1260, 1335,&
       1328,  1332,  1342, 1369, 1405, 1451, 1502, 1524, 1506, 1477, 1522, 1594, 1572, 1436, 1328, 1325, 1382, 1446 /))
  !
  REAL(SP), PARAMETER :: AL_SP(18,21) = RESHAPE(SHAPE= (/18,21 /), SOURCE= &
       (/.268, .269, .335, .377, .327, .286, .273, .270, .271, .268, .267, .264, .264, .264, .261, .256, .248, .245,&
       .261, .256, .254, .282, .327, .355, .321, .293, .287, .271, .267, .273, .270, .270, .268, .268, .264, .263,&
       .266, .264, .257, .252, .267, .289, .325, .339, .319, .301, .292, .284, .281, .281, .277, .282, .276, .274,&
       .262, .274, .258, .251, .247, .254, .273, .291, .316, .322, .320, .302, .294, .290, .287, .292, .283, .277,& 
       .322, .275, .264, .259, .250, .245, .273, .255, .271, .284, .294, .308, .296, .299, .288, .289, .282, .278,&
       .267, .300, .260, .268, .254, .242, .243, .242, .239, .246, .267, .277, .280, .290, .282, .281, .274, .271,&
       .259, .274, .275, .252, .265, .248, .249, .237, .283, .236, .247, .254, .254, .271, .268, .267, .258, .262,&
       .260, .255, .268, .268, .268, .264, .248, .239, .229, .240, .236, .234, .238, .244, .252, .251, .244, .255,&
       .255, .255, .244, .247, .317, .246, .255, .244, .237, .231, .227, .231, .235, .232, .235, .241, .237, .245,&
       .256, .254, .254, .249, .227, .319, .253, .253, .240, .237, .238, .233, .231, .230, .228, .234, .227, .241,&
       .257, .254, .252, .235, .253, .240, .284, .251, .246, .241, .235, .228, .222, .225, .225, .219, .228, .233,&
       .244, .240, .245, .238, .248, .230, .283, .252, .244, .244, .238, .235, .234, .236, .228, .224, .225, .231,&
       .244, .241, .244, .237, .237, .249, .219, .324, .239, .245, .242, .242, .232, .233, .221, .227, .231, .218,&
       .241, .245, .249, .239, .243, .250, .217, .254, .308, .237, .247, .244, .234, .228, .233, .224, .227, .226,&
       .243, .243, .232, .227, .235, .253, .227, .220, .320, .270, .243, .252, .248, .238, .234, .241, .225, .227,&
       .225, .226, .234, .230, .226, .233, .249, .225, .216, .300, .286, .237, .240, .247, .243, .234, .231, .238,&
       .268, .260, .247, .238, .233, .241, .254, .248, .207, .227, .315, .260, .226, .237, .240, .239, .239, .240,&
       .248, .246, .238, .226, .213, .221, .226, .226, .204, .194, .248, .316, .234, .216, .236, .233, .221, .230,&
       .200, .202, .198, .194, .206, .207, .227, .224, .207, .185, .198, .275, .315, .233, .229, .231, .233, .236,&
       .202, .209, .221, .226, .230, .245, .202, .257, .246, .225, .215, .246, .320, .321, .244, .239, .251, .253,& 
       .246, .248, .255, .265, .274, .285, .292, .284, .273, .250, .225, .239, .295, .352, .320, .258, .260, .269/))   
  ! 
  !  p-d data from Barklem and O'Mara 1997, MNRAS, 290, 102
  !
  REAL(SP), PARAMETER :: CS_PD(18,18) =  RESHAPE(SHAPE= (/18,18 /), SOURCE= &
       (/425,  461,  507,  566,  630,  706,  799,  889,  995, 1083, 1191, 1334, 1478, 1608, 1790, 1870, 1936, 2140,&
       429,  460,  505,  565,  633,  704,  795,  896,  985, 1082, 1199, 1340, 1487, 1611, 1795, 1872, 1937, 2136,&
       419,  451,  501,  556,  627,  700,  785,  891,  977, 1088, 1212, 1346, 1493, 1604, 1793, 1863, 1930, 2144,&
       402,  437,  489,  544,  614,  695,  779,  875,  975, 1102, 1221, 1350, 1488, 1591, 1774, 1844, 1919, 2126,&
       384,  418,  467,  529,  595,  674,  769,  856,  976, 1108, 1224, 1338, 1467, 1570, 1743, 1817, 1900, 2118,&
       366,  397,  443,  505,  576,  651,  755,  841,  973, 1095, 1210, 1308, 1435, 1545, 1702, 1786, 1878, 2081,&
       356,  387,  432,  489,  562,  635,  722,  841,  961, 1078, 1175, 1273, 1397, 1517, 1672, 1763, 1863, 2034,&
       359,  388,  431,  479,  545,  624,  707,  834,  943, 1059, 1158, 1256, 1368, 1490, 1647, 1747, 1849, 1998,&
       361,  394,  436,  483,  547,  615,  704,  817,  920, 1027, 1124, 1238, 1358, 1465, 1624, 1736, 1838, 1978,&
       400,  382,  440,  489,  546,  610,  690,  817,  897,  998, 1115, 1201, 1351, 1453, 1599, 1728, 1829, 1953,&
       474,  461,  416,  491,  549,  612,  701,  806,  883,  974, 1078, 1194, 1310, 1456, 1569, 1716, 1818, 1925,&
       531,  518,  507,  463,  547,  615,  694,  784,  881,  958, 1047, 1153, 1297, 1432, 1547, 1688, 1809, 1901,&
       594,  585,  577,  564,  513,  615,  695,  779,  879,  949, 1041, 1145, 1264, 1388, 1544, 1644, 1804, 1879,&
       675,  659,  651,  639,  632,  576,  695,  782,  879,  957, 1046, 1141, 1254, 1391, 1524, 1614, 1793, 1871,&
       739,  734,  726,  719,  715,  708,  663,  776,  901,  971, 1022, 1117, 1232, 1355, 1478, 1616, 1766, 1887,&
       819,  821,  805,  784,  773,  761,  736,  761,  888,  958, 1044, 1145, 1237, 1346, 1487, 1614, 1721, 1891,&
       899,  895,  871,  852,  856,  861,  854,  759,  883,  984, 1027, 1113, 1226, 1355, 1467, 1568, 1703, 1885,&
       973,  946,  955,  925,  939,  927,  902,  920,  870,  987, 1061, 1145, 1234, 1319, 1439, 1552, 1722, 1859/))
  !
  REAL(SP), PARAMETER :: AL_PD(18,18) = RESHAPE(SHAPE= (/18,18 /), SOURCE= &
       (/.281, .288, .283, .282, .278, .281, .272, .274, .268, .257, .251, .243, .246, .251, .254, .268, .304, .308,&  
       .290, .297, .291, .290, .286, .282, .277, .275, .267, .254, .252, .244, .250, .257, .260, .274, .308, .312,& 
       .294, .299, .293, .294, .288, .289, .281, .276, .265, .256, .251, .247, .258, .264, .268, .283, .318, .317,&  
       .297, .298, .302, .300, .289, .295, .290, .276, .264, .256, .260, .258, .268, .277, .281, .292, .330, .327,&  
       .305, .311, .313, .315, .305, .304, .299, .279, .271, .272, .273, .276, .285, .290, .293, .302, .340, .340,&  
       .292, .294, .303, .305, .301, .307, .290, .277, .274, .278, .287, .288, .295, .302, .306, .312, .343, .346,&  
       .268, .277, .279, .285, .285, .290, .279, .278, .280, .283, .295, .296, .305, .310, .313, .315, .342, .346,&  
       .288, .285, .280, .278, .278, .277, .272, .271, .279, .288, .297, .305, .310, .313, .311, .310, .335, .338,&  
       .314, .304, .292, .282, .275, .275, .262, .272, .290, .293, .299, .307, .308, .310, .303, .302, .325, .328,&  
       .346, .329, .313, .295, .283, .275, .264, .274, .288, .302, .307, .310, .306, .307, .292, .296, .315, .320,&  
       .320, .295, .326, .318, .294, .277, .275, .271, .293, .303, .305, .309, .309, .303, .294, .294, .310, .313,&  
       .304, .310, .297, .320, .317, .297, .283, .274, .298, .305, .308, .311, .313, .300, .290, .293, .305, .306,&  
       .314, .313, .308, .297, .325, .314, .293, .276, .292, .309, .314, .308, .303, .296, .286, .291, .301, .302,&  
       .308, .311, .307, .312, .288, .340, .305, .285, .294, .310, .315, .309, .296, .285, .281, .288, .298, .295,&  
       .313, .310, .315, .303, .313, .294, .331, .286, .294, .307, .320, .316, .303, .281, .278, .285, .290, .292,&  
       .315, .306, .308, .297, .295, .283, .334, .297, .280, .294, .314, .321, .313, .291, .280, .279, .287, .290,&  
       .308, .304, .305, .297, .279, .285, .251, .278, .278, .284, .297, .314, .307, .289, .274, .274, .274, .291,& 
       .301, .299, .298, .285, .265, .279, .241, .285, .260, .286, .302, .306, .302, .288, .277, .263, .271, .293/))
  
  !  
  ! d-f data from Barklem, O'Mara and Ross, 1998, MNRAS, 296, 1057
  !
  REAL(SP), PARAMETER :: CS_DF(18,18) = RESHAPE(SHAPE= (/18,18 /), SOURCE= &
       (/808,  873,  958, 1059, 1175, 1306, 1453, 1615, 1793, 1979, 2121, 2203, 2461, 2604, 2764, 2757, 2784, 3156,&
       798,  866,  953, 1052, 1172, 1299, 1450, 1606, 1776, 1967, 2114, 2196, 2451, 2601, 2763, 2767, 2783, 3142,&
       781,  848,  934, 1030, 1149, 1276, 1416, 1596, 1751, 1944, 2100, 2188, 2436, 2594, 2767, 2777, 2795, 3123,&
       766,  831,  915, 1010, 1124, 1239, 1398, 1564, 1729, 1912, 2083, 2180, 2426, 2585, 2776, 2790, 2808, 3106,&
       750,  814,  897,  987, 1097, 1201, 1355, 1530, 1718, 1875, 2060, 2171, 2414, 2575, 2779, 2809, 2820, 3103,&
       733,  797,  872,  950, 1049, 1166, 1326, 1502, 1670, 1851, 2026, 2165, 2396, 2562, 2779, 2827, 2832, 3099,&
       726,  786,  853,  936, 1011, 1128, 1303, 1472, 1649, 1844, 1979, 2159, 2371, 2548, 2778, 2840, 2848, 3103,&
       709,  783,  847,  912, 1002, 1093, 1270, 1419, 1606, 1787, 1951, 2139, 2335, 2533, 2775, 2847, 2863, 3104,&
       758,  721,  838,  907, 1010, 1066, 1211, 1401, 1600, 1774, 1972, 2098, 2313, 2528, 2781, 2857, 2892, 3121,&
       869,  882,  820,  870, 1003, 1098, 1165, 1368, 1527, 1735, 1896, 2030, 2288, 2534, 2776, 2844, 2902, 3123,&
       970,  967,  934,  938,  918, 1130, 1194, 1287, 1507, 1679, 1821, 2021, 2271, 2525, 2732, 2786, 2882, 3085,&
       1079, 1043, 1056, 1007, 1014, 1021, 1200, 1326, 1424, 1668, 1818, 1988, 2242, 2493, 2672, 2719, 2853, 3035,&
       1174, 1173, 1127, 1154, 1104, 1099, 1169, 1288, 1442, 1580, 1704, 1882, 2136, 2400, 2561, 2648, 2832, 2994,&
       1285, 1278, 1269, 1225, 1252, 1229, 1116, 1343, 1380, 1594, 1710, 1874, 2054, 2309, 2484, 2607, 2813, 2932,&
       1440, 1408, 1422, 1380, 1383, 1341, 1361, 1192, 1448, 1454, 1675, 1873, 2069, 2246, 2432, 2610, 2811, 2878,&
       1572, 1545, 1553, 1517, 1481, 1502, 1469, 1349, 1373, 1561, 1586, 1781, 2072, 2301, 2490, 2626, 2754, 2832,&
       1698, 1701, 1694, 1641, 1617, 1651, 1566, 1600, 1374, 1547, 1698, 1749, 1989, 2289, 2511, 2594, 2689, 2774,&
       1870, 1841, 1786, 1752, 1777, 1757, 1666, 1732, 1522, 1533, 1707, 1817, 1928, 2194, 2435, 2574, 2665, 2742/))
  !
  REAL(SP), PARAMETER :: AL_DF(18,18) = RESHAPE(SHAPE= (/18,18 /), SOURCE= &
       (/.295, .286, .299, .300, .307, .310, .311, .311, .316, .319, .325, .351, .364, .369, .372, .379, .373, .351,&
       .295, .295, .301, .302, .311, .316, .314, .314, .320, .321, .324, .349, .361, .365, .368, .374, .368, .349,&  
       .286, .298, .302, .304, .311, .323, .321, .319, .324, .323, .323, .345, .355, .358, .362, .367, .361, .343,& 
       .290, .295, .307, .316, .322, .329, .326, .325, .329, .324, .321, .343, .350, .351, .354, .360, .358, .337,&  
       .292, .299, .307, .321, .327, .336, .333, .330, .330, .320, .321, .338, .344, .344, .345, .352, .352, .332,&  
       .291, .299, .309, .323, .335, .339, .335, .333, .327, .323, .319, .333, .336, .336, .336, .344, .345, .329,& 
       .297, .302, .312, .321, .340, .338, .333, .327, .325, .319, .318, .324, .329, .330, .330, .336, .337, .325,& 
       .319, .314, .317, .327, .334, .344, .339, .327, .323, .318, .312, .318, .319, .322, .322, .326, .327, .316,&  
       .333, .328, .339, .325, .359, .351, .332, .325, .322, .311, .309, .310, .311, .316, .314, .317, .321, .313,& 
       .274, .273, .323, .412, .318, .339, .359, .328, .324, .311, .309, .325, .322, .315, .318, .319, .325, .314,&  
       .297, .296, .273, .302, .436, .325, .354, .335, .326, .311, .314, .330, .323, .324, .325, .323, .330, .314,& 
       .284, .295, .296, .280, .300, .438, .322, .348, .332, .318, .320, .332, .335, .334, .335, .331, .333, .309,& 
       .280, .278, .285, .297, .279, .320, .445, .319, .320, .324, .328, .338, .348, .346, .345, .336, .328, .300,& 
       .280, .273, .267, .273, .284, .268, .343, .390, .323, .308, .318, .325, .343, .348, .346, .337, .311, .286,&  
       .277, .270, .260, .266, .276, .263, .294, .408, .337, .324, .299, .308, .331, .334, .345, .327, .315, .280,& 
       .270, .262, .258, .260, .273, .273, .262, .375, .410, .298, .312, .294, .313, .331, .328, .322, .307, .270,&  
       .271, .267, .262, .264, .274, .269, .261, .323, .351, .359, .294, .325, .310, .318, .321, .315, .291, .268,& 
       .275, .276, .272, .276, .279, .270, .264, .295, .393, .340, .319, .287, .320, .330, .316, .302, .280, .261 /))
  !
  PUBLIC :: GET_DAMPING
  PUBLIC :: GET_NEFF
  PUBLIC :: GET_COLLISIONAL_PARAM
  !
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! get_damping
  ! get_neff
  ! get_collisional_param
  ! transition_sp
  ! transition_pd
  ! transition_df
  ! subroutine
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_DAMPING(I,LAMBDA,ION,EPLOW,NHYD,TEMP,SIGMA,ALPHA,DAMP,DDAMDT_PG,DDAMDT_RHO,DDAMDP_TEM)
    !
    IMPLICIT NONE
    !
    INTEGER,    INTENT(IN)              :: I, ION
    REAL(SP),   INTENT(IN)              :: TEMP
    REAL(DP),   INTENT(IN)              :: EPLOW, NHYD, SIGMA, ALPHA, LAMBDA
    REAL(DP),   INTENT(INOUT)             :: DAMP
    REAL(DP),   INTENT(INOUT), OPTIONAL   :: DDAMDT_RHO, DDAMDT_PG, DDAMDP_TEM
    ! Internal
    REAL(SP),   PARAMETER               :: RAT_POL_HE_H = 0.3181818
    !IS IT USED AT ALL?   ;REAL(SP),   PARAMETER               :: RAT_POL_H2_H = 1.2121212
    REAL(DP),   PARAMETER               :: V0 = 1D6 ! cm/s
    REAL(DP)                            :: MU_X_H, MU_X_HE, E1, E2, GAMMAF, BETA
    REAL(DP)                            :: CHYD, CHYDP, EPUPP, RINDEX, LAMBDAVAC
    REAL(DP)                            :: GAMMA_COL, GAMMA_RAD, VDOPPLER
    REAL(DP)                            :: NHE, EDIFF1, EDIFF2
    REAL(DP)                            :: DVDOPDTEMP, DGAMMACDTEMP_PG, DGAMMACDTEMP_RHO
    LOGICAL                             :: NAN, CHECKNAN
    !
    REAL(DP)                            :: T8
    !
    T8=DBLE(TEMP)
    !
    ! Determine Doppler width (in velocity units: cm/s)
    !
    ! THERE SEEMS TO BE A MISTAKE HERE:
    VDOPPLER = SQRT(2.D0*KBOL*T8/(MAMU*MATOM(I)))
    ! THE SQRT(SPI) FACTOR DOES MATCH NEITHER WITTMANN'S NOR LANDIS'S FORMALISM
    !
    !VDOPPLER = SQRT(2.*REAL(KBOL)*TEMP/(REAL(MAMU)*MATOM(I)*SPI))
    !
    !
    ! Determine radiative damping
    CALL RADIATIVE_DAMPING(LAMBDA,GAMMA_RAD)
    ! Reduced mass of atom under consideration and Hydrogen
    MU_X_H = MAMU*(MATOM(I)*MATOM(1))/(MATOM(I)+MATOM(1))
    ! Reduced mass of atom under consideration and Helium
    MU_X_HE = MAMU*(MATOM(I)*MATOM(2))/(MATOM(I)+MATOM(2))
    ! Calculate number helium atoms per cm^3
    NHE = NHYD*10D0**(ABUND(2)-12D0)
    !
    ! Calculate Collisional damping according to ABO Theory: Anstee, Barklem & O'Mara
    !
    IF (SIGMA.GT.0.AND.ALPHA.GT.0) THEN
       !
       E1 = 2.D0-ALPHA/2.D0
       E2 = (1.D0-ALPHA)/2.D0
       ! Gamma function: 5th order polynomial interpolation over the [1,2] region
       !GAMMAF=3.6924374-6.5872093*E1+6.3611214*E1**2.-3.2531344*E1**3.+0.88197419*E1**4.-0.095283853*E1**5.
       GAMMAF=1.D0+(-0.5748646D0+(0.9512363D0+(-0.6998588D0&
           +(0.4245549D0-0.1010678D0*E1)*E1)*E1)*E1)*E1
       ! Beta factor (if left like this GAMMA_COL will be in HWHM)
       BETA=(4.D0/DPI)**(ALPHA/2.D0)*GAMMAF*(V0**ALPHA)*(SIGMA*RBOHR**2.D0)&
           *(8.D0*KBOL/DPI)**E2
       ! Now it will be in FWHM
       BETA=2.D0*BETA
       ! Calculate collisional damping: this is already FWHM
       GAMMA_COL=BETA*T8**E2*(NHYD*MU_X_H**(-E2)&
           +RAT_POL_HE_H*NHE*MU_X_HE**(-E2))
       ! Units of GAMMA_COL are s^(-1)
       !
       ! Now we deal with the derivatives of the collisional damping...
       ! ... with respect to T and constant PG and RHO
       IF (MRESPFUNCT.EQV..TRUE.) THEN
         DGAMMACDTEMP_PG=GAMMA_COL*E2/T8+BETA*T8**E2&
             *DNHYDDTEMP_PG*(MU_X_H**(-E2)+RAT_POL_HE_H&
             *MU_X_HE**(-E2)*10D0**(ABUND(2)-12D0))
         DGAMMACDTEMP_RHO=GAMMA_COL*E2/T8+BETA*T8**E2&
             *DNHYDDTEMP_RHO*(MU_X_H**(-E2)+RAT_POL_HE_H&
             *MU_X_HE**(-E2)*10D0**(ABUND(2)-12D0))
       ENDIF
       !
    ELSE
       ! Change from air to vacuum wavelength
       RINDEX=1.0004
       IF (LAMBDA*1D8.GE.1800.D0) CALL REFRAX(LAMBDA*1D4,RINDEX)
       LAMBDAVAC=DBLE(RINDEX*LAMBDA)
       ! Excitation potential of the upper level in vacuum
       EPUPP = EPLOW+(HPLA*LIGHT/LAMBDAVAC)/EVOLT
       IF (ION.EQ.1) THEN
          EDIFF1=MAX(XI(I)-EPUPP-XI(I)*REAL(ION-1),1.0D0)
          EDIFF2=MAX(XI(I)-EPLOW-XI(I)*REAL(ION-1),3.0D0)
          CHYD=1.0
       ELSE IF (ION.EQ.2) THEN
          EDIFF1=MAX(XII(I)-EPUPP-XI(I)*REAL(ION-1),1.0D0)
          EDIFF2=MAX(XII(I)-EPLOW-XI(I)*REAL(ION-1),3.0D0)
          CHYD=1.741
       ELSE
          EDIFF1=0.0D0
          EDIFF2=0.0D0
          WRITE(*,*) 'ONLY NEUTRAL AND ONCE IONIZED STATES ALLOWED, STOPPING!'
          STOP
       ENDIF 
       !
       CHYD=CHYD*17.0D0*(8.0D0*KBOL/DPI)**0.3D0
       CHYD=CHYD*(10.0D0**(-30.5325D0)*(1.0D0/EDIFF1**2.0D0&
           -1.0D0/EDIFF2**2.0D0))**0.4D0
       CHYDP=CHYD*(NHYD*MU_X_H**(-0.3D0)+NHE*MU_X_HE**(-0.3D0))
       GAMMA_COL=CHYDP*TEMP**0.3D0
       !
       ! Derivatives:
       IF (MRESPFUNCT.EQV..TRUE.) THEN
         !
         DGAMMACDTEMP_RHO=CHYD*(MU_X_H**(-0.3D0)+10D0**(ABUND(2)-12D0)&
              *MU_X_HE**(-0.3D0))*(DNHYDDTEMP_RHO*T8**0.3D0&
              +0.3D0*NHYD*T8**(-0.7D0))
         DGAMMACDTEMP_PG=CHYD*(MU_X_H**(-0.3D0)+10D0**(ABUND(2)-12D0)&
              *MU_X_HE**(-0.3D0))*(DNHYDDTEMP_PG*T8**0.3D0&
              +0.3D0*NHYD*T8**(-0.7D0))
         !
       ENDIF
    END IF
    ! Finally add both to obtain total damping
    DAMP = (GAMMA_RAD+GAMMA_COL)*(LAMBDA)/(4.D0*DPI*VDOPPLER)
    ! Check NAN internally?
    !NAN=.FALSE.
    !NAN=CHECKNAN(DBLE(DAMP))
    !IF (NAN.EQV..TRUE.) THEN
    IF (DAMP.NE.DAMP) THEN
       NAN=.FALSE.
       NAN=CHECKNAN(DBLE(DAMP))
       PRINT*,'NaN detected in DAMP'
       PRINT*,DAMP, TEMP, GAMMA_RAD, GAMMA_COL, LAMBDA, VDOPPLER
       PRINT*, 'Here?'
WRITE(*,*) GAMMA_COL, BETA, T8, E2, NHYD, MU_X_H, RAT_POL_HE_H, NHE,MU_X_HE
       STOP
    ENDIF
    !


    !
    ! Derivatives:
    IF (MRESPFUNCT.EQV..TRUE.) THEN
      !
      DVDOPDTEMP = VDOPPLER/(2D0*T8)
      ! Now derivatives of teh total damping with respect to T...
      ! ...and constant PG and RHO
      DDAMDTEMP_PG=DAMP/(GAMMA_RAD+GAMMA_COL)*DGAMMACDTEMP_PG&
          -DAMP/VDOPPLER*DVDOPDTEMP
      DDAMDTEMP_RHO=DAMP/(GAMMA_RAD+GAMMA_COL)*DGAMMACDTEMP_RHO&
          -DAMP/VDOPPLER*DVDOPDTEMP
      DDAMDPG_TEMP=LAMBDA/(4.D0*DPI*VDOPPLER)*GAMMA_COL/NHYD*DNHYDDPG_TEMP
! .NEW
      DDAMDRHO_TEMP=LAMBDA/(4.D0*DPI*VDOPPLER)*GAMMA_COL/NHYD*DNHYDDRHO_TEMP
! NEW.
      !
      IF (PRESENT(DDAMDT_PG).AND.PRESENT(DDAMDT_RHO).AND.PRESENT(DDAMDP_TEM)) THEN
         DDAMDT_PG  = DDAMDTEMP_PG
         DDAMDT_RHO = DDAMDTEMP_RHO
         DDAMDP_TEM = DDAMDPG_TEMP
      ENDIF
    ENDIF
    !
    !
    !
    ! Calculate Collisional damping according to Unsold theory: following Witmann 1972
    !
    !IF (SIGMA.EQ.0.OR.ALPHA.EQ.0) THEN
    !   IF (ION.EQ.1) THEN 
    !      E1=MAXVAL(XI(I)-EPLOW-XI(I)*REAL(ION-1.),1.)
    !      E2=MAXVAL(XI(I)-EPLOW-XI(I)*REAL(ION-1.),3.)
    !   ENDIF
    !   IF (ION.EQ.2) THEN
    !      E1=MAXVAL(XII(I)-EPLOW-XI(I)*REAL(ION-1.),1.)
    !      E2=MAXVAL(XII(I)-EPLOW-XI(I)*REAL(ION-1.),3.)
    !   ENDIF
    !   CHYD = LAMBDA*10.**(0.4*ALOG10(1./E1**2.-1./E2**2.)-12.213)*5.34784E3
    !   IF (ION.EQ.2) CHYD=CHYD*1.741
    !ENDIF
  END SUBROUTINE GET_DAMPING
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_NEFF(ZN,EPLOW,L0,ION,NLOW_EFF,NUPP_EFF,ELOWLIM,EUPPLIM)
    !
USE user_mpi
    IMPLICIT NONE
    INTEGER,  INTENT(IN)               :: ZN, ION
    REAL(DP), INTENT(IN)               :: EPLOW, L0
    REAL(DP), INTENT(IN), OPTIONAL     :: ELOWLIM,EUPPLIM
    REAL(DP), INTENT(OUT)              :: NLOW_EFF, NUPP_EFF
    REAL(DP)                           :: XH, XLIMIT, XLIMITLOW, XLIMITUPP, XLOW, XUPP
    !
    XLOW=EPLOW*EVOLT                  ! Energy (erg) of the lower level
    XUPP=XLOW+HPLA*LIGHT/(L0*1E-8)    ! Energy (erg) of the upper level
    XH=XII(1)*EVOLT                   ! Ionization potential of hydrogen
    ! -------------------------------------------------------------------------
    ! Series Limit: not exactly ionization potential, but in general very close
    ! -------------------------------------------------------------------------
    IF (PRESENT(ELOWLIM) .AND. PRESENT(EUPPLIM)) THEN
       XLIMITLOW=ELOWLIM
       XLIMITUPP=EUPPLIM
    ELSE
       ! If Z > 1
       IF (ZN.GT.1) THEN
          ! If single ionized:
          IF (ION.EQ.1) THEN
             XLIMIT=DBLE(XI(ZN))*EVOLT
          ELSE IF (ION.EQ.2) THEN
             XLIMIT=DBLE(XII(ZN))*EVOLT
          ELSE
             PRINT*, 'MORE THAN ONCE IONIZED ELEMENTS ARE NOT IMPLEMENTED YET!'
             STOP
          ENDIF
       ELSE IF (ZN.EQ.1) THEN
          IF (ION.GT.1) THEN
             PRINT*,'Hydrogen cannot be ionized twice !'
             PRINT*,'Error in file lines_database.dat. STOP'
             STOP
          ELSE
             XLIMIT=DBLE(XII(ZN))*EVOLT
          ENDIF
       ELSE
          PRINT*, 'ION must be greater than 0!'
          STOP
       ENDIF
       !! Single ionized elements Z > 1
       !IF (ION.EQ.1.AND.ZN.GT.1) XLIMIT=DBLE(XI(ZN))*EVOLT
       !! Double ionized elements Z > 1
       !IF (ION.EQ.2.AND.ZN.GT.1) XLIMIT=DBLE(XII(ZN))*EVOLT
       !! Case of hydrogen Z = 1
       !!IF (ZN.EQ.1) XLIMIT=DBLE(XII(ZN))*EVOLT
       !!IF (ZN.EQ.1.AND.ION.GT.1) THEN
       !!   PRINT*,'Hydrogen cannot be ionized twice !'
       !!   PRINT*,'Error in file lines_database.dat. STOP'
       !!   STOP
       !!ENDIF
       ! Converting from ergs to cm-1
       XLIMITLOW=XLIMIT/(HPLA*LIGHT)
       XLIMITUPP=XLIMIT/(HPLA*LIGHT)
    ENDIF
    ! Converting from ergs to cm-1
    XLOW=XLOW/(HPLA*LIGHT)
    XUPP=XUPP/(HPLA*LIGHT)
    XH=XH/(HPLA*LIGHT)
    !
    ! Effective quantum numbers
    !
    NLOW_EFF=SQRT(XH/(XLIMITLOW-XLOW))
    NUPP_EFF=SQRT(XH/(XLIMITUPP-XUPP))
IF (mpi__myrank.EQ.0) THEN
  IF ( (NLOW_EFF.NE.NLOW_EFF) .OR. (NUPP_EFF.NE.NUPP_EFF) ) THEN
    PRINT*, NLOW_EFF, XH, XLIMITLOW, XLOW
    PRINT*, NUPP_EFF, XH, XLIMITUPP, XUPP
  ENDIF
ENDIF
    !
  END SUBROUTINE GET_NEFF
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_COLLISIONAL_PARAM(NLOW_EFF,NUPP_EFF,OTRANSITION,SIGMA,ALPHA,IERR)
    IMPLICIT NONE
    !
    REAL(DP),     INTENT(IN)    :: NLOW_EFF, NUPP_EFF
    INTEGER,      INTENT(IN)    :: OTRANSITION(2)
    REAL(DP),     INTENT(OUT)   :: SIGMA, ALPHA
    INTEGER,      INTENT(INOUT) :: IERR
    INTEGER                     :: TRANSITION_TYPE
    !
    ! IERR values = 1 (not a sp, pd, df transition)
    !             = 2 (lower level outside table ranges)
    !             = 3 (upper level outside table ranges)
    TRANSITION_TYPE = SUM(OTRANSITION)
    SELECT CASE(TRANSITION_TYPE)
    CASE(1)
       ! s-p
       CALL TRANSITION_SP(NLOW_EFF,NUPP_EFF,OTRANSITION,SIGMA,ALPHA,IERR)
    CASE(3)
       ! p-d
       CALL TRANSITION_PD(NLOW_EFF,NUPP_EFF,OTRANSITION,SIGMA,ALPHA,IERR)
    CASE(5)
       ! d-f
       CALL TRANSITION_DF(NLOW_EFF,NUPP_EFF,OTRANSITION,SIGMA,ALPHA,IERR)
    CASE DEFAULT
       IERR = 1
    END SELECT
    !print*,nlow_eff,nupp_eff,otransition,sigma,alpha,ierr
    !
  END SUBROUTINE GET_COLLISIONAL_PARAM
  !
  !------------------------------------------------
  !
  SUBROUTINE TRANSITION_SP(NLOW_EFF,NUPP_EFF,OTRANSITION,SIGMA,ALPHA,IERR)
    IMPLICIT NONE
    !
    REAL(DP),     INTENT(IN)    :: NLOW_EFF, NUPP_EFF
    REAL(DP),     INTENT(OUT)   :: SIGMA, ALPHA
    INTEGER,      INTENT(INOUT) :: IERR
    INTEGER,      INTENT(IN)    :: OTRANSITION(2)
    !
    REAL(DP)     :: CS(21,18), AL(21,18), CS2(21,18), AL2(21,18)
    REAL(DP)     :: NSSG(21), NSPG(18), NSS, NSP
    REAL(DP)     :: A1, A2, A3, A4, S1, S2, S3, S4, AT
    INTEGER      :: I, J, IX, IY
    !
    IF (OTRANSITION(1) .EQ. 0) NSS=NLOW_EFF
    IF (OTRANSITION(1) .EQ. 1) NSP=NLOW_EFF
    IF (OTRANSITION(2) .EQ. 0) NSS=NUPP_EFF
    IF (OTRANSITION(2) .EQ. 1) NSP=NUPP_EFF
    ! Check table limits
    IF ((NSS .GT. 3.) .OR. (NSS .LT. 1.) .OR. (NSP .GT. 3.) .OR. (NSP .LT. 1.3)) THEN
       IERR=3
    ENDIF
    !
    IF (IERR .EQ. 0) THEN
       !
       DO I=1,21
          NSSG(I)=1.0+REAL(I-1)*0.1
       ENDDO
       !
       DO I=1,18
          NSPG(I)=1.3+REAL(I-1)*0.1
       ENDDO
       !
       DO I=1,21
          DO J=1,18
             ! SIGMA-table
             CS(I,J)=CS_SP(J,I)
             ! ALPHA-table
             AL(I,J)=AL_SP(J,I)
          ENDDO
       ENDDO
       !
       IX=MINLOC(ABS(NSS-NSSG),dim=1)
       IY=MINLOC(ABS(NSP-NSPG),dim=1)
       !
       IF (NSS < NSSG(IX)) IX=IX-1
       IF (NSP < NSSG(IY)) IY=IY-1
       ! Bilinear interpolation: SIGMA
       S1=CS(IX,IY)
       S2=CS(IX+1,IY)
       S3=CS(IX+1,IY+1)
       S4=CS(IX,IY+1)
       A1=ABS((NSSG(IX+1)-NSS)*(NSPG(IY+1)-NSP))
       A2=ABS((NSSG(IX)-NSS)*(NSPG(IY+1)-NSP))
       A3=ABS((NSSG(IX)-NSS)*(NSPG(IY)-NSP))
       A4=ABS((NSSG(IX+1)-NSS)*(NSPG(IY)-NSP))
       AT=A1+A2+A3+A4
       SIGMA=(S1*A1+S2*A2+S3*A3+S4*A4)/AT
       ! Bilinear interpolation: ALPHA
       S1=AL(IX,IY)
       S2=AL(IX+1,IY)
       S3=AL(IX+1,IY+1)
       S4=AL(IX,IY+1)
       ALPHA=(S1*A1+S2*A2+S3*A3+S4*A4)/AT
       !PRINT*,'Bilinear SIGMA=',SIGMA
       !PRINT*,'Bilinear ALPHA=',ALPHA
    ENDIF
    !
  END SUBROUTINE TRANSITION_SP
  !
  !------------------------------------------------
  !
  SUBROUTINE TRANSITION_PD(NLOW_EFF,NUPP_EFF,OTRANSITION,SIGMA,ALPHA,IERR)
    IMPLICIT NONE
    REAL(DP),     INTENT(IN)    :: NLOW_EFF, NUPP_EFF
    REAL(DP),     INTENT(OUT)   :: SIGMA ,ALPHA
    INTEGER,      INTENT(IN)    :: OTRANSITION(2)
    INTEGER,      INTENT(INOUT) :: IERR
    !
    REAL(DP)     :: CS(18,18), AL(18,18), CS2(18,18), AL2(18,18)
    REAL(DP)     :: NSDG(18), NSPG(18), NSD, NSP
    REAL(DP)     :: A1, A2, A3, A4, S1, S2, S3, S4, AT
    INTEGER      :: I, J, IX, IY
    !
    IF (OTRANSITION(1) .EQ. 1) NSP=NLOW_EFF
    IF (OTRANSITION(1) .EQ. 2) NSD=NLOW_EFF
    IF (OTRANSITION(2) .EQ. 1) NSP=NUPP_EFF
    IF (OTRANSITION(2) .EQ. 2) NSD=NUPP_EFF
    ! Check table limits
    IF ((NSP .GT. 3.) .OR. (NSP .LT. 1.3) .OR. (NSD .GT. 4.) .OR. (NSD .LT. 2.3)) THEN
       IERR=3
    ENDIF
    !
    IF (IERR .EQ. 0) THEN
       !
       DO I=1,18
          NSPG(I)=1.3+REAL(I-1)*0.1
          NSDG(I)=2.3+REAL(I-1)*0.1
          DO J=1,18
             CS(I,J)=CS_PD(J,I)
             AL(I,J)=AL_PD(J,I)
          ENDDO
       ENDDO
       IX=MINLOC(ABS(NSP-NSPG), dim=1)
       IY=MINLOC(ABS(NSD-NSDG), dim=1)
       !
       IF (NSP < NSPG(IX)) IX=IX-1
       IF (NSD < NSDG(IY)) IY=IY-1
       ! Bilinear interpolation: SIGMA
       S1=CS(IX,IY)
       S2=CS(IX+1,IY)
       S3=CS(IX+1,IY+1)
       S4=CS(IX,IY+1)
       A1=ABS((NSPG(IX+1)-NSP)*(NSDG(IY+1)-NSD))
       A2=ABS((NSPG(IX)-NSP)*(NSDG(IY+1)-NSD))
       A3=ABS((NSPG(IX)-NSP)*(NSDG(IY)-NSD))
       A4=ABS((NSPG(IX+1)-NSP)*(NSDG(IY)-NSD))
       AT=A1+A2+A3+A4
       SIGMA=(S1*A1+S2*A2+S3*A3+S4*A4)/AT
       ! Bilinear interpolation: ALPHA
       S1=AL(IX,IY)
       S2=AL(IX+1,IY)
       S3=AL(IX+1,IY+1)
       S4=AL(IX,IY+1)
       ALPHA=(S1*A1+S2*A2+S3*A3+S4*A4)/AT
       !PRINT*,'Bilinear SIGMA=',SIGMA
       !PRINT*,'Bilinear ALPHA=',ALPHA   
    ENDIF
    !
  END SUBROUTINE TRANSITION_PD
  !
  !------------------------------------------------
  !
  SUBROUTINE TRANSITION_DF(NLOW_EFF,NUPP_EFF,OTRANSITION,SIGMA,ALPHA,IERR)
    IMPLICIT NONE
    REAL(DP),     INTENT(IN)    :: NLOW_EFF, NUPP_EFF
    REAL(DP),     INTENT(OUT)   :: SIGMA, ALPHA
    INTEGER,      INTENT(IN)    :: OTRANSITION(2)
    INTEGER,      INTENT(INOUT) :: IERR
    !
    REAL(DP)     :: CS(18,18), AL(18,18), CS2(18,18), AL2(18,18)
    REAL(DP)     :: NSDG(18), NSFG(18), NSD, NSF
    REAL(DP)     :: A1, A2, A3, A4, S1, S2, S3, S4, AT
    INTEGER      :: I, J, IX, IY
    !
    IF (OTRANSITION(1) .EQ. 2) NSD=NLOW_EFF
    IF (OTRANSITION(1) .EQ. 3) NSF=NLOW_EFF
    IF (OTRANSITION(2) .EQ. 2) NSD=NUPP_EFF
    IF (OTRANSITION(2) .EQ. 3) NSF=NUPP_EFF
    ! Check table limits
    IF ((NSD .GT. 4.) .OR. (NSD .LT. 2.3) .OR. (NSF .GT. 5.) .OR. (NSF .LT. 3.3)) THEN
       PRINT*, NSD, 4
       PRINT*, NSD, 2.3
       PRINT*, NSF, 5
       PRINT*, NSF, 3.3
       IERR=3
    ENDIF
    !
    IF (IERR .EQ. 0) THEN
       DO I=1,18
          NSDG(I)=2.3+REAL(I-1)*0.1
          NSFG(I)=3.3+REAL(I-1)*0.1
          DO J=1,18
             CS(I,J)=CS_DF(J,I)
             AL(I,J)=AL_DF(J,I)
          ENDDO
       ENDDO
       IX=MINLOC(ABS(NSD-NSDG), dim=1)
       IY=MINLOC(ABS(NSF-NSFG), dim=1)
       !
       IF (NSD < NSDG(IX)) IX=IX-1
       IF (NSF < NSFG(IY)) IY=IY-1
       ! Bilinear interpolation: SIGMA
       S1=CS(IX,IY)
       S2=CS(IX+1,IY)
       S3=CS(IX+1,IY+1)
       S4=CS(IX,IY+1)
       A1=ABS((NSDG(IX+1)-NSD)*(NSFG(IY+1)-NSF))
       A2=ABS((NSDG(IX)-NSD)*(NSFG(IY+1)-NSF))
       A3=ABS((NSDG(IX)-NSD)*(NSFG(IY)-NSF))
       A4=ABS((NSDG(IX+1)-NSD)*(NSFG(IY)-NSF))
       AT=A1+A2+A3+A4
       SIGMA=(S1*A1+S2*A2+S3*A3+S4*A4)/AT
       ! Bilinear interpolation: ALPHA
       S1=AL(IX,IY)
       S2=AL(IX+1,IY)
       S3=AL(IX+1,IY+1)
       S4=AL(IX,IY+1)
       ALPHA=(S1*A1+S2*A2+S3*A3+S4*A4)/AT
       !PRINT*,'Bilinear SIGMA=',SIGMA
       !PRINT*,'Bilinear ALPHA=',ALPHA
    ENDIF
    !
  END SUBROUTINE TRANSITION_DF
  !
  !------------------------------------------------
  !
  PURE SUBROUTINE RADIATIVE_DAMPING(LAMBDA,GAMMA_RAD)
    !
    IMPLICIT NONE
    !
    REAL(DP),   INTENT(IN)      :: LAMBDA
    REAL(DP),   INTENT(INOUT)     :: GAMMA_RAD
    ! Calculate radiative damping assuming a classical oscillator
    ! LAMBDA must be in cm for GAMMA_RAD to be in s^(-1)
    GAMMA_RAD=0.2223D0/LAMBDA**2.D0
    !
  END SUBROUTINE RADIATIVE_DAMPING
  !
  !================================================
  !
END MODULE DAMPING
!
