#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 00:56:33 2024

@author: Alexandra
"""


import numpy as np
from scipy.stats import chisquare
from scipy.stats import chi2_contingency 
import pandas as pd
import copy
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Helvetica'

clade_6 = ['P02LT-2109', 'P02LT-2309', 'P02LT-2209', 'P02LT-2205', 'P02LT-3114', 'P02LT-2412', 'P02LT-3513', 'P02LT-3510', 'P02LT-2416', 'P02LT-0512', 'P02LT-3106', 'P02LT-2407', 'P02LT-0713', 'P02SL-0308', 'P02BL-B308', 'P02LT-1905', 'P02LT-3601', 'P02LG-2012', 'P02LT10510', 'P02LT-1501', 'P02LT-1808', 'P02LT-1801', 'P02LT10414', 'P02LT10310', 'P02LT-2809', 'P02LT-2715', 'P02LT-2709', 'P02LT-2507', 'P02LT-2414', 'P02LT-2404', 'P02LT-2308', 'P02LT-0803', 'P02LT-3506', 'P02LT-3703', 'P02LT-0612', 'P02LT-0611', 'P02LT-0609', 'P02LT-0605', 'P02LT-0603', 'P02BL-B724', 'P02BL-B723', 'P02BL-B722', 'P02BL-B721', 'P02BL-B720', 'P02BL-B719', 'P02BL-B718', 'P02BL-B717', 'P02BL-B716', 'P02BL-B715', 'P02BL-B714', 'P02BL-B713', 'P02BL-B712', 'P02BL-B711', 'P02BL-B710', 'P02BL-B709', 'P02BL-B708', 'P02BL-B707', 'P02BL-B704', 'P02BL-B703', 'P02BL-B702', 'P02BL-B701', 'P02LG-0108', 'P02LT-1912', 'P02LT10401', 'P02LT10308', 'P02LT-2810', 'P02LT-2806', 'P02LT-2804', 'P02LT-2615', 'P02LT-2608', 'P02LT-2605', 'P02LT-2601', 'P02LT-3010', 'P02LT-1807', 'P02LT-4103', 'P02LT-2610', 'P02LT-2609', 'P02SL-0423', 'P02LT-3805', 'P02LT10509', 'P02LT-4206', 'P02SL-0401', 'P02LT10307', 'P02LT-3905', 'P02LT-3904', 'P02LT-3516', 'P02LT-3514', 'P02LT-3511', 'P02LT-3310', 'P02LT-3308', 'P02LT-3304', 'P02LT-0709', 'P02LT-4114', 'P02LT10301', 'P02LT-3303', 'P02LT-1907', 'P02LT-3107', 'P02LT-3808', 'P02LG-0101', 'P02LT10514', 'P02LT10507', 'P02LT-1507', 'P02LT10412', 'P02LT10405', 'P02LT-4106', 'P02LT-2815', 'P02LT-2814', 'P02LT-2813', 'P02LT-2811', 'P02LT-2808', 'P02LT-2807', 'P02LT-2802', 'P02LT-2801', 'P02LT-2506', 'P02LT-2413', 'P02LT-2410', 'P02LT-2401', 'P02LT-3111', 'P02LT-3110', 'P02LT-3810', 'P02LT-3711', 'P02LT-3709', 'P02LT-3707', 'P02LT-3705', 'P02LT-3701', 'P02LT-0601']
clade_2 = ['P02LT-2713', 'P02LT-2710', 'P02LT-2708', 'P02BL-B812', 'P02BL-B806', 'P02BL-B803', 'P02LT-3012', 'P02LT-3316', 'P02LT-3301', 'P02LT-2913', 'P02LT-2910', 'P02LT-2905', 'P02LT-2901', 'P02LT-1406', 'P02LG-0117', 'P02SL-0420', 'P02LT10503', 'P02LT10501', 'P02LT10309', 'P02LT-0508', 'P02LT-1115', 'P02LT-1114', 'P02LT-1113', 'P02LT-1112', 'P02LT-1110', 'P02LT-1108', 'P02LT-1107', 'P02LT-1106', 'P02LT-1105', 'P02LT-1104', 'P02LT-1103', 'P02LT-1102', 'P02LT-1101', 'P02LT-0805', 'P02LT-0714', 'P02LT-0711', 'P02BLB4e07', 'P02BLB4e05', 'P02BLB4e04', 'P02BLB4e03', 'P02BLB4e02', 'P02LG-2021', 'P02BL-B415', 'P02BL-B414', 'P02BL-B411', 'P02BL-B410', 'P02BL-B409', 'P02BL-B407', 'P02BL-B406', 'P02BL-B404', 'P02BL-B403', 'P02BL-B402', 'P02BL-B401', 'P02LT10511', 'P02LT-1513', 'P02LT-1512', 'P02LT-1502', 'P02LT-1901', 'P02LT10306', 'P02LT-2405', 'P02LT-0516', 'P02LT-1111', 'P02LT-0915', 'P02LT-0911', 'P02LT-0908', 'P02LT-0902', 'P02LT-0802', 'P02LT-0801', 'P02LT-3108', 'P02LT-3616', 'P02LT-3612', 'P02LT-3607', 'P02LT-3602', 'P02BLB2e06', 'P02BLB2e05', 'P02BLB2e04', 'P02BLB2e03', 'P02BLB2e02', 'P02BLB2e01', 'P02SL-0302', 'P02BL-B222', 'P02BL-B221', 'P02BL-B220', 'P02BL-B214', 'P02BL-B213', 'P02BL-B210', 'P02BL-B209', 'P02BL-B208', 'P02BL-B206', 'P02BL-B204', 'P02BL-B202', 'P02BL-B201', 'P02LT10502', 'P02LT-1909', 'P02LT-1903', 'P02LT-3112', 'P02LT-3109', 'P02LT-3714', 'P02LT-3708', 'P02LT-0715', 'P02LT-0712', 'P02LT-0710', 'P02BL-B303', 'P02LT10411', 'P02LT-2303', 'P02LT-0705', 'P02LT-1404', 'P02LT-1402', 'P02LT-1401', 'P02LT-1216', 'P02LT-1215', 'P02LT-1214', 'P02LT-1213', 'P02LT-1212', 'P02LT-1211', 'P02LT-1210', 'P02LT-1209', 'P02LT-1208', 'P02LT-1207', 'P02LT-1206', 'P02LT-1205', 'P02LT-1204', 'P02LT-1203', 'P02LT-1202', 'P02LT-1201', 'P02LT-3615', 'P02LT-3608', 'P02BLB3e15', 'P02BL-B816', 'P02BL-B815', 'P02BL-B814', 'P02BL-B805', 'P02BL-B802', 'P02BL-B801', 'P02BL-B624', 'P02BL-B622', 'P02BL-B621', 'P02BL-B620', 'P02BL-B619', 'P02BL-B616', 'P02BL-B614', 'P02BL-B612', 'P02BL-B611', 'P02BL-B610', 'P02BL-B609', 'P02BL-B608', 'P02BL-B607', 'P02BL-B606', 'P02BL-B604', 'P02BL-B603', 'P02BL-B602', 'P02BL-B601', 'P02SL-0323', 'P02SL-0316', 'P02SL-0312', 'P02LG-0120', 'P02LG-0112', 'P02LG-0110', 'P02LG-0107', 'P02LG-0102', 'P02LG-2017', 'P02LG-2001', 'P02SL-0419', 'P02SL-0415', 'P02SL-0413', 'P02SL-0411', 'P02SL-0406', 'P02SL-0404', 'P02LT-3016', 'P02LT-3014', 'P02LT-3009', 'P02LT-3005', 'P02LT-3004', 'P02LT-3003', 'P02LT-3002', 'P02LT-1716', 'P02LT-1714', 'P02LT-1708', 'P02LT-1707', 'P02BL-B301', 'P02BL-B223', 'P02BL-B211', 'P02LT10516', 'P02LT10506', 'P02LT-1612', 'P02LT-1610', 'P02LT-1608', 'P02LT-1607', 'P02LT-1605', 'P02LT-1515', 'P02LT-1510', 'P02LT-1509', 'P02LT-1814', 'P02LT-1809', 'P02LT10404', 'P02LT10304', 'P02LT10303', 'P02LT-4116', 'P02LT-4110', 'P02LT-4108', 'P02LT-4101', 'P02LT-2503', 'P02LT-2316', 'P02LT-2314', 'P02LT-2313', 'P02LT-2306', 'P02LT-2305', 'P02LT-2304', 'P02LT-2302', 'P02LT-2216', 'P02LT-2215', 'P02LT-2214', 'P02LT-2213', 'P02LT-2212', 'P02LT-2211', 'P02LT-2210', 'P02LT-2208', 'P02LG-2016', 'P02LT-2206', 'P02LT-2203', 'P02LT-2202', 'P02LT-2201', 'P02LT-0514', 'P02LT-0510', 'P02LT-0502', 'P02LT-4903', 'P02LT-0914', 'P02LT-0905', 'P02LT-0901', 'P02LT-0807', 'P02LT-3507', 'P02LT-3505', 'P02LT-3312', 'P02LT-3116', 'P02LT-3103', 'P02LT-3813', 'P02LT-3806', 'P02LT-3803', 'P02LT-0708', 'P02LT-0704', 'P02LT-0703', 'P02LT-0701', 'P02LT-0615', 'P02LT-0614', 'P02LT-0610']
clade_3 = ['P02SL-0405', 'P02LT-1006', 'P02LG-0116', 'P02LG-0104', 'P02LG-0118', 'P02SL-0310', 'P02SL-0421', 'P02BLB3e16', 'P02BLB3e14', 'P02BLB3e12', 'P02BLB3e11', 'P02BLB3e10', 'P02BLB3e09', 'P02BLB3e08', 'P02BLB3e07', 'P02BLB3e04', 'P02LG-0115', 'P02BL-B310', 'P02BL-B302', 'P02LT-1002', 'P02LG-0105', 'P02SL-0318', 'P02LT-1015', 'P02LT-0706', 'P02LT-0504', 'P02SL-0414', 'P02SL-0324', 'P02SL-0417', 'P02LT-1014', 'P02LG-0111', 'P02LT-0904', 'P02SL-0407', 'P02SL-0315', 'P02LT-4215', 'P02LT-4210', 'P02LT-0604', 'P02SL-0416', 'P02LT-1806', 'P02LT-0616', 'P02BL-B521', 'P02LT-1005', 'P02LT-0707', 'P02LT-1010', 'P02LT-1004', 'P02LT10505', 'P02SL-0424', 'P02LT-1013', 'P02LG-0106', 'P02LT10305', 'P02LT-0505', 'P02LT-0613', 'P02LT-3102', 'P02LT-1504', 'P02LT-0702', 'P02LG-0114', 'P02LT-3614', 'P02SL-0301', 'P02LG-0113', 'P02SL-0418', 'P02SL-0412', 'P02BL-B509', 'P02LT-1508', 'P02LT-0513', 'P02LT-0503', 'P02LT-1016', 'P02LT-1012', 'P02LT-1011', 'P02LT-1009', 'P02LT-1001', 'P02SL-0303', 'P02LT-2402', 'P02LT-2415']
clade_4 = ['P02LT-1405', 'P02LT-3011', 'P02LT-1611', 'P02LT10512', 'P02LT-0501', 'P02LT-4913', 'P02LT-4912', 'P02LT-1503', 'P02LT-1505', 'P02LT-1701', 'P02LT-1615', 'P02LT-4911', 'P02LT-1711', 'P02LT-1710', 'P02LT-1812', 'P02LT-1409', 'P02LG-0109', 'P02LT-4916', 'P02LT-1614', 'P02LT-1506', 'P02LG-2014', 'P02LG-2010', 'P02LG-2008', 'P02LG-2004', 'P02LT-0507', 'P02LT-3807', 'P02LT-3804', 'P02LG-2019', 'P02SL-0409', 'P02BL-B224', 'P02BL-B218', 'P02BL-B212', 'P02LT-2116', 'P02LT-2115', 'P02LT-2114', 'P02LT-2113', 'P02LT-2112', 'P02LT-2111', 'P02LT-2110', 'P02LT-2108', 'P02LT-2107', 'P02LT-2106', 'P02LT-2104', 'P02LT-2103', 'P02LT-2101', 'P02LT-1911', 'P02LT-4115', 'P02LT-2513', 'P02LT-2510', 'P02LT-2315', 'P02LT-2105', 'P02LT-2307', 'P02BL-B217', 'P02BL-B216', 'P02BL-B207', 'P02BL-B203', 'P02LT-2310', 'P02LT-2204', 'P02LT-1314', 'P02LT-3502', 'P02LT-3315', 'P02LT-3307', 'P02LT-3305', 'P02LT-2102', 'P02LT-2207', 'P02LT-2914', 'P02LT-0814', 'P02LT-0810', 'P02LT-4902', 'P02LT-4905', 'P02LT-4901', 'P02LT-2703', 'P02LT-2311', 'P02LT-3509', 'P02LT-3508', 'P02LT-3503', 'P02LT-2702', 'P02LT-2706', 'P02LT-2411', 'P02LT-2714', 'P02LT-2712', 'P02LT-4910', 'P02LT10409', 'P02LT-2515', 'P02BL-B807', 'P02BL-B811', 'P02BL-B809', 'P02LT-1713', 'P02LT-1706', 'P02LT-1603', 'P02BL-B813', 'P02LT-2504', 'P02BL-B808', 'P02BL-B810', 'P02SL-0410', 'P02LG-0122', 'P02LT-4915']
clade_5 = ['P02LT-1408', 'P02LT-3710', 'P02LT-1410', 'P02LT-3605', 'P02LT-3604', 'P02BLB3e13', 'P02BLB3e06', 'P02BLB3e05', 'P02BLB3e03', 'P02BLB3e02', 'P02BLB3e01', 'P02SL-0307', 'P02LG-2023', 'P02LG-2018', 'P02LG-2005', 'P02SL-0408', 'P02SL-0422', 'P02SL-0403', 'P02SL-0402', 'P02BL-B311', 'P02BL-B307', 'P02BL-B304', 'P02LT10515', 'P02LT10508', 'P02LT-1913', 'P02LT10403', 'P02LT10302', 'P02LT-4109', 'P02LT-4102', 'P02LT-3613', 'P02LT-1914', 'P02LT10416', 'P02LT-4104', 'P02LT-2805', 'P02LT-2707', 'P02LT-2701', 'P02LT-2516', 'P02LT-1702', 'P02LT-3603', 'P02LT-1712', 'P02LT-1704', 'P02LT-1606', 'P02LT-1601', 'P02LT-1803', 'P02LT-2501', 'P02LT-2616', 'P02LT-2614', 'P02LT-2607', 'P02LT-2606', 'P02LT-4914', 'P02LT-4909', 'P02LT-4907', 'P02LT-4906', 'P02LT-4904', 'P02LT-1311', 'P02LT-1308', 'P02LT-1305', 'P02LT-1302', 'P02LT-1301', 'P02SL-0304', 'P02LT-0509', 'P02LT-3313', 'P02LT-3611', 'P02LT-4212', 'P02LT-4209', 'P02LT-4208', 'P02LT-3906', 'P02LT-1309', 'P02LT-3311', 'P02LT-4214', 'P02LT-4211', 'P02LT-4203', 'P02BL-B309', 'P02LT-3916', 'P02LT-3915', 'P02LT-3914', 'P02LT-3913', 'P02LT-3902', 'P02LT-3309', 'P02LT-3113', 'P02LT10408', 'P02LT-3104', 'P02LT-3101', 'P02LT-2912', 'P02LT-2911', 'P02LT-2903', 'P02LT-2909', 'P02LT-2908', 'P02LT-2906', 'P02LT-2902', 'P02LT-3815', 'P02LT-3814', 'P02LT-3811', 'P02LT-3006', 'P02BL-B416', 'P02LT-1916', 'P02LT-3802', 'P02LT-3716', 'P02LT-3910', 'P02LT-3907', 'P02LT-4908', 'P02LT-3001', 'P02LT10513', 'P02LG-2013', 'P02LT-1816', 'P02LT10415', 'P02LT-3809', 'P02LT-3713', 'P02LT-1415', 'P02BLB4e08', 'P02LG-0123', 'P02LT-1715', 'P02BL-B419', 'P02BL-B418', 'P02BL-B417', 'P02BL-B413', 'P02BL-B405', 'P02LT-1910', 'P02LT-4105', 'P02SL-0320', 'P02LT-1609', 'P02LT-1802', 'P02LT-2502', 'P02LT-2406', 'P02LT-3115', 'P02LT-1414', 'P02LT-3609', 'P02LT-3606', 'P02SL-0322', 'P02LT-1705', 'P02LT-1703', 'P02BL-B524', 'P02BL-B523', 'P02BL-B522', 'P02BL-B520', 'P02BL-B519', 'P02BL-B518', 'P02BL-B517', 'P02BL-B516', 'P02BL-B514', 'P02BL-B513', 'P02BL-B512', 'P02BL-B511', 'P02BL-B510', 'P02BL-B508', 'P02BL-B507', 'P02BL-B501', 'P02BL-B506', 'P02BL-B505', 'P02BL-B504', 'P02BL-B502', 'P02LT-1604', 'P02LT-4107', 'P02LT-3712', 'P02BL-B804', 'P02LG-2024', 'P02LG-2009', 'P02LG-2007', 'P02LG-2002', 'P02LT-1904', 'P02LT-1810', 'P02LT10407', 'P02LT10406', 'P02LT10316', 'P02LT10315', 'P02LT10313', 'P02LT10311', 'P02LT-4111', 'P02LT-2803', 'P02LT-2603', 'P02LT-2602', 'P02LT-3706', 'P02LG-2020', 'P02LT-1915', 'P02LT-1908', 'P02LT-1906', 'P02LT-1902', 'P02LT10314', 'P02LT-2613', 'P02LT-3702', 'P02LT-1815', 'P02LT-3816', 'P02LT-3812', 'P02LT-0716', 'P02LT-1413', 'P02LT-1403', 'P02LT-4202', 'P02BL-B605', 'P02SL-0321', 'P02LG-2022', 'P02LG-2006', 'P02LG-2003', 'P02LT-1811', 'P02LT-1804', 'P02LT-0511', 'P02LT-1312', 'P02LT-1310', 'P02LT-1307', 'P02LT-1304', 'P02LT-1303', 'P02LT-0816', 'P02LT-0812', 'P02LT-0804', 'P02LT-4216', 'P02LT-4213', 'P02LT-4207', 'P02LT-4205', 'P02LT-4204', 'P02LT-4201', 'P02LT-3912', 'P02LT-3911', 'P02LT-3909', 'P02LT-3908', 'P02LT-3515', 'P02LT-3512', 'P02LT-3504', 'P02LT-3501', 'P02LT-2907', 'P02LT-1412', 'P02BL-B623', 'P02BL-B618', 'P02BL-B617', 'P02BL-B615', 'P02BL-B613', 'P02LT-1613', 'P02LT-1602', 'P02LT-1813', 'P02LT-1411', 'P02SL-0306', 'P02LT-3105', 'P02LT-3704', 'P02LT-1008', 'P02LT-0608', 'P02LG-0103', 'P02LT-0515', 'P02LT-0606', 'P02SL-0319', 'P02LT-0815', 'P02LT-0811', 'P02LT-0813', 'P02LT-0809', 'P02LT-0806', 'P02LT-0602']
clade_1 = ["P02LT-0913", "P02LT-0907", "P02LT-2512", "P02LT-0506", "P02LT-0912", "P02LT-0906", "P02LT-2312", "P02LT-2301", "P02LT-2916", "P02LT-3007", "P02LT-2604"]
clade_4_5_outliers = ["P02LT-3306", "P02LT-3610", "P02LT-2704", "P02LT-1805", "P02LT-1407", "P02LT-1514", "P02LT-1116", "P02LT-1109", "P02LT-1007", "P02LT-1003", "P02LT-2915"]
clade_Q_R = ["P02LG-2011", "P02LT-3008","P02LT-3801","P02SL-0309"]

blod_time_dictionary= {"BL-B2": "-8", "BLB2e": "-7", "B3": "-6", "B4": "-5", "B5": "-3", "B6": "-2", "B7": "-1", "B8": "0"};
#note BL-B2 is actually -7, but its coded as -8 for plotting purposes

O_ant_neg = clade_1
O_ant_long = clade_2 + clade_3
O_ant_med = clade_4_5_outliers + clade_5 + clade_4 + clade_Q_R
O_ant_short = clade_6

LPS_to_analyze= [[O_ant_neg],[O_ant_long],[O_ant_med],[O_ant_short]]

lung_loc=pd.read_excel("isolate_metadata_table.xlsx",sheet_name='S2 Isolate Metadata',skiprows=[0,1])
lung_locations = lung_loc['Source Label']
lung_lobes = lung_loc['Lung Lobe/Side']

samples_for_analysis = clade_1 + clade_2 + clade_3 + clade_4 + clade_5 + clade_6 + clade_4_5_outliers + clade_Q_R

unique_site = np.asarray(list(set([x[3:len(x)-2] for x in samples_for_analysis])))

removed_sites = ['BL','SL','LG']
removed_sites_blood_issue = ['SL','LG']

unique_site_bl_lg_sl = [x for x in unique_site if any(y in x for y in removed_sites_blood_issue)]
unique_site_bl_lg_sl = unique_site_bl_lg_sl + ["BLB2e", "BL-B2", "B3", "B4", "B5", "B6","B7", "B8"]
# blood sites are hardcoded in to deal with extra blood plate issue

unique_site = [x for x in unique_site if not any(y in x for y in removed_sites)]

redone_unique_site_nums = copy.deepcopy(unique_site)
for item in ['LT-0','LT-','LT','SL-0','LG-0','LG-']:
    redone_unique_site_nums = [x.replace(item, '') for x in redone_unique_site_nums]
    
redone_unique_site_nums = [int(x) for x in redone_unique_site_nums]
unique_site_total = unique_site + unique_site_bl_lg_sl


#######


observed_LPS_counts = []
for ii in range(0, len(unique_site_total)):
    obs_calc = [sum(unique_site_total[ii] in s for s in x[0]) for x in LPS_to_analyze]
    observed_LPS_counts.append(obs_calc)
    
    
observed_LPS_counts = np.asarray(observed_LPS_counts)
observed_LPS_counts_freq = observed_LPS_counts.T/np.sum(observed_LPS_counts,1)
    
#############

lung_loc=pd.read_excel("combined_table_for_reviewers_final_final_meta_sunday.xlsx",sheet_name='S2 Isolate Metadata',skiprows=[0,1])
lung_loc = lung_loc.iloc[0:42,:]
lung_locations = lung_loc['Source Label'].astype(int)
lung_lobes = lung_loc['Lung Lobe/Side']

lobes = ['LLL', 'LUL', 'RLL', 'RML', 'RUL']

sites_in_lobe = [list(lung_locations[lung_lobes==x]) for x in lobes]
index_of_site_by_lobe = [[redone_unique_site_nums.index(x) for x in y] for y in sites_in_lobe]


# Add on blood, spleen, and lymph node catagories 
body_sites = ["B","LG","SL"]
body_site_indexes = [[i for i in range(len(unique_site_total)) if y in unique_site_total[i]] for y in body_sites]
body_site_nums = [[unique_site_total[i] for i in range(len(unique_site_total)) if y in unique_site_total[i]] for y in body_sites]



index_of_site_by_lobe = index_of_site_by_lobe + body_site_indexes
sites_in_lobe = sites_in_lobe + body_site_nums

# Reorder blood by date
blood_dates = [[val for key,val in blod_time_dictionary.items() if key in x] for x in  sites_in_lobe[5]]
blood_dates = [int(x[0]) for x in blood_dates]
sorted_blood_dates = np.argsort(blood_dates)

sites_in_lobe[5] = [str(x) for x in np.asarray(blood_dates)[sorted_blood_dates]]
sites_in_lobe[5][sites_in_lobe[5].index("-7")] = "-7b"
sites_in_lobe[5][sites_in_lobe[5].index("-8")] = "-7a"

index_of_site_by_lobe[5] = [index_of_site_by_lobe[5][x] for x in sorted_blood_dates]



colors=  ['#1f78b4','#a6cee3','#33a02c','#b2df8a']
hfont = {'fontname':'Helvetica'}





#############################################
fig_size = (25, 10)
font_size_set = 16
ax = plt.figure(figsize=fig_size,dpi=300)

# Stacked bar chart with loop
all_x_vals = []
site_iso_counts = []
counter = 3
for ij in range(0,len(sites_in_lobe)):
    x_vals = list(range(counter, len(sites_in_lobe[ij])+counter))
    if ij<len(sites_in_lobe)-4:
        counter = max(x_vals) + 2
    else: 
        counter = max(x_vals) + 4
    indexes = index_of_site_by_lobe[ij]
    all_x_vals += x_vals
    index_LPS_counts = np.take(observed_LPS_counts_freq,indexes,1)
    num_LPS_counts = np.take(observed_LPS_counts.T,indexes,1)

    site_iso_counts.append(sum(num_LPS_counts))
    
    for i in range(0,4):
        plt.bar(x_vals, index_LPS_counts[i], bottom = np.sum(index_LPS_counts[:i], axis = 0),color=colors[i])
      
plt.yticks(fontsize = 30) 

site_labels = []
for j in range(0,len(sites_in_lobe)):
    site_labels += [str(x) for x in sites_in_lobe[j]]
    #site_labels += [str(sites_in_lobe[j][x]) + "\n" + str(site_iso_counts[j][x]) for x in range(0, len(sites_in_lobe[j]))]



new_site_replaceables = ["LG-","LG","SL-","SL"]
for y in new_site_replaceables:
    site_labels = [x.replace(y,"") for x in site_labels]

plt.xticks(ticks=all_x_vals, labels=site_labels,fontsize=font_size_set,**hfont,rotation=360-45,fontweight='heavy')
plt.tick_params(axis="x",direction="out", pad=20,which='major')

iso_flatten = [item for row in site_iso_counts for item in row]
plt.xticks(ticks=np.asarray(all_x_vals)+.001, labels=iso_flatten,fontsize=font_size_set,minor=True)
plt.tick_params(axis="x",direction="out", pad=55,which='minor')

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=True) # labels along the bottom edge are off

for axis in ['bottom','left']:
    plt.gca().spines[axis].set_linewidth(4)
for axis in ['top','right']:
    plt.gca().spines[axis].set_linewidth(0)
    


##############################################
### combine counts by lobe
lobe_comp_freq = []
lobe_comp = []
for x in index_of_site_by_lobe:
    index_LPS_counts = np.take(observed_LPS_counts.T,x,1)
    norm_counts = np.sum(index_LPS_counts,1)/np.sum(index_LPS_counts)
    unnorm_counts= np.sum(index_LPS_counts,1)
    lobe_comp_freq.append(norm_counts)
    lobe_comp.append(unnorm_counts)
    
LPS_freq_by_lobe = np.vstack(lobe_comp_freq).T

fig = plt.figure(figsize=fig_size,dpi=300)
for i in range(0,4):
    plt.bar(range(0,8), LPS_freq_by_lobe[i], bottom = np.sum(LPS_freq_by_lobe[:i], axis = 0),color=colors[i])

plt.yticks(fontsize = 30) 

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=True) # labels along the bottom edge are off

for axis in ['bottom','left']:
    plt.gca().spines[axis].set_linewidth(4)
for axis in ['top','right']:
    plt.gca().spines[axis].set_linewidth(0)

    
plt.xticks(ticks=range(0,8), labels=[sum(x) for x in lobe_comp],fontsize=30,**hfont)
plt.tick_params(axis="x",direction="out", pad=20)

plt.show()

lobe_comp_for_chi = np.asarray(lobe_comp)[0:5,:]
lobe_composition_for_chi = np.asarray(lobe_comp_for_chi).squeeze().T
lobe_O_ant_chi_square = chi2_contingency(lobe_composition_for_chi)







