\ENCODING=ISO-8859-1
\Problem name: ilog.cplex

Minimize
 obj: CmaxUpper
Subject To
 constraint1_2_1#0:    X11 + X10  = 1
 constraint1_2_2#1:    X21 + X20  = 1
 constraint1_2_4#2:    X41 + X40  = 1
 constraint1_2_5#3:    X51 + X50  = 1
 constraint1_2_6#4:    X61 + X60  = 1
 constraint1_2_9#5:    X91 + X90  = 1
 constraint1_3_0,0#6:  alpha0 + 2999 X00 <= 3000
 constraint1_3_0,1#7:  alpha1 + 2999 X01 <= 3000
 constraint1_3_1,0#8:  2998 X10 + alpha0 <= 3000
 constraint1_3_1,1#9:  2998 X11 + alpha1 <= 3000
 constraint1_3_2,0#10: 2997 X20 + alpha0 <= 3000
 constraint1_3_2,1#11: 2997 X21 + alpha1 <= 3000
 constraint1_3_3,0#12: alpha0 + 2996 X30 <= 3000
 constraint1_3_3,1#13: alpha1 + 2996 X31 <= 3000
 constraint1_3_4,0#14: 2995 X40 + alpha0 <= 3000
 constraint1_3_4,1#15: 2995 X41 + alpha1 <= 3000
 constraint1_3_5,0#16: 2994 X50 + alpha0 <= 3000
 constraint1_3_5,1#17: 2994 X51 + alpha1 <= 3000
 constraint1_3_6,0#18: 2993 X60 + alpha0 <= 3000
 constraint1_3_6,1#19: 2993 X61 + alpha1 <= 3000
 constraint1_3_7,0#20: alpha0 + 2992 X70 <= 3000
 constraint1_3_7,1#21: alpha1 + 2992 X71 <= 3000
 constraint1_3_8,0#22: alpha0 + 2991 X80 <= 3000
 constraint1_3_8,1#23: alpha1 + 2991 X81 <= 3000
 constraint1_3_9,0#24: 2990 X90 + alpha0 <= 3000
 constraint1_3_9,1#25: 2990 X91 + alpha1 <= 3000
 constraint1_4_0,0#26: - 3001 X00 + beta0 >= -3000
 constraint1_4_0,1#27: - 3001 X01 + beta1 >= -3000
 constraint1_4_1,0#28: - 3002 X10 + beta0 >= -3000
 constraint1_4_1,1#29: - 3002 X11 + beta1 >= -3000
 constraint1_4_2,0#30: - 3003 X20 + beta0 >= -3000
 constraint1_4_2,1#31: - 3003 X21 + beta1 >= -3000
 constraint1_4_3,0#32: - 3004 X30 + beta0 >= -3000
 constraint1_4_3,1#33: - 3004 X31 + beta1 >= -3000
 constraint1_4_4,0#34: - 3005 X40 + beta0 >= -3000
 constraint1_4_4,1#35: - 3005 X41 + beta1 >= -3000
 constraint1_4_5,0#36: - 3006 X50 + beta0 >= -3000
 constraint1_4_5,1#37: - 3006 X51 + beta1 >= -3000
 constraint1_4_6,0#38: - 3007 X60 + beta0 >= -3000
 constraint1_4_6,1#39: - 3007 X61 + beta1 >= -3000
 constraint1_4_7,0#40: - 3008 X70 + beta0 >= -3000
 constraint1_4_7,1#41: - 3008 X71 + beta1 >= -3000
 constraint1_4_8,0#42: - 3009 X80 + beta0 >= -3000
 constraint1_4_8,1#43: - 3009 X81 + beta1 >= -3000
 constraint1_4_9,0#44: - 3010 X90 + beta0 >= -3000
 constraint1_4_9,1#45: - 3010 X91 + beta1 >= -3000
 constraint1_8_0#46:   alpha0 - alpha1 <= -2
 constraint1_9_0#47:   beta0 - beta1 <= -2
Bounds
 0 <= X11 <= 1
 0 <= X10 <= 1
 0 <= X21 <= 1
 0 <= X20 <= 1
 0 <= X41 <= 1
 0 <= X40 <= 1
 0 <= X51 <= 1
 0 <= X50 <= 1
 0 <= X61 <= 1
 0 <= X60 <= 1
 0 <= X91 <= 1
 0 <= X90 <= 1
 1 <= alpha0 <= 10
 0 <= X00 <= 1
 1 <= alpha1 <= 10
 0 <= X01 <= 1
 0 <= X30 <= 1
 0 <= X31 <= 1
 0 <= X70 <= 1
 0 <= X71 <= 1
 0 <= X80 <= 1
 0 <= X81 <= 1
 1 <= beta0 <= 10
 1 <= beta1 <= 10
Generals
 X11  X10  X21  X20  X41  X40  X51  X50  X61  X60  X91  X90  X00  X01  X30 
 X31  X70  X71  X80  X81 
End
