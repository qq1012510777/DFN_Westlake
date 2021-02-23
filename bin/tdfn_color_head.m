clc;
close all;
clear all;

JXY_1=[-9.86, -11.896, 12.467; -9.896, -17.999, 18; -9.731, -17.999, 10.8; -9.602, -10.152, -0.002; -9.567, -17.999, 3.6; -9.403, -17.999, -3.599; -9.993, -1.55, 11.506; -10.092, -4.841, 18; -9.824, -8.127, 8.387; -9.714, -13.648, 7.198; -10.066, 7.98, 8.491; -10.179, 16.964, 7.556; -10.283, 15.929, 12.778; -9.912, 11.647, -0.654; -9.823, 1.857, 1.853; -9.721, 7.974, -6.637; -9.84, -2.795, 5.639; -10.288, 8.315, 18; -10.386, 14.893, 18; -9.963, 2.453, 7.589; -9.652, 2.358, -5.986; -9.533, 2.961, -11.606; -9.939, 18, -3.648; -10.075, 18, 2.334; -9.802, 18, -9.632; -9.548, 11.979, -16.807; -9.665, 18, -15.616; -9.601, 8.588, -12.268; -9.431, 5.958, -17.999; -9.683, -3.926, -0.502; -9.163, -12.009, -18; -9.238, -17.999, -10.799; -9.074, -17.999, -17.999; -9.392, -13.18, -7.201; -9.532, -3.029, -7.696; -9.355, -9.017, -11.542; -9.253, -6.02, -18; -9.342, -0.03, -17.999; -10.133, 5.027, 13.356; -9.994, -11.42, 18; -10.19, 1.736, 18; -9.878, -14.948, 15.234; -9.814, -17.999, 14.4; -9.796, -14.948, 11.634; -9.584, -14.076, 1.798; -9.485, -17.999, 0; -9.502, -14.076, -1.8; -10.042, -3.196, 14.753; -9.958, -6.484, 13.194; -9.908, -4.839, 9.947; -9.723, -15.823, 8.999; -9.787, -12.772, 9.833; -9.658, -11.9, 3.598; -9.713, -9.14, 4.192; -9.769, -10.888, 7.793; -10.123, 12.472, 8.024; -10.231, 16.447, 10.167; -10.174, 11.954, 10.635; -9.868, 6.752, 0.599; -9.772, 4.915, -2.391; -9.816, 9.81, -3.646; -9.917, -2.172, 8.572; -9.832, -5.461, 7.013; -10.285, 12.122, 15.389; -10.334, 15.411, 15.389; -10.337, 11.604, 18; -9.902, -0.17, 6.614; -9.832, -0.468, 3.746; -9.893, 2.155, 4.721; -9.592, 2.66, -8.796; -9.627, 5.467, -9.122; -9.686, 5.166, -6.312; -9.925, 14.824, -2.151; -10.007, 18, -0.657; -9.994, 14.824, 0.839; -9.738, 2.108, -2.066; -9.675, 14.99, -13.22; -9.607, 14.99, -16.212; -9.734, 18, -12.624; -9.989, 9.814, 3.918; -9.945, 4.919, 5.172; -9.516, 7.273, -15.134; -9.49, 8.969, -17.403; -9.575, 10.284, -14.538; -9.762, -3.36, 2.568; -9.721, -6.474, 2.818; -9.643, -7.039, -0.252; -9.201, -15.004, -14.399; -9.156, -17.999, -14.399; -9.119, -15.004, -18; -9.462, -8.105, -7.448; -9.567, -6.591, -3.849; -9.497, -11.666, -3.601; -9.304, -7.518, -14.771; -9.393, -4.524, -12.848; -9.444, -6.023, -9.619; -9.567, 5.775, -11.937; -9.482, 4.46, -14.803; -10.045, 14.306, 3.45; -9.374, -11.099, -9.371; -9.297, -3.025, -18; -9.437, -1.53, -12.848; -9.702, 13.294, -10.95; -9.761, 12.987, -8.134; -9.661, 8.281, -9.452; -9.87, 18, -6.64; -9.83, 12.987, -5.142; -10.177, 8.147, 13.245; -10.1, 6.503, 10.924; -10.211, 6.671, 15.678; -9.976, -8.369, 15.234; -10.043, -8.131, 18; -9.927, -11.658, 15.234; -9.842, -10.012, 10.427; -9.315, -15.59, -9; -9.297, -13.508, -11.171; -9.608, -3.477, -4.099; -9.259, -10.513, -14.771; -9.208, -9.014, -18; -9.533, -0.033, -9.651; -9.437, 1.465, -14.803; -9.387, 2.963, -17.999; -10.127, 17.482, 4.945; -10.048, 3.74, 10.473; -10.063, 1.738, 12.431; -9.978, 0.451, 9.548; -10.091, 0.092, 14.753; -10.141, -1.552, 18; -10.239, 5.025, 18; -10.162, 3.381, 15.678; -9.398, -15.59, -5.4; -9.321, -17.999, -7.199; -10.015, 5.216, 8.04; -9.668, -0.783, -3.244; -9.592, -0.335, -6.841; -9.753, -1.034, 0.675; -9.649, -17.999, 7.2; -9.641, -15.823, 5.399; -9.945, -14.709, 18; ];
JXY_2D_1=[-10.596, -8.547, 0; -16.038, -14.73, 0; -8.837, -14.625, 0; 1.85, -6.621, 0; -1.636, -14.52, 0; 5.564, -14.415, 0; -9.788, 1.813, 0; -16.234, -1.572, 0; -6.571, -4.718, 0; -5.299, -10.221, 0; -6.915, 11.387, 0; -6.114, 20.386, 0; -11.321, 19.274, 0; 2.177, 15.188, 0; -0.184, 5.362, 0; 8.215, 11.602, 0; -3.902, 0.653, 0; -16.43, 11.584, 0; -16.528, 18.162, 0; -5.931, 5.874, 0; 7.649, 5.977, 0; 13.261, 6.662, 0; 5.077, 21.585, 0; -0.907, 21.497, 0; 11.061, 21.672, 0; 18.328, 15.756, 0; 17.046, 21.759, 0; 13.838, 12.299, 0; 19.609, 9.752, 0; 2.257, -0.387, 0; 19.877, -8.215, 0; 12.765, -14.31, 0; 19.967, -14.205, 0; 9.095, -9.544, 0; 9.439, 0.614, 0; 13.375, -5.317, 0; 19.788, -2.226, 0; 19.699, 3.763, 0; -11.737, 8.364, 0; -16.136, -8.151, 0; -16.332, 5.005, 0; -13.317, -11.638, 0; -12.437, -14.677, 0; -9.716, -11.586, 0; 0.107, -10.57, 0; 1.964, -14.467, 0; 3.707, -10.518, 0; -13.011, 0.12, 0; -11.402, -3.145, 0; -8.18, -1.452, 0; -7.068, -12.423, 0; -7.948, -9.384, 0; -1.724, -8.421, 0; -2.36, -5.67, 0; -5.935, -7.47, 0; -6.514, 15.887, 0; -8.717, 19.83, 0; -9.118, 15.331, 0; 0.996, 10.275, 0; 4.015, 8.482, 0; 5.196, 13.395, 0; -6.845, 1.233, 0; -5.236, -2.032, 0; -13.876, 15.429, 0; -13.925, 18.718, 0; -16.479, 14.873, 0; -4.916, 3.264, 0; -2.043, 3.008, 0; -3.058, 5.618, 0; 10.454, 6.319, 0; 10.738, 9.132, 0; 7.932, 8.789, 0; 3.627, 18.386, 0; 2.085, 21.541, 0; 0.635, 18.343, 0; 3.732, 5.669, 0; 14.695, 18.714, 0; 17.687, 18.757, 0; 14.054, 21.715, 0; -2.369, 13.288, 0; -3.55, 8.375, 0; 16.724, 11.025, 0; 18.968, 12.754, 0; 16.083, 14.027, 0; -0.822, 0.133, 0; -1.026, -2.983, 0; 2.053, -3.504, 0; 16.321, -11.263, 0; 16.366, -14.257, 0; 19.922, -11.21, 0; 9.267, -4.465, 0; 5.644, -3.003, 0; 5.472, -8.083, 0; 16.581, -3.771, 0; 14.614, -0.805, 0; 11.407, -2.351, 0; 13.549, 9.48, 0; 16.435, 8.207, 0; -1.968, 17.787, 0; 11.235, -7.43, 0; 19.743, 0.768, 0; 14.569, 2.188, 0; 12.45, 16.985, 0; 9.638, 16.637, 0; 11.027, 11.95, 0; 8.069, 21.628, 0; 6.646, 16.593, 0; -11.673, 11.486, 0; -9.326, 9.876, 0; -14.084, 9.974, 0; -13.415, -5.059, 0; -16.185, -4.862, 0; -13.366, -8.349, 0; -8.583, -6.632, 0; 10.93, -11.927, 0; 13.07, -9.813, 0; 5.848, 0.113, 0; 16.626, -6.766, 0; 19.833, -5.221, 0; 11.35, 3.638, 0; 16.48, 5.212, 0; 19.654, 6.757, 0; -3.51, 20.942, 0; -8.834, 7.119, 0; -10.763, 5.088, 0; -7.859, 3.843, 0; -13.06, 3.409, 0; -16.283, 1.716, 0; -16.381, 8.294, 0; -14.035, 6.684, 0; 7.33, -11.979, 0; 9.165, -14.362, 0; -6.423, 8.631, 0; 4.953, 2.794, 0; 8.544, 3.295, 0; 1.036, 2.487, 0; -5.236, -14.572, 0; -3.467, -12.371, 0; -16.087, -11.44, 0; ];
JM_1=[1, 42, 2, 43, 3, 44; 4, 45, 5, 46, 6, 47; 7, 48, 8, 49, 9, 50; 1, 44, 3, 51, 10, 52; 10, 53, 4, 54, 9, 55; 11, 56, 12, 57, 13, 58; 14, 59, 15, 60, 16, 61; 17, 62, 7, 50, 9, 63; 18, 64, 13, 65, 19, 66; 20, 67, 17, 68, 15, 69; 21, 70, 22, 71, 16, 72; 14, 73, 23, 74, 24, 75; 21, 72, 16, 60, 15, 76; 25, 77, 26, 78, 27, 79; 14, 80, 11, 81, 15, 59; 28, 82, 29, 83, 26, 84; 30, 85, 17, 86, 4, 87; 31, 88, 32, 89, 33, 90; 34, 91, 35, 92, 4, 93; 36, 94, 37, 95, 35, 96; 28, 97, 22, 98, 29, 82; 14, 99, 12, 56, 11, 80; 34, 100, 36, 96, 35, 91; 35, 95, 37, 101, 38, 102; 28, 103, 25, 104, 16, 105; 16, 104, 25, 106, 23, 107; 11, 58, 13, 64, 18, 108; 39, 109, 11, 108, 18, 110; 1, 111, 8, 112, 40, 113; 1, 52, 10, 55, 9, 114; 34, 115, 32, 116, 36, 100; 30, 87, 4, 92, 35, 117; 36, 116, 32, 88, 31, 118; 36, 118, 31, 119, 37, 94; 22, 120, 35, 102, 38, 121; 22, 121, 38, 122, 29, 98; 14, 61, 16, 107, 23, 73; 14, 75, 24, 123, 12, 99; 20, 124, 39, 125, 7, 126; 7, 127, 41, 128, 8, 48; 39, 110, 18, 129, 41, 130; 39, 130, 41, 127, 7, 125; 34, 93, 4, 47, 6, 131; 34, 131, 6, 132, 32, 115; 17, 63, 9, 54, 4, 86; 20, 69, 15, 81, 11, 133; 21, 134, 30, 117, 35, 135; 30, 136, 15, 68, 17, 85; 20, 133, 11, 109, 39, 124; 20, 126, 7, 62, 17, 67; 28, 84, 26, 77, 25, 103; 28, 105, 16, 71, 22, 97; 21, 76, 15, 136, 30, 134; 21, 135, 35, 120, 22, 70; 10, 51, 3, 137, 5, 138; 10, 138, 5, 45, 4, 53; 1, 114, 9, 49, 8, 111; 1, 113, 40, 139, 2, 42; ];
Data_1=[88.1156; 100; 83.9153; 60.9367; 69.1366; 52.4831; 86.1068; 100; 79.6252; 76.5243; 80.1798; 76.2276; 87.6294; 62.1909; 65.8369; 51.0466; 73.6497; 100; 100; 77.8742; 50.1558; 37.9697; 57.9106; 66.2865; 48.4481; 40.9389; 45.1523; 41.5744; 20; 60.3517; 20; 36.8339; 20; 44.9532; 44.1303; 35.2232; 20; 20; 89.3815; 100; 100; 92.9553; 91.6198; 86.4152; 64.2815; 60.8097; 57.4646; 91.8897; 89.3424; 83.2508; 79.7095; 83.3161; 69.9695; 70.7463; 77.7814; 78.7512; 81.9287; 83.3581; 63.7857; 57.848; 57.5256; 79.0749; 76.0581; 94.8661; 94.0466; 100; 76.0675; 69.8; 72.8466; 45.1234; 44.6703; 50.4589; 60.415; 62.1505; 64.0174; 56.6635; 45.0268; 43.168; 46.8685; 69.5577; 73.1178; 33.1244; 31.6991; 41.2655; 65.9817; 68.2375; 60.7698; 27.1142; 27.8703; 20; 44.5876; 52.6771; 51.5874; 28.9795; 33.5931; 39.318; 40.0853; 30.4314; 68.0265; 40.9542; 20; 33.6472; 44.7561; 50.1888; 45.493; 53.1662; 53.957; 88.5786; 84.9896; 93.8078; 93.1089; 100; 93.062; 84.59; 41.5615; 35.6632; 53.1427; 28.8837; 20; 40.9781; 30.5282; 20; 70.9827; 82.6576; 87.8358; 81.2011; 92.1172; 100; 100; 93.6626; 48.1627; 44.7391; 78.9463; 54.486; 47.8378; 62.8799; 76.5391; 73.3569; 100; ];
figure(1)
P_1 = patch('Vertices', JXY_1, 'Faces', JM_1, 'FaceVertexCData', Data_1, 'FaceColor', 'interp', 'EdgeAlpha', 1);
hold on;
view(3);
colorbar;
hold on;
figure(1)
plot3([-18 18],[-18 -18],[18 18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([18 18],[-18 18],[18 18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([18 -18],[18 18],[18 18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([-18 -18],[18 -18],[18 18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([-18 18],[-18 -18],[-18 -18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([18 18],[-18 18],[-18 -18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([18 -18],[18 18],[-18 -18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([-18 -18],[18 -18],[-18 -18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([-18 18],[-18 -18],[-18 -18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([18 18],[-18 -18],[-18 18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([18 -18],[-18 -18],[18 18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([-18 -18],[-18 -18],[18 -18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([-18 18],[18 18],[-18 -18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([18 18],[18 18],[-18 18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([18 -18],[18 18],[18 18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([-18 -18],[18 18],[18 -18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([-18 -18],[-18 18],[-18 -18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([-18 -18],[18 18],[-18 18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([-18 -18],[18 -18],[18 18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([-18 -18],[-18 -18],[18 -18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([18 18],[-18 18],[-18 -18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([18 18],[18 18],[-18 18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([18 18],[18 -18],[18 18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
plot3([18 18],[-18 -18],[18 -18],'color',[1 0 0],'Linewidth',3);
grid on; hold on;
axis([-18 18 -18 18 -18 18])
hold on;
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
title('DFN head distribution');
hold on;

%%*****flux rate vector***
figure(1)
q_vector_1=[-9.83022, -15.9656, 13.7556, 0.0444088, 0.0139031, -1.95478; -9.5249, -15.3844, -0.00100535, 0.0466283, -0.0161449, -2.03238; -9.97046, -4.84067, 12.631, 0.0429494, 0.070166, -1.9276; -9.76967, -14.5154, 10.1552, 0.0582705, -0.38778, -2.29963; -9.71418, -10.6436, 5.19429, 0.0447279, 0.0559798, -1.99625; -10.1769, 13.6244, 9.60839, 0.0402424, 0.164559, -1.87067; -9.81975, 7.15943, -1.81344, 0.0446005, -0.284525, -1.76817; -9.8866, -4.15849, 8.51073, 0.0515178, 0.190432, -2.38159; -10.32, 13.0456, 16.2592, 0.0446817, 0.166749, -2.0666; -9.87663, 0.504832, 5.02709, 0.052318, -0.249925, -2.12891; -9.63608, 4.43092, -8.07751, 0.0511091, -0.189627, -2.11535; -9.97624, 15.8824, -0.656807, 0.0286797, 0.0186262, -1.26872; -9.733, 4.06295, -3.59078, 0.0530432, -0.560202, -1.95795; -9.6729, 15.9931, -14.0193, 0.021425, -0.522476, -0.597303; -9.93488, 7.16153, 3.22952, 0.0434988, 0.103786, -1.97363; -9.52806, 8.84193, -15.6925, 0.0904729, -2.49901, -2.33101; -9.70952, -5.62543, 1.71119, 0.0499935, 0.322238, -2.40093; -9.15964, -16.0033, -15.6, 0.0453138, 0.168981, -2.09575; -9.50993, -8.78821, -4.96721, 0.0522239, -0.220067, -2.1443; -9.38121, -6.02268, -12.4135, 0.0499729, -0.198771, -2.05959; -9.5228, 5.836, -13.9589, 0.0851631, -1.19735, -2.9489; -10.0534, 12.1972, 5.13063, 0.0470226, 0.149833, -2.15811; -9.42777, -8.40958, -8.81407, 0.0490676, 0.133386, -2.23696; -9.37676, -3.02725, -14.5657, 0.058186, -0.0111397, -2.54204; -9.70905, 11.5206, -9.51326, 0.0487296, -0.379033, -1.88733; -9.82141, 14.6578, -6.63999, 0.034403, -0.216248, -1.366; -10.2133, 10.7411, 13.0896, 0.0544051, -0.123639, -2.30287; -10.1635, 7.10716, 13.2824, 0.0411453, 0.206508, -1.93764; -9.98311, -9.38696, 16.1557, 0.0543315, -0.00870673, -2.37475; -9.80039, -11.2249, 9.35102, 0.0557504, 0.037564, -2.46715; -9.32975, -13.3996, -9.84832, 0.0608178, -0.0662272, -2.62135; -9.6069, -5.70336, -2.73426, 0.0558142, -0.250857, -2.28148; -9.25341, -13.0092, -13.4478, 0.0599161, -0.39785, -2.36515; -9.25819, -9.01617, -15.8478, 0.0603089, -0.0204519, -2.62896; -9.47005, -0.0332203, -12.4349, 0.0580979, -0.43776, -2.25941; -9.43635, 2.96279, -15.8691, 0.0709028, 0.0226029, -3.12124; -9.85814, 12.5403, -3.64743, 0.0469874, 0.00270183, -2.06043; -10.0563, 15.5373, 3.0783, 0.0387975, -0.0497356, -1.66734; -10.0307, 1.97629, 10.8172, 0.0493047, 0.180882, -2.27838; -10.0926, -1.5525, 15.8351, 0.0543973, -0.0453022, -2.35372; -10.2048, 5.02589, 16.452, 0.0588517, -0.028563, -2.55982; -10.1064, 1.7373, 14.2873, 0.0433219, 0.0682738, -1.94268; -9.46668, -13.7782, -3.60147, 0.0549359, 0.233819, -2.5597; -9.34553, -16.3937, -7.20054, 0.0452254, -0.0544076, -1.94592; -9.75627, -7.02597, 4.67455, 0.0454661, -0.0191908, -1.97947; -9.9519, 4.09676, 5.97777, 0.0451821, 0.133574, -2.06685; -9.62356, -1.5329, -4.72912, 0.0465339, -0.0721047, -1.99168; -9.78335, -1.62178, 2.32968, 0.0540852, -0.313362, -2.16489; -10.0552, 5.15325, 9.81232, 0.048897, -0.393123, -1.88546; -9.93313, -0.631335, 8.24477, 0.0439748, -0.120783, -1.84775; -9.65158, 12.8558, -12.9035, 0.0199983, -0.649295, -0.451932; -9.61926, 6.50774, -10.1714, 0.048154, -0.609172, -1.71174; -9.72052, 0.0961612, -1.54581, 0.0472095, 0.0345977, -2.091; -9.57332, 0.763177, -8.43069, 0.064555, -0.630338, -2.41649; -9.67187, -16.5495, 7.19929, 0.0424158, 0.00193881, -1.85964; -9.62869, -13.9341, 3.59844, 0.0607579, -0.229213, -2.51222; -9.92635, -8.28943, 12.9515, 0.0428604, 0.0293408, -1.89702; -9.9177, -13.7727, 16.1557, 0.0553325, -0.0207741, -2.41072; ];
quiver3(q_vector_1(:,1), q_vector_1(:,2), q_vector_1(:,3), q_vector_1(:,4), q_vector_1(:,5), q_vector_1(:,6));

hold on;
%show frac pressure contour and flow rate vectors in 2D
figure(2)
x_tr1 = JXY_2D_1(:, 1);
y_tr1 = JXY_2D_1(:, 2);
z_tr1 = Data_1(:);
nx1 = 500;
ny1 = 500;
[X1,Y1] = meshgrid(linspace(min(x_tr1),max(x_tr1),nx1),linspace(min(y_tr1),max(y_tr1),ny1)) ;
Z1 =griddata(x_tr1,y_tr1,z_tr1,X1,Y1) ;
contour(X1,Y1,Z1) ;
hold on;
S_1= patch('Vertices', JXY_2D_1, 'Faces', JM_1, 'FaceVertexCData', Data_1, 'FaceColor', 'interp', 'EdgeAlpha', 0.2, 'facealpha', 0);
hold on;
q_vector_2d1=[-11.8243, -12.6347, 0, 1.95487, 0.0423951, 0; 1.92583, -11.8529, 0, 2.03294, 0.0134783, 0; -10.8653, -1.49333, 0, 1.92685, 0.0982619, 0; -8.24487, -11.1319, 0, 2.30577, -0.354261, 0; -3.34097, -7.1879, 0, 1.99572, 0.0850763, 0; -8.11767, 17.0158, 0, 1.86851, 0.191826, 0; 3.40233, 10.7173, 0, 1.77269, -0.258753, 0; -6.75463, -0.7511, 0, 2.37912, 0.225145, 0; -14.7609, 16.3401, 0, 2.06444, 0.196871, 0; -3.34, 3.963, 0, 2.13297, -0.218895, 0; 9.70807, 8.0801, 0, 2.1185, -0.158795, 0; 2.11543, 19.4234, 0, 1.26864, 0.0371185, 0; 5.22613, 7.64673, 0, 1.9666, -0.531663, 0; 15.4784, 19.7288, 0, 0.605187, -0.51377, 0; -1.64143, 10.6459, 0, 1.97239, 0.132552, 0; 17.2584, 12.6021, 0, 2.36864, -2.46504, 0; 0.0678333, -2.11893, 0, 2.39649, 0.357232, 0; 17.5365, -12.2444, 0, 2.09356, 0.199528, 0; 6.79443, -5.18437, 0, 2.14791, -0.188813, 0; 14.2007, -2.3103, 0, 2.06287, -0.168752, 0; 15.5694, 9.5709, 0, 2.96721, -1.15437, 0; -3.61793, 15.6539, 0, 2.15621, 0.181289, 0; 10.6362, -4.74967, 0, 2.23532, 0.165991, 0; 16.3085, 0.7165, 0, 2.54259, 0.0259119, 0; 11.0383, 15.1907, 0, 1.89328, -0.351524, 0; 8.1178, 18.286, 0, 1.36944, -0.196337, 0; -11.5565, 14.0818, 0, 2.30507, -0.0900733, 0; -11.6951, 10.445, 0, 1.93486, 0.234751, 0; -14.3228, -6.091, 0, 2.37524, 0.0259066, 0; -7.48967, -7.82977, 0, 2.46697, 0.0735241, 0; 11.7451, -9.7246, 0, 2.62274, -0.0280196, 0; 4.51513, -2.13207, 0, 2.28557, -0.217603, 0; 15.3393, -9.28173, 0, 2.37145, -0.363377, 0; 17.6801, -5.25373, 0, 2.62967, 0.0178667, 0; 14.1327, 3.67947, 0, 2.26629, -0.404828, 0; 17.5228, 6.72553, 0, 3.12139, 0.0680969, 0; 5.15637, 16.1249, 0, 2.06071, 0.0327338, 0; -1.6151, 19.0238, 0, 1.66834, -0.0254332, 0; -9.15293, 5.35007, 0, 2.27604, 0.214091, 0; -14.119, 1.74813, 0, 2.35476, -0.0109954, 0; -14.8341, 8.31753, 0, 2.56064, 0.00874776, 0; -12.62, 5.0605, 0, 1.94196, 0.0965895, 0; 5.5029, -10.1943, 0, 2.5566, 0.271128, 0; 9.14153, -12.7573, 0, 1.94703, -0.0260447, 0; -2.8751, -3.56267, 0, 1.98007, 0.00966121, 0; -4.3444, 7.54107, 0, 2.06518, 0.1637, 0; 6.4481, 2.06747, 0, 1.99307, -0.0430748, 0; -0.610467, 1.8757, 0, 2.1699, -0.281808, 0; -8.1953, 8.54167, 0, 1.89161, -0.365641, 0; -6.54123, 2.77993, 0, 1.84984, -0.0938508, 0; 14.4091, 16.5753, 0, 0.461685, -0.642708, 0; 11.7713, 10.1874, 0, 1.72108, -0.584222, 0; 3.24, 3.65013, 0, 2.09081, 0.0650753, 0; 10.116, 4.4175, 0, 2.42626, -0.595117, 0; -5.2582, -13.123, 0, 1.8599, 0.0290441, 0; -1.6958, -10.4551, 0, 2.51603, -0.192596, 0; -11.1345, -4.94677, 0, 1.89687, 0.0569909, 0; -14.2574, -10.4768, 0, 2.4114, 0.0143635, 0; ];
quiver3(q_vector_2d1(:,1), q_vector_2d1(:,2), q_vector_2d1(:,3), q_vector_2d1(:,4), q_vector_2d1(:,5), q_vector_2d1(:,6));

hold on;
xlabel('x (m)');
ylabel('y (m)');
title('head contour and flow rate vector (Fig. 1)');
hold on;