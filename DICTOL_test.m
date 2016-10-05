clc;
clear all;
close all;
fprintf('================= SRC ====================\n');
acc_src = SRC_top;
fprintf('================= LCKSVD =================\n');
acc_lcksvd = LCKSVD_top;
fprintf('================= DLSI ===================\n');
acc_dlsi = DLSI_top;
fprintf('================= FDDL ===================\n');
acc_fddl = FDDL_top;
fprintf('================= COPAR ================\n');
acc_COPAR = COPAR_top;
fprintf('================= D2L2R2 =================\n');
acc_d2l2r2 = D2L2R2_top;
fprintf('================= LRSDL ==================\n');
acc_lrsdl = LRSDL_top;
fprintf('================= Summaray =================\n');
fprintf('+--------------------------+\n')
fprintf('|  Method    |   Accuray   |\n')
fprintf('+------------+-------------+\n')
fprintf('|  SRC       |   %2.2f%%    |\n', 100*acc_src);
fprintf('|  LCKSVD1   |   %2.2f%%    |\n', 100*acc_lcksvd(1));
fprintf('|  LCKSVD2   |   %2.2f%%    |\n', 100*acc_lcksvd(2));
fprintf('|  DLSI      |   %2.2f%%    |\n', 100*acc_dlsi);
fprintf('|  FDDL      |   %2.2f%%    |\n', 100*acc_fddl);
fprintf('|  COPAR   |   %2.2f%%    |\n', 100*acc_COPAR);
fprintf('|  D2L2R2    |   %2.2f%%    |\n', 100*acc_d2l2r2);
fprintf('|  LRSDL     |   %2.2f%%    |\n', 100*acc_lrsdl);
fprintf('+--------------------------+\n')

