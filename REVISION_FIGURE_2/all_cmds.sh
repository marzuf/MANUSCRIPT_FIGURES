#!/bin/bash

#./all_cmds.sh



Rscript revision_fig5_pvalPlot_v2.R

exit 0


Rscript convert_to_matrix.R RWPE1 20 chr12
Rscript convert_to_matrix.R RWPE1 20 chr7
Rscript convert_to_matrix.R RWPE1 20 chr17
Rscript convert_to_matrix.R RWPE1 10 chr12
Rscript convert_to_matrix.R RWPE1 10 chr7
Rscript convert_to_matrix.R RWPE1 10 chr17

Rscript convert_to_matrix.R 22Rv1 20 chr12
Rscript convert_to_matrix.R 22Rv1 20 chr7
Rscript convert_to_matrix.R 22Rv1 20 chr17
Rscript convert_to_matrix.R 22Rv1 10 chr12
Rscript convert_to_matrix.R 22Rv1 10 chr7
Rscript convert_to_matrix.R 22Rv1 10 chr17
exit 0


java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar dump observed KR /mnt/ndata/Yuanlong/2.Results/1.Juicer/mega_GSE118514_RWPE1/mega/aligned/inter.hic 12 12 BP 10000 RWPE1_chr12_obs_KR_10kb.txt
java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar dump observed KR /mnt/ndata/Yuanlong/2.Results/1.Juicer/mega_GSE118514_RWPE1/mega/aligned/inter.hic 12 12 BP 20000 RWPE1_chr12_obs_KR_20kb.txt
java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar dump observed KR /mnt/ndata/Yuanlong/2.Results/1.Juicer/mega_GSE118514_RWPE1/mega/aligned/inter.hic 12 12 BP 100000 RWPE1_chr12_obs_KR_100kb.txt

java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar dump observed KR /mnt/ndata/Yuanlong/2.Results/1.Juicer/mega_GSE118514_RWPE1/mega/aligned/inter.hic 17 17 BP 10000 RWPE1_chr17_obs_KR_10kb.txt
java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar dump observed KR /mnt/ndata/Yuanlong/2.Results/1.Juicer/mega_GSE118514_RWPE1/mega/aligned/inter.hic 17 17 BP 20000 RWPE1_chr17_obs_KR_20kb.txt
java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar dump observed KR /mnt/ndata/Yuanlong/2.Results/1.Juicer/mega_GSE118514_RWPE1/mega/aligned/inter.hic 17 17 BP 100000 RWPE1_chr17_obs_KR_100kb.txt

java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar dump observed KR /mnt/ndata/Yuanlong/2.Results/1.Juicer/mega_GSE118514_RWPE1/mega/aligned/inter.hic 7 7 BP 10000 RWPE1_chr7_obs_KR_10kb.txt
java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar dump observed KR /mnt/ndata/Yuanlong/2.Results/1.Juicer/mega_GSE118514_RWPE1/mega/aligned/inter.hic 7 7 BP 20000 RWPE1_chr7_obs_KR_20kb.txt
java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar dump observed KR /mnt/ndata/Yuanlong/2.Results/1.Juicer/mega_GSE118514_RWPE1/mega/aligned/inter.hic 7 7 BP 100000 RWPE1_chr17_obs_KR_100kb.txt

java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar dump observed KR /mnt/ndata/Yuanlong/2.Results/1.Juicer/mega_GSE118514_RWPE1/mega/aligned/inter.hic 12 12 BP 10000 RWPE1_chr12_obs_KR_10kb.txt
java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar dump observed KR /mnt/ndata/Yuanlong/2.Results/1.Juicer/mega_GSE118514_RWPE1/mega/aligned/inter.hic 12 12 BP 20000 RWPE1_chr12_obs_KR_20kb.txt
java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar dump observed KR /mnt/ndata/Yuanlong/2.Results/1.Juicer/mega_GSE118514_RWPE1/mega/aligned/inter.hic 12 12 BP 100000 RWPE1_chr12_obs_KR_100kb.txt


java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar dump observed KR /mnt/ptemp/Yuanlong/2.Results/1.Juicer/GSE118514_22Rv1/rep1/run1/aligned/inter.hic 17 17 BP 10000 22Rv1_chr17_obs_KR_10kb.txt
java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar dump observed KR /mnt/ptemp/Yuanlong/2.Results/1.Juicer/GSE118514_22Rv1/rep1/run1/aligned/inter.hic 17 17 BP 20000 22Rv1_chr17_obs_KR_20kb.txt
java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar dump observed KR /mnt/ptemp/Yuanlong/2.Results/1.Juicer/GSE118514_22Rv1/rep1/run1/aligned/inter.hic 17 17 BP 100000 22Rv1_chr17_obs_KR_100kb.txt

java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar dump observed KR /mnt/ptemp/Yuanlong/2.Results/1.Juicer/GSE118514_22Rv1/rep1/run1/aligned/inter.hic 7 7 BP 10000 22Rv1_chr7_obs_KR_10kb.txt
java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar dump observed KR /mnt/ptemp/Yuanlong/2.Results/1.Juicer/GSE118514_22Rv1/rep1/run1/aligned/inter.hic 7 7 BP 20000 22Rv1_chr7_obs_KR_20kb.txt
java -Xms512m -Xmx2048m -jar /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar dump observed KR /mnt/ptemp/Yuanlong/2.Results/1.Juicer/GSE118514_22Rv1/rep1/run1/aligned/inter.hic 7 7 BP 100000 22Rv1_chr7_obs_KR_100kb.txt
