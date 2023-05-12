#!/bin/bash
#SBATCH --job-name=CfelMapping
#SBATCH --output=/scratch/tkay/SNG/code/%x.o.log
#SBATCH --error=/scratch/tkay/SNG/code/%x.e.log
#SBATCH --mail-user=tkay@unil.ch
#SBATCH --mail-type=FAIL,END
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3G
#SBATCH --time=23:00:00


module add gcc/9.3.0;
module add star/2.7.8a;

# Generate genome indices
cd /scratch/tkay/SNG/data
mkdir ./genomes/STARgenome

STAR --runThreadN 20 --runMode genomeGenerate --genomeDir ./genomes/STARgenome --genomeFastaFiles ./genomes/Cfel.assembly2.0.fasta  --genomeSAindexNbases 12

# Mapping
mkdir ./aligned_transcriptomes

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A100_L1_R1_001_P0rheoFEWw0d.fastq.gz ./transcriptomes/A100_L1_R2_001_pApCLraXdVMu.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A100

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A101_L1_R1_001_bh0nLmwPFu09.fastq.gz ./transcriptomes/A101_L1_R2_001_uw55T9FK8nRK.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A101

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A102_L1_R1_001_qhmfhj4Vwep5.fastq.gz ./transcriptomes/A102_L1_R2_001_G6XmXoL8xSmM.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A102

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A103_L1_R1_001_4Lm1iGIAf4f9.fastq.gz ./transcriptomes/A103_L1_R2_001_wFDe0zIi3ElS.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A103

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A104_L1_R1_001_RTjsagx3TyH8.fastq.gz ./transcriptomes/A104_L1_R2_001_QalIUKiIb7Z2.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A104

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A105_L1_R1_001_62yJsCboWPZF.fastq.gz ./transcriptomes/A105_L1_R2_001_Aco5Mkerwa2z.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A105

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A106_L1_R1_001_htexHGv3q0wb.fastq.gz ./transcriptomes/A106_L1_R2_001_DxG3ANz53fQ4.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A106

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A107_L1_R1_001_AsKS9PlU2htX.fastq.gz ./transcriptomes/A107_L1_R2_001_zJmBNA4dpDoV.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A107

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A108_L1_R1_001_MxNNhbVgHkBE.fastq.gz ./transcriptomes/A108_L1_R2_001_RuFt8vgihJPP.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A108

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A109_L1_R1_001_Zw5CJLViyd2L.fastq.gz ./transcriptomes/A109_L1_R2_001_v9QPLQIcmjwH.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A109

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A110_L1_R1_001_eQpXm8bYP149.fastq.gz ./transcriptomes/A110_L1_R2_001_qkrBqIOgJU63.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A110

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A111_L1_R1_001_I81QXlF2J7QQ.fastq.gz ./transcriptomes/A111_L1_R2_001_tnZHL9xqfzfm.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A111

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A112_L1_R1_001_OuQ2MDeappMZ.fastq.gz ./transcriptomes/A112_L1_R2_001_EyrIpANeUKNx.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A112

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A113_L1_R1_001_CWeIVEqCdE7L.fastq.gz ./transcriptomes/A113_L1_R2_001_5iqGu1JErU6N.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A113

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A115_L1_R1_001_xyywSoB2N6ZG.fastq.gz ./transcriptomes/A115_L1_R2_001_3pmbnSloDAVz.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A115

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A116_L1_R1_001_OjcWRKm9iwq8.fastq.gz ./transcriptomes/A116_L1_R2_001_LvjUMv58AHVE.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A116

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A117_L1_R1_001_ZdpDGpyaHTLn.fastq.gz ./transcriptomes/A117_L1_R2_001_6kOitAsKe0yw.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A117

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A118_L1_R1_001_qsTdkafQjVCX.fastq.gz ./transcriptomes/A118_L1_R2_001_5ljQQDucDIsY.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A118

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A119_L1_R1_001_LIKm9fF4WqRq.fastq.gz ./transcriptomes/A119_L1_R2_001_3hES2KCWiAVx.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A119

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A120_L1_R1_001_a8rnk8lK8jiU.fastq.gz ./transcriptomes/A120_L1_R2_001_uBd2KZ4MNrdy.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A120

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A121_L1_R1_001_Qjbgcfjdr7OX.fastq.gz ./transcriptomes/A121_L1_R2_001_rIohjdTM0Gpt.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A121

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A122_L2_R1_001_nKgF9xdESp12.fastq.gz ./transcriptomes/A122_L2_R2_001_bz5gLtFlGOGP.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A122

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A124_L2_R1_001_olnqH4hEOjSI.fastq.gz ./transcriptomes/A124_L2_R2_001_xvJ4Te3xKnRx.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A124

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A125_L2_R1_001_lpU2KCijQb7W.fastq.gz ./transcriptomes/A125_L2_R2_001_G2moWKAvnlsO.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A125

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A126_L2_R1_001_JpIa1fmpmHCy.fastq.gz ./transcriptomes/A126_L2_R2_001_uRaaKQIcmbtS.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A126

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A127_L2_R1_001_N0VBqOuyL8D9.fastq.gz ./transcriptomes/A127_L2_R2_001_DXf4SfLaDiF1.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A127

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A128_L2_R1_001_paFoa4fHYwQR.fastq.gz ./transcriptomes/A128_L2_R2_001_097F6u5IPiIj.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A128

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A12_L1_R1_001_eoOxZQEyS98O.fastq.gz ./transcriptomes/A12_L1_R2_001_Ph4mIPh9QBY2.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A012

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A130_L2_R1_001_IH3Aej2gqOQJ.fastq.gz ./transcriptomes/A130_L2_R2_001_jCOFCYwcmK85.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A130

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A132_L2_R1_001_To1QrWULkWqs.fastq.gz ./transcriptomes/A132_L2_R2_001_AhqKaK47GtAU.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A132

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A133_L2_R1_001_RxsbwfAkYjgz.fastq.gz ./transcriptomes/A133_L2_R2_001_FwRxcxSyaLTH.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A133

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A134_L2_R1_001_NIz9bdQLgNGJ.fastq.gz ./transcriptomes/A134_L2_R2_001_4aG7QbFmz5TS.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A134

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A135_L2_R1_001_47FMKdt5iM5g.fastq.gz ./transcriptomes/A135_L2_R2_001_nDaEp7lv5PPv.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A135

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A136_L2_R1_001_GH6C9ptrAfa5.fastq.gz ./transcriptomes/A136_L2_R2_001_I1BhdgffEsuT.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A136

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A137_L2_R1_001_ofZoCXjyZnTv.fastq.gz ./transcriptomes/A137_L2_R2_001_7zDdVlhCOorl.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A137

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A138_L2_R1_001_KrC68jnkqZXG.fastq.gz ./transcriptomes/A138_L2_R2_001_m5dA6gkZwTYu.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A138

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A140_L2_R1_001_RFcgOAhlivST.fastq.gz ./transcriptomes/A140_L2_R2_001_IyfJS5zdI63c.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A140

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A14_L1_R1_001_vCSotH72J6b2.fastq.gz ./transcriptomes/A14_L1_R2_001_i9TRSQnzpWep.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A014

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A15_L1_R1_001_cP96lbOaav6i.fastq.gz ./transcriptomes/A15_L1_R2_001_2dDzI5GF9gTd.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A015

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A16_L1_R1_001_fT2Ui89BW7Ub.fastq.gz ./transcriptomes/A16_L1_R2_001_Sm5K0cefCYYy.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A016

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A17_L1_R1_001_1LTjtlR19C5J.fastq.gz ./transcriptomes/A17_L1_R2_001_6OHP6HqrE9cl.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A017

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A19_L1_R1_001_cBRXgEY5yiSk.fastq.gz ./transcriptomes/A19_L1_R2_001_w2IHs9CxhW1t.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A019

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A1_L2_R1_001_uZqzhdDvSJJ7.fastq.gz ./transcriptomes/A1_L2_R2_001_w4EanZDBJegq.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A001

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A20_L1_R1_001_i6RMkRbNHv1w.fastq.gz ./transcriptomes/A20_L1_R2_001_HGgzA3nm9gL5.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A020

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A21_L1_R1_001_Pk0znlWB1wvH.fastq.gz ./transcriptomes/A21_L1_R2_001_LK5LIMTt67OS.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A021

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A22_L1_R1_001_M3Gf4vo24VoU.fastq.gz ./transcriptomes/A22_L1_R2_001_QINsCqgBR3sP.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A022

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A23_L1_R1_001_WuNA2ByaMWgo.fastq.gz ./transcriptomes/A23_L1_R2_001_CH4Yu3H6TO0c.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A023

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A24_L1_R1_001_PqpgWSLdYzYl.fastq.gz ./transcriptomes/A24_L1_R2_001_QVDCaXBR8XXt.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A024

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A25_L1_R1_001_yDI2mGueYQwP.fastq.gz ./transcriptomes/A25_L1_R2_001_ZwUmH8f26CCa.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A025

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A26_L1_R1_001_zKBkYxigG6iC.fastq.gz ./transcriptomes/A26_L1_R2_001_Pky5qTWwlHpR.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A026

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A27_L1_R1_001_3vWF3nU9lokb.fastq.gz ./transcriptomes/A27_L1_R2_001_4BUla74RW13j.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A027

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A2_L2_R1_001_BKPrUf5BLsOU.fastq.gz ./transcriptomes/A2_L2_R2_001_g6R9TJ8rNltI.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A002

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A4_L2_R1_001_AUuR7oAb3XFa.fastq.gz ./transcriptomes/A4_L2_R2_001_Lq0a91QNBSlw.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A004

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A59_L1_R1_001_IivEZSYtLFRm.fastq.gz ./transcriptomes/A59_L1_R2_001_MBC4jh89jasK.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A059

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A5_L2_R1_001_jSfQ6NofHlxE.fastq.gz ./transcriptomes/A5_L2_R2_001_Qd8AtUcRYCdu.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A005

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A60_L1_R1_001_88uy0fgx9inV.fastq.gz ./transcriptomes/A60_L1_R2_001_vjkAeCGaJKnI.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A060

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A61_L1_R1_001_Bg5w3OfK4cf0.fastq.gz ./transcriptomes/A61_L1_R2_001_pdLHyFoqqmER.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A061

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A62_L1_R1_001_zOLmzbsIDfw8.fastq.gz ./transcriptomes/A62_L1_R2_001_0V8yzWQI2pdP.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A062

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A63_L1_R1_001_dMula3i8IMp5.fastq.gz ./transcriptomes/A63_L1_R2_001_IJYkLEkn9Y1u.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A063

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A64_L1_R1_001_i3riqHLRoW3M.fastq.gz ./transcriptomes/A64_L1_R2_001_ioIBYHHaMA15.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A064

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A65_L1_R1_001_fHs4J2NWTaRG.fastq.gz ./transcriptomes/A65_L1_R2_001_hCcEiugWIGQf.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A065

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A66_L1_R1_001_JU1g2acIYbjA.fastq.gz ./transcriptomes/A66_L1_R2_001_rYR8kcHN10kf.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A066

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A68_L1_R1_001_fn2bXKsW8vVW.fastq.gz ./transcriptomes/A68_L1_R2_001_Nj3AKEpdNVVt.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A068

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A69_L1_R1_001_7TtjhoRAmcjk.fastq.gz ./transcriptomes/A69_L1_R2_001_vjcoJelSrjYc.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A069

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A6_L2_R1_001_eneTuJrnp43A.fastq.gz ./transcriptomes/A6_L2_R2_001_1JY20OQM9Eeh.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A006

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A70_L1_R1_001_FPfZjTXTHVDK.fastq.gz ./transcriptomes/A70_L1_R2_001_VwbLhZ69eyQq.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A070

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A71_L1_R1_001_vPfp1pNvHZZZ.fastq.gz ./transcriptomes/A71_L1_R2_001_KxM2m0VX2f5r.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A071

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A73_L1_R1_001_ffi7UG3A19TP.fastq.gz ./transcriptomes/A73_L1_R2_001_O143V45ym7cE.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A073

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A74_L1_R1_001_oTBn00PxJ6OE.fastq.gz ./transcriptomes/A74_L1_R2_001_lcplhJlnfEMh.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A074

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A75_L1_R1_001_TCz1eKPTtD6f.fastq.gz ./transcriptomes/A75_L1_R2_001_2oNuSA445b9S.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A075

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A76_L1_R1_001_2AvhwL7wCthD.fastq.gz ./transcriptomes/A76_L1_R2_001_FhoFYSXiUiYp.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A076

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A77_L1_R1_001_rWX4IxoWdpCW.fastq.gz ./transcriptomes/A77_L1_R2_001_gnt7boCJIvbZ.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A077

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A79_L1_R1_001_4UbGAbyIabal.fastq.gz ./transcriptomes/A79_L1_R2_001_NG8NmmJkliWx.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A079

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A81_L1_R1_001_y9typkR7Tv8P.fastq.gz ./transcriptomes/A81_L1_R2_001_4OOBQ41Ogj5e.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A081

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A82_L1_R1_001_RCw6jFM8jtGJ.fastq.gz ./transcriptomes/A82_L1_R2_001_uQmG9vV89SNt.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A082

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A83_L1_R1_001_oae79L6ur0Ve.fastq.gz ./transcriptomes/A83_L1_R2_001_toGfjpQTmkt7.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A083

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A85_L1_R1_001_JYODEFce46QR.fastq.gz ./transcriptomes/A85_L1_R2_001_Zf6jASbWMv2Z.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A085

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A86_L1_R1_001_kFwIQTKcEBfw.fastq.gz ./transcriptomes/A86_L1_R2_001_jhOozHgl5HUC.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A086

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A8_L2_R1_001_DxOReMKCXFqE.fastq.gz ./transcriptomes/A8_L2_R2_001_8F5movD6cW4z.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A008

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A90_L1_R1_001_p5l52AZmqd1a.fastq.gz ./transcriptomes/A90_L1_R2_001_IvCzKpPdtAYn.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A090

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A91_L1_R1_001_Kvz4lUnxCEVI.fastq.gz ./transcriptomes/A91_L1_R2_001_ol3zQImCw6Rh.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A091

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A92_L1_R1_001_1fv2cRFYYaS0.fastq.gz ./transcriptomes/A92_L1_R2_001_t0tldanpxm23.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A092

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A93_L1_R1_001_fzMmkB9cSfR5.fastq.gz ./transcriptomes/A93_L1_R2_001_W9DwzRKOtQHB.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A093

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A94_L1_R1_001_drPcVeatdH5O.fastq.gz ./transcriptomes/A94_L1_R2_001_L30C7AYKJT68.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A094

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A95_L1_R1_001_zBwjLd8qsOTq.fastq.gz ./transcriptomes/A95_L1_R2_001_NvOi3aqX68Xh.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A095

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A96_L1_R1_001_UmE4asElfOxd.fastq.gz./transcriptomes/A96_L1_R2_001_olg2qO0V25Rz.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A096

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A97_L1_R1_001_OKyfjTmujuX2.fastq.gz ./transcriptomes/A97_L1_R2_001_yUC61Nrq2FP7.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A097

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A98_L1_R1_001_Xn6djU9o0mkG.fastq.gz ./transcriptomes/A98_L1_R2_001_FSnePwIZiE4S.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A098

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A99_L1_R1_001_eI8PjVEqM4Ub.fastq.gz ./transcriptomes/A99_L1_R2_001_8wDThfeLmYYx.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A099

STAR --runThreadN 20 --genomeDir ./STARgenome --readFilesIn ./transcriptomes/A9_L2_R1_001_WJZBKusFvsk6.fastq.gz ./transcriptomes/A9_L2_R2_001_wkH3x7WincWS.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ./aligned_transcriptomes/SNG_A009
