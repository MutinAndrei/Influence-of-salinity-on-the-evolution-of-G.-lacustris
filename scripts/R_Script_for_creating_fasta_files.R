#Packages
library(xlsx)
library(reutils)
#Table from the supplement of Hou et al., 2022 (DOI: 10.1111/gc.160)
S1 <- read.xlsx("salt_table.xlsx", 2)
S1 <- S1[-c(88,89,29, 79, 62, 68, 82), ]
S1 <- S1[!duplicated(S1), ]



#Table transformations for the gene 28S
S28 <- S1
S28 <- S28[-which(S28$X28S == "-"), ]
S28 <- S28[, c(4,7,8,11,17)]
S28 <- S28[!is.na(S28$X28S), ]
S28 <- S28[!duplicated(S28), ]
S28[is.na(S28)] <- "-"
S28[S28$Voucher.number == "-" , 4] <- c(1:3)
#Receiving .fasta file by ID for the 28S gene
sum(S28$Voucher.number == "-")
S28_for_seq <- S28[ , 5]
S28fastas <- efetch(uid = S28_for_seq, db = "nucleotide", rettype = "fasta", retmode = "text")
write(content(S28fastas), "S28_lacustris_HOU_table_salt.fa")
#Preparing the vector for renaming .fasta file 28S
S28 <- S28[, -c(5)]
S28_un <- transform(S28, newcol=paste(Voucher.number, salinity..salty.brackish.fresh., GMYC.entity, Group, sep="_"))
S28_un <- S28_un[ , -c(1:4)]
#Table transformations for the gene COI
COI <- S1
COI <- COI[-which(COI$COI == "-"), ]
COI <- COI[, c(4,7,8,11,18)]
COI <- COI[!is.na(COI$COI), ]
COI <- COI[!duplicated(COI), ]
COI[is.na(COI)] <- "-"
COI[COI$Voucher.number == "-" , 4] <- c(4:23)
#Receiving .fasta file by ID for the COI gene
COI_for_seq <- COI[ , 5]
COIfastas <- efetch(uid = COI_for_seq, db = "nucleotide", rettype = "fasta", retmode = "text")
write(content(COIfastas), "COI_lacustris_HOU_table_salt.fa")
#Preparing the vector for renaming .fasta file COI
COI <- COI[, -c(5)]
COI_un <- transform(COI, newcol=paste(Voucher.number, salinity..salty.brackish.fresh., GMYC.entity, Group, sep="_"))
COI_un <- COI_un[ , -c(1:4)]
#Table transformations for the gene COIII
COIII <- S1
COIII <- COIII[-which(COIII$COIII == "-"), ]
COIII <- COIII[, c(4,7,8,11,19)]
COIII <- COIII[!is.na(COIII$COIII), ]
COIII <- COIII[!duplicated(COIII), ]
COIII[is.na(COIII)] <- "-"
#Receiving .fasta file by ID for the COIII gene
COIII_for_seq <- COIII[ , 5]
COIIIfastas <- efetch(uid = COIII_for_seq, db = "nucleotide", rettype = "fasta", retmode = "text")
write(content(COIIIfastas), "COIII_lacustris_HOU_table_salt.fa")
#Preparing the vector for renaming .fasta file COIII
COIII <- COIII[, -c(5)]
COIII_un <- transform(COIII, newcol=paste(Voucher.number, salinity..salty.brackish.fresh., GMYC.entity, Group, sep="_"))
COIII_un <- COIII_un[ , -c(1:4)]
#Table transformations for the gene S12-16
S12_16 <- S1
S12_16 <- S12_16[-which(S12_16$X12.16S == "-"), ]
S12_16 <- S12_16[, c(4,7,8,11,20)]
S12_16 <- S12_16[!is.na(S12_16$X12.16S), ]
S12_16 <- S12_16[!duplicated(S12_16), ]
S12_16[is.na(S12_16)] <- "-"
S12_16[S12_16$Voucher.number == "-" , 4] <- c(25:32)
#Receiving .fasta file by ID for the S12-16 gene
sum(S12_16$Voucher.number == "-")
S12_16_for_seq <- S12_16[ , 5]
S12_16fastas <- efetch(uid = S12_16_for_seq, db = "nucleotide", rettype = "fasta", retmode = "text")
write(content(S12_16fastas), "S12_16_lacustris_HOU_table_salt.fa")
#Preparing the vector for renaming .fasta file S12-16
S12_16 <- S12_16[, -c(5)]
S12_16_un <- transform(S12_16, newcol=paste(Voucher.number, salinity..salty.brackish.fresh., GMYC.entity, Group, sep="_"))
S12_16_un <- S12_16_un[ , -c(1:4)]
#Table transformations for the gene EF1a
EF1a <- S1
EF1a <- EF1a[-which(EF1a$EF1a == "-"), ]
EF1a <- EF1a[, c(4,7,8,11,21)]
EF1a <- EF1a[!is.na(EF1a$EF1a), ]
EF1a <- EF1a[!duplicated(EF1a), ]
EF1a[is.na(EF1a)] <- "-"
#Receiving .fasta file by ID for the EF1a gene
sum(EF1a$Voucher.number == "-")
EF1a_for_seq <- EF1a[ , 5]
EF1afastas <- efetch(uid = EF1a_for_seq, db = "nucleotide", rettype = "fasta", retmode = "text")
write(content(EF1afastas), "EF1a_lacustris_HOU_table_salt.fa")
#Preparing the vector for renaming .fasta file EF1a
EF1a <- EF1a[, -c(5)]
EF1a_un <- transform(EF1a, newcol=paste(Voucher.number, salinity..salty.brackish.fresh., GMYC.entity, Group, sep="_"))
EF1a_un <- EF1a_un[ , -c(1:4)]
#Table transformations for the gene AK
AK <- S1
AK <- AK[-which(AK$AK == "-"), ]
AK <- AK[, c(4,7,8,11,22)]
AK <- AK[!is.na(AK$AK), ]
AK <- AK[!duplicated(AK), ]
AK[is.na(AK)] <- "-"
#Receiving .fasta file by ID for the AK gene
sum(AK$Voucher.number == "-")
AK_for_seq <- AK[ , 5]
AKfastas <- efetch(uid = AK_for_seq, db = "nucleotide", rettype = "fasta", retmode = "text")
write(content(AKfastas), "AK_lacustris_HOU_table_salt.fa")
#Preparing the vector for renaming .fasta file AK
AK <- AK[, -c(5)]
AK_un <- transform(AK, newcol=paste(Voucher.number, salinity..salty.brackish.fresh., GMYC.entity, Group, sep="_"))
AK_un <- AK_un[ , -c(1:4)]
#Table transformations for the gene PEPCK
PEPCK <- S1
PEPCK <- PEPCK[-which(PEPCK$PEPCK == "-"), ]
PEPCK <- PEPCK[, c(4,7,8,11,23)]
PEPCK <- PEPCK[!is.na(PEPCK$PEPCK), ]
PEPCK <- PEPCK[!duplicated(PEPCK), ]
PEPCK[is.na(PEPCK)] <- "-"
#Receiving .fasta file by ID for the PEPCK gene
sum(PEPCK$Voucher.number == "-")
PEPCK_for_seq <- PEPCK[ , 5]
PEPCKfastas <- efetch(uid = PEPCK_for_seq, db = "nucleotide", rettype = "fasta", retmode = "text")
write(content(PEPCKfastas), "PEPCK_lacustris_HOU_table_salt.fa")
#Preparing the vector for renaming .fasta file PEPCK
PEPCK <- PEPCK[, -c(5)]
PEPCK_un <- transform(PEPCK, newcol=paste(Voucher.number, salinity..salty.brackish.fresh., GMYC.entity, Group, sep="_"))
PEPCK_un <- PEPCK_un[ , -c(1:4)]
#Table transformations for the gene HSP70
HSP70 <- S1
HSP70 <- HSP70[-which(HSP70$HSP70 == "-"), ]
HSP70 <- HSP70[, c(4,7,8,11,24)]
HSP70 <- HSP70[!is.na(HSP70$HSP70), ]
HSP70 <- HSP70[!duplicated(HSP70), ]
HSP70[is.na(HSP70)] <- "-"
#Receiving .fasta file by ID for the HSP70 gene
sum(HSP70$Voucher.number == "-")
HSP70_for_seq <- HSP70[ , 5]
HSP70fastas <- efetch(uid = HSP70_for_seq, db = "nucleotide", rettype = "fasta", retmode = "text")
write(content(HSP70fastas), "HSP70_lacustris_HOU_table_salt.fa")
#Preparing the vector for renaming .fasta file HSP70
HSP70 <- HSP70[, -c(5)]
HSP70_un <- transform(HSP70, newcol=paste(Voucher.number, salinity..salty.brackish.fresh., GMYC.entity, Group, sep="_"))
HSP70_un <- HSP70_un[ , -c(1:4)]


library(Biostrings)
#renaming .file fasta
S28_seq <- readDNAStringSet("S28_lacustris_HOU_table_salt.fa")
names(S28_seq) <- S28_un
writeXStringSet(S28_seq, "S28_lacustris_HOU_table_salt.fa")
#renaming .file fasta
COI_seq <- readDNAStringSet("COI_lacustris_HOU_table_salt.fa")
names(COI_seq) <- COI_un
writeXStringSet(COI_seq, "COI_lacustris_HOU_table_salt.fa")
#renaming .file fasta
COIII_seq <- readDNAStringSet("COIII_lacustris_HOU_table_salt.fa")
names(COIII_seq) <- COIII_un
writeXStringSet(COIII_seq, "COIII_lacustris_HOU_table_salt.fa")
#renaming .file fasta
S12_16_seq <- readDNAStringSet("S12_16_lacustris_HOU_table_salt.fa")
names(S12_16_seq) <- S12_16_un
writeXStringSet(S12_16_seq, "S12_16_lacustris_HOU_table_salt.fa")
#renaming .file fasta
EF1a_seq <- readDNAStringSet("EF1a_lacustris_HOU_table_salt.fa")
names(EF1a_seq) <- EF1a_un
writeXStringSet(EF1a_seq, "EF1a_lacustris_HOU_table_salt.fa")
#renaming .file fasta
AK_seq <- readDNAStringSet("AK_lacustris_HOU_table_salt.fa")
names(AK_seq) <- AK_un
writeXStringSet(AK_seq, "AK_lacustris_HOU_table_salt.fa")
#renaming .file fasta
PEPCK_seq <- readDNAStringSet("PEPCK_lacustris_HOU_table_salt.fa")
names(PEPCK_seq) <- PEPCK_un
writeXStringSet(PEPCK_seq, "PEPCK_lacustris_HOU_table_salt.fa")
#renaming .file fasta
HSP70_seq <- readDNAStringSet("HSP70_lacustris_HOU_table_salt.fa")
names(HSP70_seq) <- HSP70_un
writeXStringSet(HSP70_seq, "HSP70_lacustris_HOU_table_salt.fa")
