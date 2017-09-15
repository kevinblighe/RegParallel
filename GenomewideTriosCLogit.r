
R

cpucores <- 16

#MICROSOFT R OPEN (MRO) ONLY
###########################
#MKL libraries are installed alongside the MRO instalation
#Check the number of CPU cores being used
getMKLthreads()

#Set the # of CPU cores to use
setMKLthreads(cpucores)
###########################

#Install the parallel processing packges required for certain functions not covered by MRO's parallel processing capabilites
require(parallel)

#Set number of CPU cores 
options("mc.cores"=cpucores)

require(doParallel)
cores <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cpucores)
Sys.setenv("MC_CORES"=cpucores)



#########################################
#Read the PLINK data into R
#########################################

#Create output directories
dir.create(file.path("/MyDir", "EnvironmentalTobaccoSmoke"), showWarnings=FALSE)
dir.create(file.path("/MyDir", "MaternalSmokingPregnancy"), showWarnings=FALSE)
dir.create(file.path("/MyDir", "AlleleInfo"), showWarnings=FALSE)

#Save the output directories and interaction terms as vectors
directories <- c("/MyDir/EnvironmentalTobaccoSmoke","/MyDir/MaternalSmokingPregnancy")
interactions <- c("EnvironmentalTobaccoSmoke","MaternalSmokingPregnancy")

#Loop through each chromosome
for (i in c(1:24)))
{
	print(paste("Starting analysis for chromosome ", i, sep=""))

	#Read in the plain text pseudocontrols data
	ped <- read.table(paste("MyTrios.Chr", i, ".ped", sep=""), header=FALSE, stringsAsFactors=FALSE)
	ped <- ped[,7:ncol(ped)]
	map <- read.table(paste("MyTrios.Chr", i, ".map", sep=""), header=FALSE, stringsAsFactors=FALSE)
	fam <- read.table(paste("MyTrios.Chr", i, ".fam", sep=""), header=FALSE, stringsAsFactors=FALSE)
	colnames(fam) <- c("FID", "IID", "PID", "MID", "gender", "CaseControl")

	#Merge each genotype together (genotypes are pairwise columns)
	for (col in seq(1, ncol(ped), 2))
	{
		if (col==1)
		{
			ped.new <- data.frame(factor(paste(ped[,col], ped[,col+1], sep="")))
		}
		else
		{
			ped.new <- cbind(ped.new, factor(paste(ped[,col], ped[,col+1], sep="")))
		}
	}
	rownames(ped.new) <- fam$IID
	colnames(ped.new) <- map[,2]



	#########################################
	#Read in phenotype data for the trios
	#########################################

	trio.metadata <- read.csv("Metadata.csv", stringsAsFactors=FALSE)

	#Extract phenotype data of interest (listed above)
	#NB - first two columns, FID and IID also extracted
	phenotypes <- data.frame(trio.metadata[,1:2], trio.metadata[,c("age","MaternalSmokingPregnancy","EnvironmentalTobaccoSmoke")])

	#Read in eigenvectors from PCA
	eigenvectors <- read.csv(paste("MyTrios.Chr", i, ".eigenvec", sep=""), header=FALSE, sep=" ")

	#Create a new data-frame that will form that used in conditional logistic regression
	conditional <- data.frame(fam, eigenvectors[,3:4], ped.new)
	colnames(conditional)[c(7,8)] <- c("PC1", "PC2")

	#Only keep in the phenotypes file any samples that also have genotype data
	phenotypes <- phenotypes[which(phenotypes[,2] %in% rownames(conditional)),]

	#Merge the covariates data to the genotype data
	conditional <- data.frame(phenotypes[match(rownames(conditional), phenotypes[,2]), 3:ncol(phenotypes)], conditional, row.names=rownames(conditional))

	#Check all variables
	conditional$age	#continuous
	conditional$gender #Categorical: M, male; F, female
	conditional$PC1	#Continuous
	conditional$PC2	#Continuous
	conditional$FID #Categorical
	conditional$CaseControl #Numeric-categorical: 1, case; 0, pseudo-control
	conditional$EnvironmentalTobaccoSmoke #Categorical: 1, exposed; 0, not exposed
	conditional$MaternalSmokingPregnancy #Categorical: 1, smoked; 0, not smoked

	#Convert 00 genotypes (missing) to NA
	conditional[conditional=="00"] <- NA



	#################################
	#Convert genotypes into numerical for additive modeling
	#################################
	#Read in the numerical encoding for genotypes (A1 is minor allele)
	AlleleInfo <- read.table("MajorMinor.tsv", header=TRUE, sep="\t")

	print(paste("Converting genotypes to 0 (A2|A2), 1 (A1|A2), 2 (A1|A1) for additive disease model...", sep=""))

	#For writing
	dfAlleleInfo <- data.frame(row.names=colnames(conditional[,12:ncol(conditional)]))

	#For double-checking that the true minor allele is A1 in the original data (saved in MajorMinor.tsv)
	write.table("SNP\tGenotype\tTally", paste("/MyDir/AlleleInfo/chr", i, ".AlleleCheck.tsv", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE)

	#Count number of cases and controls
	nCases <- nrow(conditional[conditional$CaseControl==1,])
	nControls <- nrow(conditional[conditional$CaseControl==0,])

	#Lookup through each SNP and convert
	for (m in 12:ncol(conditional))
	{
		genotypes <- sort(table(conditional[,m]), decreasing=FALSE)
		genotypes <- genotypes[names(genotypes)!="00"]
		genotypes <- data.frame(colnames(conditional)[m], genotypes)
		write.table(genotypes, paste("/MyDir/AlleleInfo/chr", i, ".AlleleCheck.tsv", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)

		#Determine minor (A1) and major (A2) alleles per genotype
		minor.allele <- AlleleInfo[which(AlleleInfo[,1] %in% colnames(conditional)[m]),2]
		major.allele <- AlleleInfo[which(AlleleInfo[,1] %in% colnames(conditional)[m]),3]

		#Configure the genotypes as:
		#	Major allele, homozygote = A2|A2
		#	Het, A1|A2
		#	Het, A2|A1
		#	Minor allele, homozygote = A1|A1
		major.hom <- paste(major.allele, major.allele, sep="")
		het1 <- paste(minor.allele, major.allele, sep="")
		het2 <- paste(major.allele, minor.allele, sep="")
		minor.hom <- paste(minor.allele, minor.allele, sep="")

		#Substitute the genotypes in the dataframe
		conditional[,m] <- gsub(major.hom, 0, conditional[,m])
		conditional[,m] <- gsub(het1, 1, conditional[,m])
		conditional[,m] <- gsub(het2, 1, conditional[,m])
		conditional[,m] <- gsub(minor.hom, 2, conditional[,m])

		#Convert to numeric
		conditional[,m] <- as.numeric(conditional[,m])

		#Determine total missingness
		missingness <- table(is.na(conditional[,m]))[2]/(nCases + nControls)
		if (is.na(missingness)) {missingness <- 0}

		#Determine missingness in cases
		missingnessCases <- table(is.na(conditional[conditional$CaseControl==1,m]))[2]/(nCases)
		if (is.na(missingnessCases)) {missingnessCases <- 0}

		#Determine missingness in controls
		missingnessControls <- table(is.na(conditional[conditional$CaseControl==0,m]))[2]/(nControls)
		if (is.na(missingnessControls)) {missingnessControls <- 0}

		#How many cases have genotype?
		nCasesWithoutMissingness <- table(is.na(conditional[conditional$CaseControl==1,m]))[1]
		nControlsWithoutMissingness <- table(is.na(conditional[conditional$CaseControl==0,m]))[1]

		#How many controls have genotype?
		nMinorAlleleCases <- sum(conditional[conditional$CaseControl==1,m], na.rm=TRUE)
		nMinorAlleleControls <- sum(conditional[conditional$CaseControl==0,m], na.rm=TRUE)

		#Frequency of minor allele in cases with genotype
		freqMinorAlleleCases <- nMinorAlleleCases / ((nCases-(nCases*missingnessCases))*2)
		if (is.na(freqMinorAlleleCases)) {freqMinorAlleleCases <- nMinorAlleleCases / (nCases*2)}

		#Frequency of minor allele in controls with genotype
		freqMinorAlleleControls <- nMinorAlleleControls / ((nControls-(nControls*missingnessControls))*2)
		if (is.na(freqMinorAlleleControls)) {freqMinorAlleleControls <- nMinorAlleleControls / (nControls*2)}

		#Merge results
		dfAlleleInfo <- rbind(dfAlleleInfo, data.frame(colnames(conditional)[m], minor.allele, major.allele, missingness, missingnessCases, missingnessControls, nCasesWithoutMissingness, nControlsWithoutMissingness, nMinorAlleleCases, nMinorAlleleControls, freqMinorAlleleCases, freqMinorAlleleControls))
	}

	#Output the results fr each chromosome
	colnames(dfAlleleInfo) <- c("SNP", "MinorAllele", "MajorAllele", "Missingness", "MissingnessCases", "MissingnessControls", "NumCasesWithGenotype", "NumControlsWithGenotype", "TallyMinorAlleleCases", "TallyMinorAlleleControls", "FreqMinorAlleleCases", "FreqMinorAlleleControls")
	write.table(dfAlleleInfo, paste("/MyDir/AlleleInfo/chr", i, ".AlleleInfo.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE)



	#################################
	#Conditional logistic regression
	#################################
	require(survival)
	require(MASS)

	print(paste("Starting conditional logistic regression for ", ncol(conditional)-11, " genotypes...", sep=""))

	#Open up output files for writing
	for (k in 1:length(interactions))
	{
		print(paste("Interaction term to be studied will be: ", interactions[k], sep=""))
		print(paste("Output direction will be: ", directories[k], sep=""))

		#Output the file header
		write.table(c("chr:SNP:bp\tTerm\tp.value"), paste(directories[k], "/chr", i, ".LRT.tsv", sep=""), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=FALSE)
		write.table(c("chr:SNP:bp\tTerm\tCoef\texp(Coef)\tCoef.SE\tZ.score\tp.value"), paste(directories[k], "/chr", i, ".tsv", sep=""), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=FALSE)

		#Create en empty list to hold the formulae
		formula.list <- list()

		#Store each possible formula in the list
		for (j in 12:ncol(conditional))
		{
			formula.list[[j-11]] <- as.formula(paste("CaseControl ~ strata(FID) + age + gender + PC1 + PC2 + ", interactions[k], "*", colnames(conditional)[j], sep=""))
		}

		#Run a base model without the SNP - this will replace later models that fail
		formula <- as.formula(paste("CaseControl ~ strata(FID) + age + gender + PC1 + PC2 + ", interactions[k], sep=""))
		basemodel <- clogit(formula, conditional, method="breslow", ties="breslow", singular.ok=TRUE)

		#varnames <- colnames(conditional[,12:ncol(conditional)])
		varnames <- paste(paste("chr", map[,1], sep=""), map[,2], map[,4], sep=":")
		blocks <- floor((length(varnames))/5000)+1

		for (l in 1:blocks)
		{
			#Run the conditional logistic regression using mclappy, i.e., multiple cores
			#First block
			if (l==1)
			{
				print(paste("Processing 5,000 variants, batch ", l, " of ", blocks, sep=""))
				print(paste("Index1: 1; ", "Index2: ", (5000*l), sep=""))
				models <- mclapply(formula.list[1:(5000*l)], function(x) clogit(x, conditional, method="breslow", ties="breslow", singular.ok=TRUE))
				names(models) <- varnames[1:(5000*l)]
			}

			#Final block
			if (l==blocks)
			{
				print(paste("Processing final batch ", l, " of ", blocks, sep=""))
				print(paste("Index1: ", (1+(5000*(l-1))), "; ", "Index2: ", length(formula.list), sep=""))
				models <- mclapply(formula.list[(1+(5000*(l-1))):length(formula.list)], function(x) clogit(x, conditional, method="breslow", ties="breslow", singular.ok=TRUE))
				names(models) <- varnames[(1+(5000*(l-1))):length(formula.list)]
			}

			#Any other blocks
			if (l>1 && l<blocks)
			{
				print(paste("Processing 5,000 variants, batch ", l, " of ", blocks, sep=""))
				print(paste("Index1: ", (1+(5000*(l-1))), "; ", "Index2: ", (5000*l), sep=""))
				models <- mclapply(formula.list[(1+(5000*(l-1))):(5000*l)], function(x) clogit(x, conditional, method="breslow", ties="breslow", singular.ok=TRUE))
				names(models) <- varnames[(1+(5000*(l-1))):(5000*l)]
			}

			#If any models failed, detect them and replace with the base model
			wObjects <- mclapply(models, function(x) if (is.null(names(x))) { x <- basemodel } else {x} )

			LRT <- mclapply(wObjects, function(x) 1-pchisq(summary(x)$logtest[1], summary(x)$logtest[2],lower.tail=TRUE))
			LRT <- do.call(rbind, lapply(LRT, data.frame, stringsAsFactors=FALSE))
			LRT <- data.frame(rownames(LRT), rep(paste("SNP * ", interactions[k], sep=""), nrow(LRT)), LRT[,1], row.names=rownames(LRT))

			#Exract information from the model and store in a new list
			wObjects <- mclapply(wObjects, function(x) data.frame(rownames(summary(x)$coefficients), summary(x)$coefficients))

			#Remove age, gender, PC1, and PC2 information
			wObjects <- mclapply(wObjects, function(x) x[grep("age|gender|PC1|PC2", rownames(x), invert=TRUE),])

			#Convert the list into a data frame for writing
			wObject <- do.call(rbind, lapply(wObjects, data.frame, stringsAsFactors=FALSE))
			wObject <- data.frame(gsub("\\.[a-zA-Z0-9_():]*", "", rownames(wObject)), wObject)
			wObject[,2] <- gsub("factor", "", wObject[,2])

			#Output the results to disk
			write.table(LRT, paste(directories[k], "/chr", i, ".LRT.tsv", sep=""), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)
			write.table(wObject, paste(directories[k], "/chr", i, ".tsv", sep=""), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)
		}
	}
}

q()

