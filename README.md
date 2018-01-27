# GwasTriosCLogit
Reads PLINK GWAS trios data into R and tests each SNP independently with interaction terms through an additive genotype model using conditional logistic regression (parallel processing enabled). Produces summary stats for each genotype, ANOVA and clogit P values, and odds ratios.

NB - the original trios data for this project was converted to pseudocontrols using pseudocons: https://www.staff.ncl.ac.uk/richard.howey/pseudocons/index.html

<h1>Assumptions:</h1>
<ul>
  <li>PLINK data is uncompressed (<i>MyTrios.Chr[X].ped</i>, <i>MyTrios.Chr[X].map</i>, <i>MyTrios.Chr[X].fam</i>)</li>
  <li>PLINK data is divided up by chromosome, with sex chromosomes X and Y encoded as 23 and 24</li>
  <li>PLINK data is trios, with at least 1 case and control per family</li>
  <li>Genotype alleles are encoded as 1|2|3|4, i.e., using PLINK's <i>--allele1234</i> parameter</li>
  <li>PCA has been performed on each chromosome to produce at least 2 eigenvectors, saved as <i>MyTrios.Chr[X].eigenvec</i></li>
  <li>Minor (A1) and major (A2) alleles across the dataset are already identified and stored in columns #2 and #3 of <i>MajorMinor.tsv</i>, respectively, with column #1 being the genotype identifier</li>
  <li>Extra phenotype data is stored in <i>Metadata.csv</i>, which has as the first 2 columns the IID and FID. In our example, extra phenotype data includes age, MaternalSmokingPregnancy, and EnvironmentalTobaccoSmoke</li>
  </ul>

<h1>Extra notes:</h1>
<ul>
  <li>Number of CPU cores to use is set at the beginning of the script with <i>cpucores <- 16</i></li>
  <li>Adjustments in the model are made using age, gender, PC1, and PC2</li>
</ul>

<h1>Required packages:</h1>
<ul>
  <li>parallel</li>
  <li>doParallel</li>
  <li>survival</li>
  <li>MASS</li>
</ul>
  
<h1>Summary stats produced:</h1>
<ul>
  <li>SNP</li>
<li>MinorAllele</li>
<li>MajorAllele</li>
<li>Missingness</li>
<li>MissingnessCases</li>
<li>MissingnessControls</li>
<li>NumCasesWithGenotype</li>
<li>NumControlsWithGenotype</li>
<li>TallyMinorAlleleCases</li>
<li>TallyMinorAlleleControls</li>
<li>FreqMinorAlleleCases</li>
  <li>FreqMinorAlleleControls</li>
</ul>
