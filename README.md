# GenomewideCLogit
Reads in genomewide PLINK trios data to R and tests each SNP independently with interaction terms through conditional logistic regression (parallel processing enabled)

Assumptions
<ul>
  <li>PLINK data is uncompressed (MyTrios.Chr[X].ped, MyTrios.Chr[X].map, MyTrios.Chr[X].fam)</li>
  <li>PLINK data is divided up by chromosome, with sex chromosomes X and Y encoded as 23 and 24</li>
  <li>PLINK data is trios, with at least 1 case or control per family</li>
  <li>Genotype alleles are encoded as 1|2|3|4, i.e., using PLINK's <i>--allele1234</i> parameter</li>
  <li>PCA has been performed on each chromosome to produce at least 2 eigenvectors<, saved as <i>MyTrios.Chr[X].eigenvec</i></li>
  <li>Major and Minor alleles across the dataset are already identified and stored in columns #2 and #3 of <i>MajorMinor.tsv</i>, respectively, with column #1 being the genotype identifier</li>
  <li>Extra phenotype data is stored in <i>Metadata.csv</i>, which has as the first 2 columns the IID and FID. In our example, extra phenotype data includes age,MaternalSmokingPregnancy, and EnvironmentalTobaccoSmoke</li>
  </ul>