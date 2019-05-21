# Line1 CRC V2

## Code summary

### 1. Quality Control
* Quality control is primarily done through package PTXQC

### 2. Data preparation:
* Remove contaminants and reverse proteins
* Log transformation
* Separate different experiment
* Separate LFQ and iBAQ intensity

### 3. Imputating zero intensities - LFQ
#### Method_6 (imputing when all replicates in either case or control is zero):
# Remove entries that all case and control replicates are zeros.
* For each replica mean and standard deviation of non-zero proteins intensities are counted. 
* New intensity for each missing value in each replica is sampled from uniform distribution with parameters:
*  start = mu - 3*sd, end = mu - 2*sd
* Intnew=Unif(mean(Intreplica)-3*sd(Intreplica),mean(Intreplica)-2*sd(Intreplica))

####Method_6 (imputing partial NA values):
* For outliers (proteins, which have zero values in all but one replica) this method implies one of the methods for imputation of all zero replicas.
* For other proteins use following method:
  * Impute values for proteins, which have zero intensities only in some replicates:
  *  Build distribution of deltas for all non zero proteins, where 
    * delta =(Intrep1-Intrep2) mean(Intrep1,Intrep2)
    * Calculate mudelta , sddelta
    * Calculate new delta and new Intensity:
    * deltanew=rnorm(mu=mudelta, sd=sddelta*sqrt(2)*mean(correlations))
    * Inew=mean(Intother)*abs(1+deltanew)

###	4. Variance analysis test on imputed LFQ intensity
* Perform t-test between cases and controls
* Adjust p-values with Benjamin-Hochberg correction test
* Calculate log2fold change
* Select significant proteins with p.adj <0.05 & logfol d> 1

###	5. Normalize iBAQ intensity of significant proteins by LORF1 and integrate data