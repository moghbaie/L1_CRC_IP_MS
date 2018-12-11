# Line1 CRC V2

## Code summary

### 1. Quality Control
		Quality control is primarily done through package PTXQC

### 2. Data preparation:
  	 * Calculate iBAQ intensity
  	 * Remove contaminants and reverse proteins
  	 * Log transformation
  	 * Separate different experiment

### 3. Imputating zero intensities
  	 * Divinding proteins to four groups:
  	 	1. Proteins with zero in all replicas (cases or controls)
  	 	2. Proteins with only one non-zero replica in both cases and controls (only two values available)
  	 	3. Proteins with only one non-zero value in either cases or controls (At least three available values)
  	 		i. abs(log2fold(cases/controls)) < 1
  	 		ii. log2fold(two missing side/the other side) > 1 
  	 		iii. log2fold(two missing side/the other side) < 1
  	 	4. Proteins with two non-zero value in either cases or controls (Up to two zero value in total cases and controls)
  	 	5. Proteins with all values non-zero
  	 * Imputing zero values depends on their groups:
  	 	* Imputing methods description:
  	 		a. method 1 - small intensitye based on the mean and standard value of that condition(specific case or control) across all observed proteins
  	 		b. method 2	- relatively average intensity based on the  
  	 	* group 1, 2, 3.iii are imputed using method 1
  	 	* group 3.i, 3.ii, 4 are imputed using method 2
  	 * Log transformation
  	 * Separate different experiment

###	4. Variance analysis test
	* Filter proteins with less than two peptides count
	* Perform t-test between cases and controls
	* Adjust p-values with Benjamin-Hochberg correction test
	* Calculate log2fold change
	* Select significant proteins with p.adj <0.05 & logfol d> 1
###	5. Normalize and integrate data