# Line 1 CRC IP MS
##### contact moghbaie@rockefeller.edu for any questions
##### 12/10/2018


We had 8 groups as follow:

|	Site of Tumor	|	Condition	|	Tumor	| 	Diagnosis	| Comment	|
| ------------- | ------------- |------------- | ------------- |------------- |
| Ovary  | Igg  |	A	|	Krukenberg Carcinoma	|	Very very strong L1 expresse	|
| Ovary  | Normal  |	A	|	Krukenberg Carcinoma	|	Very very strong L1 expresse	|
|	Liver	|	Igg	|	B	|	Metastatic Rectal Adenocarcinoma	|	Strong 3/3 in liver and node	|
|	Liver	|	Normal	|	B	|	Metastatic Rectal Adenocarcinoma	|	Strong 3/3 in liver and node	|
|	Liver	|	Tumor	|	B	|	Metastatic Rectal Adenocarcinoma	|	Strong 3/3 in liver and node	|
|	Colon	|	Igg	|	C	|	Adenocarcinoma	|	2/3 expression	|
|	Colon	|	Normal	|	C	|	Adenocarcinoma	|	2/3 expression	|
|	Colon	|	Tumor	|	C	|	Adenocarcinoma	|	2/3 expression	|




Project code and data is organized in following order:

* [src](https://github.com/moghbaie/L1_CRC_v2/tree/master/src)
* [Input_data](https://github.com/moghbaie/L1_CRC_v2/tree/master/Input_data/ReadMe_Input.md)
* [Image](https://github.com/moghbaie/L1_CRC_v2/tree/master/Image)


## Project Pipeline:
<img src="https://github.com/moghbaie/L1_CRC_IP_MS/blob/master/CRC_pipeline.png" alt="Pipeline" width="1000"></img>


## QC:

## Imputation:
### Comparison before - after imputation:
<img src="Image/Imputation/comparison_before_after_imputation_Colon_6_2.png" alt="Pipeline" width="900"></img>
<br>
<img src="Image/Imputation/comparison_before_after_imputation_Liver_6_2.png" alt="Pipeline" width="900"></img>
<br>
<img src="Image/Imputation/comparison_before_after_imputation_Ovary_6_2.png" alt="Pipeline" width="600"></img>

### Average intensities before - after imputation:
![Average intensities before imputation](Image/Imputation/average_intensities_before_imputation.png)
![Average intensities after imputation](Image/Imputation/average_intensities_after_imputation.png)

## ANOVA result: 

### Venndiagram
![VennDiagram of significant proteins ](Image/Volcano_plot/Significant_VennDiagram.png)

### Volcano plots and mutations
#### comparison between Colon Tumor tissue against ORF1 and IgG
![Volcano plot compare Colon Tumor tissue to Igg ](Image/Volcano_plot/Volcano_plot_Colon_Tumor_163T_IgG_7.png)
<br>
#### comparison between Colon Tumor tissue and normal tissue against ORF1 
![Volcano plot compare Colon Tumor tissue to Normal ](Image/Volcano_plot/Volcano_plot_Colon_Tumor_163N_ORF1_8.png)
![Observed significant mutation colon tissue to normal ](Image/Phylogenetic_tree/Colon_Tumor_163N_ORF1_8.png)
<br>
#### comparison between Liver Tumor tissue against ORF1 and IgG 
![Volcano plot compare Liver Tumor tissue to Igg ](Image/Volcano_plot/Volcano_plot_Liver_Tumor_159T_IgG_4.png)
![Observed significant Liver colon tissue against ORF1 and IgG ](Image/Phylogenetic_tree/Liver_Tumor_159T_IgG_4.png)
<br>
#### comparison between Liver Tumor tissue and normal tissue against ORF1
![Volcano plot compare Liver Tumor tissue to Normal ](Image/Volcano_plot/Volcano_plot_Liver_Tumor_159N_ORF1_5.png)
![Observed significant Liver Tumor tissue and normal tissue against ORF1 ](Image/Phylogenetic_tree/Liver_Tumor_159N_ORF1_5.png)
<br>
#### comparison between Ovary Tumor tissue against ORF1 and IgG 
![Volcano plot compare Ovary Tumor tissue to Igg - phase1](Image/Volcano_plot/Volcano_plot_Ovary_Tumor_144T_IgG_2.png)
![Observed significant Ovary Tumor tissue against ORF1 and IgG - phase1](Image/Phylogenetic_tree/Ovary_Tumor_144T_IgG_2.png)
<br>
#### comparison between Ovary Tumor tissue and normal tissue against ORF1 
![Volcano plot compare Ovary Tumor tissue to Igg - phase2](Image/Volcano_plot/Volcano_plot_Ovary_Tumor_144T_IgG_10.png)
![Observed significant Ovary Tumor tissue against ORF1 and IgG - phase2](Image/Phylogenetic_tree/Ovary_Tumor_144T_IgG_10.png)

## Integrated plot:
### Venndiagram:
![VenDiagram of expressed proteins ](Image/Integrated_plot/VennDiagram.png)

### MDS plot:
![MDS plot ](Image/Integrated_plot/movie.gif)
<br>
<br>
### Polar map of significant proteins in Colon cancer in different tissue:
![Polar map - comparison of Colon significant proteins in different tissue](Image/Integrated_plot/heatmap_colon_significant.png)

### Heat map of all significant proteins across different tissues:
![Heatmap - all significant proteins](Image/Integrated_plot/heatmap_all_significant.png)