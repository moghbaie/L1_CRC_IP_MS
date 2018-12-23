# Line1 CRC V2
##### contact moghbaie@rockefeller.edu for any questions
##### 12/10/2018


We had 8 groups as follow:

|	Site of Tumor	|	Condition	|	Collection #	| Age	|	Sex	|	Diagnosis	| Comment	|
| ------------- | ------------- |------------- | ------------- |------------- | ------------- |------------- |
| Ovary  | Igg  |	144	|	56	|	F 	|	Krukenberg Carcinoma	|	Very very strong L1 expresse	|
| Ovary  | Normal  |	144	|	56	|	F 	|	Krukenberg Carcinoma	|	Very very strong L1 expresse	|
|	Liver	|	Igg	|	159	|	35	|	M 	|	Metastatic Rectal Adenocarcinoma	|	Strong 3/3 in liver and node	|
|	Liver	|	Normal	|	159	|	35	|	M 	|	Metastatic Rectal Adenocarcinoma	|	Strong 3/3 in liver and node	|
|	Liver	|	Tumor	|	159	|	35	|	M 	|	Metastatic Rectal Adenocarcinoma	|	Strong 3/3 in liver and node	|
|	Colon	|	Igg	|	163	|	75	|	F 	|	Adenocarcinoma	|	2/3 expression	|
|	Colon	|	Normal	|	163	|	75	|	F 	|	Adenocarcinoma	|	2/3 expression	|
|	Colon	|	Tumor	|	163	|	75	|	F 	|	Adenocarcinoma	|	2/3 expression	|




This project is done following previous project(Line1 CRC V1)
Detail report exists in following link:
https://docs.google.com/document/d/1Rdn2QTG43JNEqdcNk980U1I6ZvwONEP2scBosMZBJk0/edit

Project code and data is organized in following order:

* [Code](https://github.com/moghbaie/L1_CRC_v2/tree/master/Code/ReadMe_Code.md)
* [Input_data](https://github.com/moghbaie/L1_CRC_v2/tree/master/Input_data/ReadMe_Input.md)
* [Result](https://github.com/moghbaie/L1_CRC_v2/tree/master/Result)
* [Image](https://github.com/moghbaie/L1_CRC_v2/tree/master/Image)


## Project Pipeline:
<img src="https://github.com/moghbaie/L1_CRC_v2/blob/master/CRC_pipeline.png" alt="Pipeline" width="1000"></img>


## QC:

## Imputation:
### Comparison before - after imputation:
![Colon tissue - before and after imputation with method 7_1](Image/Imputation/comparison_before_after_imputation_Colon_7_1.png|width=600)
![Liver tissue - before and after imputation with method 7_1](Image/Imputation/comparison_before_after_imputation_Liver_7_1.png|width=600)
![Ovary tissue - before and after imputation with method 7_1](Image/Imputation/comparison_before_after_imputation_Ovary_7_1.png|width=400)

### Average intensities before - after imputation:
![Average intensities before imputation](Image/Imputation/average_intensities_before_imputation.png)
![Average intensities after imputation](Image/Imputation/average_intensities_after_imputation.png)

## ANOVA result: 

### Venndiagram
![VennDiagram of significant proteins ](Image/Volcano_plot/Significant_VennDiagram.png)

### Volcano plots
![Volcano plot compare Colon Tumor tissue to Igg ](Image/Volcano_plot/Volcano_plot_Colon_Tumor_Igg.png)
![Volcano plot compare Colon Tumor tissue to Normal ](Image/Volcano_plot/Volcano_plot_Colon_Tumor_Normal.png)
![Volcano plot compare Liver Tumor tissue to Igg ](Image/Volcano_plot/Volcano_plot_Liver_Tumor_Igg.png)
![Volcano plot compare Liver Tumor tissue to Normal ](Image/Volcano_plot/Volcano_plot_Liver_Tumor_Normal.png)
![Volcano plot compare Ovary Tumor tissue to Igg ](Image/Volcano_plot/Volcano_plot_Ovary_Tumor_Igg.png)

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