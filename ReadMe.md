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
<img src="https://github.com/moghbaie/L1_CRC_v2/blob/master/CRC_pipeline.png" alt="Vega Examples" width="1000"></img>


## QC

## Imputation
![Colon tissue - before and after imputation with method 7_1](https://github.com/moghbaie/L1_CRC_v2/tree/master/Image/Imputation/comparison_before_after_imputation_Colon_7_1.png)
<img src="https://github.com/moghbaie/L1_CRC_v2/tree/master/Image/Imputation/comparison_before_after_imputation_Colon_7_1.png" alt="Vega Examples" width="800"></img>
<img src="https://github.com/moghbaie/L1_CRC_v2/tree/master/Image/Imputation/comparison_before_after_imputation_Liver_7_1.png" alt="Vega Examples" width="800"></img>
<img src="https://github.com/moghbaie/L1_CRC_v2/tree/master/Image/Imputation/comparison_before_after_imputation_Ovary_7_1.png" alt="Vega Examples" width="800"></img>
<img src="https://github.com/moghbaie/L1_CRC_v2/tree/master/Image/Imputation/average_intensities_before_imputation.png" alt="Vega Examples" width="800"></img>
<img src="https://github.com/moghbaie/L1_CRC_v2/tree/master/Image/Imputation/average_intensities_after_imputation.png" alt="Vega Examples" width="800"></img>

## Anova 
### Venndiagram
<img src="https://github.com/moghbaie/L1_CRC_v2/tree/master/Image/Volcano_plot/Significant_VennDiagram.png" alt="Vega Examples" width="1000"></img>

### Volcano plots
<img src="https://github.com/moghbaie/L1_CRC_v2/tree/master/Image/Volcano_plot/Volcano_plot_Colon_Tumor_Igg.png" alt="Vega Examples" width="800"></img>
<img src="https://github.com/moghbaie/L1_CRC_v2/tree/master/Image/Volcano_plot/Volcano_plot_Colon_Tumor_Normal.png" alt="Vega Examples" width="800"></img>
<img src="https://github.com/moghbaie/L1_CRC_v2/tree/master/Image/Volcano_plot/Volcano_plot_Liver_Tumor_Igg.png" alt="Vega Examples" width="800"></img>
<img src="https://github.com/moghbaie/L1_CRC_v2/tree/master/Image/Volcano_plot/Volcano_plot_Liver_Tumor_Normal.png" alt="Vega Examples" width="800"></img>
<img src="https://github.com/moghbaie/L1_CRC_v2/tree/master/Image/Volcano_plot/Volcano_plot_Ovary_Tumor_Igg.png" alt="Vega Examples" width="800"></img>

## Integrated plot
### Venndiagram
<img src="https://github.com/moghbaie/L1_CRC_v2/tree/master/Image/Integrated_plot/VennDiagram.png" alt="Vega Examples" width="1000"></img>

### MDS plot
<img src="https://github.com/moghbaie/L1_CRC_v2/tree/master/Image/Integrated_plot/movie.gif" alt="Vega Examples" width="1000"></img>

### Polar map of significant proteins in Colon cancer in different tissue
<img src="https://github.com/moghbaie/L1_CRC_v2/tree/master/Image/Integrated_plot/heatmap_colon_significant.png" alt="Vega Examples" width="1000"></img>

### Heat map of all significant proteins across different tissues
<img src="https://github.com/moghbaie/L1_CRC_v2/tree/master/Image/Integrated_plot/heatmap_all significant.png" alt="Vega Examples" width="400"></img>