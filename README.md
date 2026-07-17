# Gap Selection Index (GSI)

The Gap Index (GSI) identifies the areas in the country with missing information (species records) and, therefore, the sites where additional sampling will improve biodiversity knowledge. This analysis follows the proposal by [García Márquez et al., 2012](http://www.biodiversity-plants.de/biodivers_ecol/article_meta.php?DOI=10.7809/b-e.00057) and modifications for the environmental dimension by [Aguiar et al., 2020](https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.13137) to identify the spatial coverage of biological information in databases based on density analysis, the climatic representativeness of these records, as well as the taxonomical complementarity of the species in the records. The GSI is quantified using values ​​that range between 0 and 1, with 0 being a well-represented sector and 0 being the underrepresented areas or areas with higher values of information gaps.

## Prerequisites

The index is calculated using the records of species present both in data portals ([SiB Colombia](https://sibcolombia.net/), [SpeciesLink](http://splink.cria.org.br/), [eBird](https://ebird.org/home)) and the information that the Humboldt Institute has compiled in recent years ([Ceiba](http://i2d.humboldt.org.co/ceiba/)). The GSI represents three dimensions as follows: i) quantification of the biological records per square kilometer, ii) the environmental representativeness of each occurrence following the methodology proposed by [Aguiar et al., 2020](https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.13137) and iii) estimated complementary of species richness based on the first-order Jackknife non-parametric estimator. 

For this reason, you need the scripts **1_Record_dimension**, **2_Ambiental_dimension** and **3_Complementarity dimension** to obtain the three GSI´s components, and then, you must run the the **4_GSI** script to obtain the GSI raster. Also, the **GAPfunctions** file is necessary because it contains some of the functions used in the analysis.


### Base Data

1) Records: A dataframe with records of the species and geographical coordinates. The file structure requires the following names in the columns: ID (or gbifID if data downloaded directly from GBIF), species, lat (or decimalLatitude), and lon (or decimalLongitude).

2) Study area: a shapefile of the Area of Insterest (AoI) in .shp format.

### Dependencies

To obtain the results you require.

* [R](https://cran.r-project.org/mirrors.html)
* [RStudio](https://www.rstudio.com/products/rstudio/download/#download)

### Libraries

1. _Record dimension_
   - sf version 1.0-2 
   - maptools version 1.1-1
   - spatstat version 2.2-0
   - raster version 3.4-13

2. _Ambiental dimension_
   - raster version 3.4-13
   - sf version 1.0-2
   - fmsb version 0.7.1
   - dismo version 1.3-3
   - rgdal version 1.5.23
   - MASS version 7.3.53
   - ROCR version 1.0-11
   - rgeos version 0.5-5

3. _Complementarity dimension_
   - raster version 3.4-13
   - rgdal version 1.5.23
   - janitor version 2.1.0
   - dplyr version 1.0.7

4. _GSI_
   - raster version 3.4-13
   - rgdal version 1.5.23
   - hyperSpec version 0.100.0


## How to run

We suggest running the routines step by step, following the order of each script. Nevertheless, you can obtain the result for each dimension independently. 

The database must be stored in a root folder to be read throughout the process.

## Repository structure

In general, this analysis is based on the construction of three independent layers following [García Márquez et al., 2012](http://www.biodiversity-plants.de/biodivers_ecol/article_meta.php?DOI=10.7809/b-e.00057), with modifications for the environmental dimension proposed by [Aguiar et al., 2020](https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.13137). These three dimensions are later integrated into a single final layer representing the Geographic Survey Index (GSI).

Each script in the repository generates one component of the analysis as follows:

| Script | Description |
|---|---|
| `1_Record_dimension.R` | Generates a layer representing the density of biological records as an indicator of sampling effort concentration across the study region. |
| `2_Environmental_dimension.R` | Implements the methodology proposed by Aguiar et al. (2020) to identify regions with environmental conditions that have been poorly sampled compared to the rest of the environmental space. |
| `3_Complementarity_dimension.R` | Builds the complementarity layer, which quantifies sampling completeness for each pixel based on the relationship between observed records and estimated species richness using two non-parametric estimators. |
| `4_GSI.R` | Integrates the three previous layers into a single Geographic Survey Index (GSI) metric. |
| `GAPfunctions.R` | Contains the auxiliary functions used throughout the analysis pipeline. This file is loaded at the beginning of each script. |

![Workflow](Figures/GSI_workflow.png)

 ## Authors and contact

* **[Elkin Alexander Tenorio Moreno](https://github.com/Elkin01)** - *Investigador Adjunto I.Humboldt* -  [Contact](etenorio@humboldt.org.co)
* **[Cristian Alexander Cruz-Rodríguez](https://github.com/crcruzr)** - *Investigador Asistente I.Humboldt* -  [Contact](ccruz@humboldt.org.co)
* **[Elkin Alexi Noguera Urbano](https://github.com/elkalexno)** - *Investigador Titular I.Humboldt* - [Contact](enoguera@humboldt.org.co)
* **Iván gonzález**
* **Laura Carolina Bello**
* **Maria Cecilia Londoño** - *Investiadora Titular I.Humboldt*  

## License

This project is licensed under the MIT License. For details, see the [LICENSE](https://github.com/PEM-Humboldt/gsi_analysis/blob/sf/terra-update/README.md) file.

## Citation

For citation please use the following DOI: https://doi.org/10.5281/zenodo.17228342

## Final considerations

This product contributes to the Annual Operational Plan to the [Instituto Humboldt](http://www.humboldt.org.co/es/) for the year 2021. Specifically to the activity associated with generating a repository with the codes used for the standardization of processes for raising baselines and monitoring biodiversity.

## References

Aguiar, L. M., Pereira, M. J. R., Zortéa, M., & Machado, R. B. (2020). Where are the bats? An environmental complementarity analysis in a megadiverse country. Diversity and Distributions, 26(11), 1510-1522.

Márquez, J. R. G., Dormann, C. F., Sommer, J. H., Schmidt, M., Thiombiano, A., Da, S. S., ... & Barthlott, W. (2012). A methodological framework to quantify the spatial quality of biological databases. Biodiversity and Ecology, 4, 25-39.



