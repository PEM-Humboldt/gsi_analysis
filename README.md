# Gap-selection-Index--GSI

The Gap Index (GSI) allows identifying the areas in the country with missing information and therefore the sites in which additional sampling will improve knowledge on biodiversity. This approach is carried out following what was proposed by [García Márquez et al., 2012](http://www.biodiversity-plants.de/biodivers_ecol/article_meta.php?DOI=10.7809/b-e.00057) to identify the spatial coverage of biological information in databases, the environmental representativeness of these records, as well as the complementarity of the species in the records. The GSI is quantified using values ​​that range between 0 and 1, with 0 being a well-represented sector, and 0 being the underrepresented areas or areas with greater information gaps.

## Prerequisites

To obtain it, the index is calculated using the records of species present both in data portals and the information that the Humboldt Institute has compiled in recent years. To calculate it, three dimensions must be estimated: i) a quantification of the biological records per square kilometer, ii) the environmental representativeness of each occurrence following the methodology proposed by [Aguiar et al., 2020](https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.13137) and iii) The identification of the complementarity of the species that have records per square kilometer based on first-order Jackknife nonparametric estimates. 

For this reason you need the scripts **1_Record_dimension**, **2_Ambiental_dimension** and **3_Complementarity dimesion** to obtain the three elemenatry elements and the **4_GSI** script to obtain the GSI raster. Also it's necessary the **GAPfunctions** file, because it contains some of the functions used in the analysis.


### Base Data

Data frame with records of the species and geographical coordinates. The structure of the file requies the following names in the columns: _ID_, _species_, _lat_ and _lon_.

### Dependencies

To obtain the results you requires 

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

We suggest that run the rutines following the numbers in the scrips. Nevertheless, you can obtain the result for each dimenson independently. 

It is important that the base data will be used throughout the process, and that this data is stored in the root of the folder.

 ## Authors and contact
 
* **[Cristian Alexander Cruz-Rodríguez](https://github.com/crcruzr)** - *Investigador Asistente I.Humboldt* -  [Contact](ccruz@humboldt.org.co)
* **[Elkin Alexi Noguera Urbano](https://github.com/elkalexno)** - *Investiador Adjunto I.Humboldt* - [Contact](enoguera@humboldt.org.co)
* **Maria Cecilia Londoño** - *Investiadora Titular I.Humboldt* 
* **Iván gonzález**
* **Laura Carolina Bello**
 

## License

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/PEM-Humboldt/indicadores-sinap/blob/master/LICENSE.md) file for details



## Final considerations

This product contributes to the Annual Operational Plan to the [Instituto Humboldt](http://www.humboldt.org.co/es/) for the year 2021. Specifically to the activity associated with generating a repository with the codes used for the standardization of processes for raising baselines and monitoring biodiversity.
