Content:

1) Biodiversity raw original field collection data from Interface Project
2) Dung beetle species classification into habitat association
3) Landscape level forest cover percent at 3 and 5km scale

#########
1) Biodiversity raw original field collection data from Interface Project
File name "beetles_biodiversity.csv"

Contact	Elizabeth Nichols
	USP Sao Paulo
	lizsnichols@gmail.com

Collectors: Elizabeth Nichols, Viviana Alarcon, Gabriel Rocha Oliveira, Bruna Chites
Dates: December 1 2014 - February 28 2015
Financing: FAPESP #  13/23457-6, FAPESP # 14/11676-8, NSF # 1158817
	
Data structure:	
- Landscape: numeric codes for each of 12 landscapes.
- Point: Numeric codes for each of 8 sampling points with each landscape
- Trap_num: Pitfall trap count/identification from 01 to 10
- Year: year of collection
- Month: month of collection, roman
- Day: day of collection, julian
- Municipality: municipality of collection
- Lat: Latitude in decimal degrees
- Long: Longitude in decimal degrees
- Alt: altitude in meters
Columns from 11 to 67 contain species abundances by species latin binomial.
	
	
#########
2) Dung beetle species classification into habitat association
File name: "2022nov_spp_trait_revFB.csv"
Dataset containing trait information of each species to classify into habitat association based on literature data and personal  data from  an expert taxonomist (Prof. Dr. Fernando Augusto B. Silva, coauthor)

Contact	Fernando Augusto B. Silva 
	UFRP, Universidade Federal Rural de Pernambuco
	fernandoabsilva@yahoo.com.br

Data structure:
- Reference: reference study where trait information of species can be found, including personal records from the taxonomist
- Genus
- spp
- latin_bin: latin binomial fo species, concatenated genus and species
- Code: latin binomial without space
- biogeo2: biogeographic association in which species has been collected: AF, Atlantic Forest; WIDE, widely distributed; TRANS, Atlantic Forest transition to Cerrado; CE, Cerrado.
- hab2: classification on habitat of preference: "F", forest; "O", open area specialist; and, "G", generalist
- Classification: resulting two categories of habitat association defined by combining information of the biogeographical distribution and habitat preference: forest specialists (FS), species of range restricted to native forests within the Atlantic Forest biome (= AF + F); and, non-forest specialists (NFS), species with a rather wider range of occurrence, not restricted to forests and/or to the Atlantic Forest region.


#########
3) Landscape level forest cover percent at 3 and 5km scale
File name: "perc_fc_JB.csv"
Data on forest cover percentage of the 12 focal landscapes measured at a 3 and 5 km radius measured based in MapBiomas maps in June 2019. We used the Sampling Design tool in ArcGIS 10.1 to calculate the percentage of native forest cover in circular buffers of 3km and 5km radii around the landscape’s centroids.

Contact	Julia Rodrigues Barreto
	USP Sao Paulo
	barretoj@usp.br

Data structure:
- Landscape: numeric codes for each of 12 landscapes.
- area_forest: forest cover area, we defined as forest those remnants at an intermediate (ca 10 years) or advanced successional stages.
- sum_classes: sum of land cover classes within that landscape
- perfc_3km or perfc_5km: forest cover percent calculated from the two entries above
Columns 5 to 7 contain the same columns, but for 5km scale.
- x and y are geographical coordinates of landscape centroids