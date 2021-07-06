# alMIghTO
A curated database of genomic features of all the available and suitable complete mitochondrial genomes from GenBank.

# allMIghTO R functions

Jointly with the dataset allMIghTO, we provide some interactive R functions, which are tailored for the dataset and they allow to easily retrive different plots (box-plots, correlograms, PCAs). These functions require the following R libraries:
* [ggplot2](https://ggplot2.tidyverse.org/)
* [ggsignif](https://const-ae.github.io/ggsignif/)
* [RColorBrewer](https://www.r-graph-gallery.com/38-rcolorbrewers-palettes.html)
* [dunn.test](https://cran.r-project.org/web/packages/dunn.test/dunn.test.pdf)
* [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html)
* [Hmisc](https://hbiostat.org/R/Hmisc/)
* [corrplot](https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html)
* [ade4](https://cran.r-project.org/web/packages/ade4/index.html)

## Loading functions and dataset

In order to load the functions on your R environment, use the comand:
```
source("Statics_GH.R")
```
Load the dataset on R, you can easily load the dataset in a variable using the function `dataset(path)`. Example:
```
mito=dataset("alMIghTO-v1.0.csv")
```
Now it will briefly explained how to build each plot:

## Box-plot

To create a boxplot, the function is `boxplot(column,dataset,list of groups,significance,return)`.
An example of the function:
```
boxplot(3,mito,list(Arthropoda="Arthropoda",Vertebrata=c("Mammalia","Aves"),Nematoda="Nematoda",Platyhelmithes="Platyhelminthes"),TRUE,FALSE)
```
[See output example](https://github.com/AlessandroFormaggioni/stunning-fiesta/blob/main/images_gh/Boxplot1_gh.pdf)

1. The first argument is the number of the **column** where the statistic is placed. Downloading the dataset, the statistics will be ordered as follows:

| Static | Column |
| ------ | ------ |
| URs | 3 |
| SUskew | 4 |
| AT | 5 |
| ATskew | 6 |
| GCskew | 7 |
| Genes | 8 |
| Length | 9 |
| CAI | 10 |
| URs_AT | 12 |
| URs\_Med\_Len | 13 |

2. The variable that contains the **dataset**

3. A **list of clades**, each component in the list can be a string with the name of a clade or a vector of clade, if you want to merge different clades in a single boxplot (as for the "Vertebrata" group in the example) 
4. **TRUE** to display the pairwise significance between groups, **FALSE** to not display it. The significance is calculated trough the Kruskal-Wallis test, followed by the Dunn's test with Bonferroni correction for pairwise comparisons. 

5. **FALSE** to display the plot and decide whether to save it or not. With **TRUE** the function returns a ggplot variable, useful if you want to concatenate more plots. Example: 
```
p3=boxplot(3,mito,list(Arthropoda="Arthropoda",Vertebrata=c("Mammalia","Aves"),Nematoda="Nematoda",Platyhelmithes="Platyhelminthes"),FALSE,TRUE)
p4=boxplot(5,mito,list(Arthropoda="Arthropoda",Vertebrata=c("Mammalia","Aves"),Nematoda="Nematoda",Platyhelmithes="Platyhelminthes"),FALSE,TRUE)
plot=grid.arrange(p3,p4, ncol=2, nrow=1)
ggsave("Boxplots.pdf",plot=plot)
```
[See output example](https://github.com/AlessandroFormaggioni/stunning-fiesta/blob/main/images_gh/Boxplot2_gh.pdf)

## Correlogram
A correlogram graphically reports all the pairwise correlations between variables (in our case the 11 statistics). The function calculates the Spearman's ρ between all comparisons. The p-value treshold is corrected with the Bonferroni method. The function is `corr(dataset,clade)`
Example:
```
corr(mito,"Metazoa")
```
[See output example](https://github.com/AlessandroFormaggioni/stunning-fiesta/blob/main/images_gh/Corrplot_gh.pdf) 

The arguments of the function are the following ones:
1. The variable that contains the **dataset**

2. With the **clade** we indicate the subset on which the correlations will be calculated

## PCA

To create a PCA, the function is `pca(dataset,clade,sub rank, number of ellipses)`. The data are scaled and centered. Just the first 4 PCs can be displayed. 
Example:
```
pca(mito,"Metazoa","phylum",8)
```
[See output example](https://github.com/AlessandroFormaggioni/stunning-fiesta/blob/main/images_gh/PCA1_gh.pdf) <br />
The arguments of the function are the following ones:
1. The whole **dataset**

2. With the **clade** we indicate the subset on which the PCs will be calculated

3. Indicating a **sub rank** we determine how the samples will be clustered, in the example samples are clustered according to the phyla they belong

4. The **number of ellipses** we want to display, deciding to calculate and display the ellipses of the n most numerous clusters. N.B. the ellipses of clusters with few points could fail to be calculated, if the ellipses are not displayed try to reduce the number of them. 
Inside the function it is also asked if you want to display the correlation circle, if so, you have to provide 4 coordiantes in order to place the correlation circle on the PCA plot (Help yourself with the two axis of the PCA plot to decide the 4 coordinates and find a blank spot). <br />

[See output example](https://github.com/AlessandroFormaggioni/stunning-fiesta/blob/main/images_gh/PCA2_gh.pdf)

These functions do not allow to exploit all the ecological information present in the database but they will be added in the future. If you need to obtain these plots divided by ecological features, feel free to contact the author (ale.formaggioni@gmail.com), as well as for any question, concern or advice.
