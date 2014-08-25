
metagene: A package to produce Metafeature plots
========================================================


This repository contains R functions used to do multiple ChIP-Seq experiments comparisons.

This package produces Metagene-like plots to compare the behavior of DNA-interacting proteins at selected groups of features. A typical analysis can be done in viscinity of transcription start sites (TSS) of genes or at any regions of interest (such as enhancers). Multiple combinations of group of features and/or group of bam files can be compared in a single analysis. Bootstraping analysis is used to compare the groups and locate regions with statistically different enrichment profiles. In order to increase the sensitivity of the analysis, alignment data is used instead of peaks produced with peak callers (i.e.: MACS2 or PICS). The metagene package uses bootstrap to obtain a better estimation of the mean enrichment and the confidence interval for every group of samples.

Currently supported species are **human** and **mouse**.

## Authors ##

[Charles Joly Beauparlant](http://ca.linkedin.com/pub/charles-joly-beauparlant/89/491/3b3 "Charles Joly Beauparlant"), [Fabien Claude Lamaze](http://ca.linkedin.com/in/fabienlamaze/en "Fabien Claude Lamaze"), [Rawane Samb](http://ca.linkedin.com/in/rawanesamb "Rawane Samb"), [Astrid Louise Deschenes](http://ca.linkedin.com/in/astriddeschenes "Astrid Louise Deschenes") and [Arnaud Droit](http://ca.linkedin.com/in/drarnaud "Arnaud Droit").

See [Arnaud Droit Lab](http://bioinformatique.ulaval.ca/home/ "Arnaud Droit Lab") website.

## License ##

This package and the underlying metagene code are distributed under the Artistic license 2.0. You are free to use and redistribute this software. 

For more information on Artistic 2.0 License see [http://opensource.org/licenses/Artistic-2.0](http://opensource.org/licenses/Artistic-2.0)

## Bugs/Feature requests ##

If you have any bugs or feature requests, [let us know](https://github.com/CharlesJB/metagene/issues). Thanks!
